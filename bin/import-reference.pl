#!/usr/bin/perl -w
#
# source file is an fna file
# also need genomes_proks.txt
# check for changes and do updates
# database should be germs_browser
# table should be ReferenceGenomes
#
use warnings;
use strict;
use slchen;
use Archive::Tar;
use DBI;
use DateTime::Format::DBI;
use DateTime::Format::DateParse;
use Data::Dumper;
use File::Spec;
use File::Basename;
use File::Fetch;
use File::Temp;
use File::Copy;
use Getopt::Long;
use GERMS;

my $show_help = 0;
my $USE_DB = 0;
my $verbose = 0;
my $download = 0;
my $species = "";
my $set_default = 0;
my $force = 0;	# update no matter what - ignore date stamps
my $OUTBREAK_BASE = $ENV{'OUTBREAK_BASE'}

&Getopt::Long::Configure("pass_through");
GetOptions (
  'species=s' => \$species,
  'default!' => \$set_default,
  'db!' => \$USE_DB,
  'verbose!' => \$verbose,
  'force!' => \$force,
  'help' => \$show_help,
  'download!' => \$download
);

my $data;
my $parsed_data;
my @f;
my @g;
my $class;
my $data_keys;
my @out;
my $id;
my $i;
my $j;
my $sth;
my $sql;
my @array;
my $param_number;
my $param_map;
my @return;
my $operation;
my $status;
my $accession;
my $linfo;
my $tempdir;
my $ff;
my $ffout;

sub print_usage {
  print "Usage: $0 <reference.fna> -species <species> [ -nodefault|default ][ -db ] [ -verbose ] [ -force ] [ -download ]\n";
  print "  Will import/update the data to the ReferenceGenomes table.\n";
  print "  Species should be like Ecoli, Senterica, Kpneumoniae.\n";
  print "  If -default, make this the default reference genome\n";
  print "  With -db do database import / update (default just show what would be done)\n";
  print "  With -force then always update\n";
  print "  With -verbose print out all the information collected\n";
  print "  With -download go download the .fna and .gff files\n";
}

if ($show_help || $species eq "") {
  &print_usage;
  exit;
}
if (!defined $ARGV[0] || (!$download && !-f $ARGV[0])) {
  &print_usage;
  exit;
}
if (!defined $OUTBREAK_BASE || !-d $OUTBREAK_BASE) {
  &print_usage;
  print "\n-- Error, OUTBREAK_BASE not set or doesn't exist\n";
  exit;
}

$data = ();
$data->{ReferenceFile} = File::Spec->rel2abs($ARGV[0]);
$accession = File::Basename::basename($data->{ReferenceFile});
$accession =~ s/_genomic\.fna$//;
$linfo = parse_lproks($accession);
if ($download && !-f $ARGV[0]) {
  print "Downloading $linfo->{fnaurl}\n" if $verbose;
  $tempdir = File::Temp::tempdir ( CLEANUP => 1 );
  $ff = File::Fetch->new( uri => $linfo->{fnaurl} );
  $ff->fetch( to => $tempdir );
  $ffout = $ff->file;
  if ($ffout =~ /\.gz$/) {
    system("gunzip $tempdir/$ffout");
    $ffout =~ s/\.gz$//;
    File::Copy::move("$tempdir/$ffout", File::Basename::dirname($data->{ReferenceFile}));
  }
}
$data->{GFFFile} = $data->{ReferenceFile};
$data->{GFFFile} =~ s/\.fna$/.gff/;
if ($download && !-f $data->{GFFFile}) {
  print "Downloading $linfo->{gffurl}\n" if $verbose;
  $tempdir = File::Temp::tempdir ( CLEANUP => 1 );
  $ff = File::Fetch->new( uri => $linfo->{gffurl} );
  $ff->fetch( to => $tempdir );
  $ffout = $ff->file;
  if ($ffout =~ /\.gz$/) {
    system("gunzip $tempdir/$ffout");
    $ffout =~ s/\.gz$//;
    File::Copy::move("$tempdir/$ffout", File::Basename::dirname($data->{ReferenceFile}));
  }
}

$data->{Species} = $species;
$data->{ReferenceName} = $accession;
$data->{ReferenceMD5} = `md5sum $data->{ReferenceFile} | cut -f1 -d' '`;
chomp $data->{ReferenceMD5};
$data->{GFFMD5} = `md5sum $data->{GFFFile} | cut -f1 -d' '`;
chomp $data->{GFFMD5};

if (defined $linfo->{shortname} && length($linfo->{shortname})) {
  $data->{ReferenceName} = $linfo->{shortname}
}
$data->{URL} = $linfo->{fnaurl};
$data->{DefaultReference} = $set_default;

open F, $data->{ReferenceFile};
my $sequence = ();
my $name;
while (<F>) {
  chomp;
  next if /^$/;
  if (s/^>//) {
    $name = $_;
    if (!defined $sequence->{$name}) {
      $sequence->{$name} = "";
    } else {
      die "Found multiple fasta headers for $name in $ARGV[0]\n";
    }
  } else {
    $sequence->{$name} .= $_;
  }
}
$data->{TotalLength} = 0;
$data->{NumReplicons} = 0;
$data->{ChromosomeLength} = 0;
foreach $name (keys %$sequence) {
  my $shortname = (split (/\s+/, $name))[0];
  $data->{TotalLength} += length($sequence->{$name});
  $data->{NumReplicons}++;
  if (defined $linfo->{name}->{$name} &&
      $linfo->{type}->{$name} eq "chromosome") {
    $data->{ChromosomeLength} += length($sequence->{$name});
  } elsif (defined $linfo->{name}->{$shortname} &&
           $linfo->{type}->{$shortname} eq "chromosome") {
    $data->{ChromosomeLength} += length($sequence->{$name});
  }
}

my $DBH = GERMS::dbconnect("germs_browser");

# need to manually handle default
if ($data->{DefaultReference} && $USE_DB) {
  if ($verbose) {
    print "Setting this default. Unsetting default for other references for $data->{Species}\n";
  }
  $sql = "UPDATE ReferenceGenomes SET DefaultReference = 0 WHERE Species = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($data->{Species});
}

# mangle the paths to eliminate $OUTBREAK_BASE in front
$data->{ReferenceFile} =~ s/^$OUTBREAK_BASE\///;
$data->{GFFFile} =~ s/^$OUTBREAK_BASE\///;

# do the database checking and inserting/updating
$status->{GERMS::do_db($data, "ReferenceGenomes", $force, $DBH, $USE_DB)}++;

foreach $i ( qw(INSERT UPDATE NOTHING) ) {
  $status->{$i} = 0 if !defined $status->{$i};
}
if ($verbose) {
  print "For file ", File::Basename::basename($ARGV[0]), "\n";
  print "  INSERT for $status->{INSERT} lines\n";
  print "  UPDATE for $status->{UPDATE} lines\n";
  print "  NOTHING for $status->{NOTHING} lines\n";
}

sub parse_lproks {
  my ($accession) = @_;
  my $lproks = "$ENV{GERMS_DATA}/genomes_proks.txt";
  my @f;
  my @g;
  my @h;
  my $i;
  my $j;
  my $r = ();
  my @replicons;
  $r->{urlbase} = "";
  open F, $lproks;
  while (<F>) {
    chomp;
    next if /^#/;
    next if /^$/;
    @f = split /\t/, $_;
    if ($f[18] =~ /$accession/) {
      $r->{urlbase} = $f[18];
    } elsif ($f[19] =~ /$accession/) {
      $r->{urlbase} = $f[19];
    }
    if ($r->{urlbase} ne "") {
      @g = split /\//, $r->{urlbase};
      $r->{fnaurl} = "$r->{urlbase}/$g[$#g]_genomic.fna.gz";
      $r->{gffurl} = "$r->{urlbase}/$g[$#g]_genomic.gff.gz";
      $r->{longname} = $f[0];
      $r->{shortname} = $f[1];
      # $f[10] is for ex:
      # chromosome:NC_002488.3/AE003849.1; plasmid pXF1.3:NC_002489.3/AE003850.3; plasmid pXF51:NC_002490.1/AE003851.1
      @replicons = split /;/, $f[10];
      foreach $i (0..$#replicons) {
        # now each should be like:
        # chromosome:NC_002488.3/AE003849.1
        # plasmid pXF1.3:NC_002489.3/AE003850.3
        $replicons[$i] =~ s/^\s+//;
        $replicons[$i] =~ s/\s+$//;
        @g = split /:/, $replicons[$i];
        @h = split /\//, $g[1];
        foreach $j (0..$#h) {
          $r->{name}->{$h[$j]} = $g[0];
          $r->{type}->{$h[$j]} = $g[0];
          $r->{type}->{$h[$j]} =~ s/\s+.*$//;
        }
      }
      last;
    }
  }
  return $r;
}
