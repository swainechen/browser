#!/usr/bin/perl -w
#
# need a full set of files - look for these based on read
# .srst2.gz - goes into Resistance and Genes
# .gcov.gz - some info for Fastq
# .lofreq.gz - not really needed here
# .tgz - info for Fastq and Assembly
#
use warnings;
use strict;
use slchen;
use Archive::Tar;
use DBI;
use DateTime::Format::DBI;
use DateTime::Format::DateParse;
use Data::Dumper;
use Compress::Zlib;
use File::Spec;
use Getopt::Long;
use GERMS;

my $GET_FILES = `which get_files.pl`; chomp $GET_FILES;
my $RUN_KRAKEN = `which run-kraken.pl`; chomp $RUN_KRAKEN;
#my $SRADB = $ENV{"SRADB"};
my %ATTR;	# needed for SQLite connection - ?
my $USE_DB = 0;
my $verbose = 0;
my $force = 0;	# update no matter what - ignore date stamps
&Getopt::Long::Configure("pass_through");
GetOptions (
  'db!' => \$USE_DB,
  'verbose!' => \$verbose,
  'force!' => \$force
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
my $fastq_data;
my $mlst_data;
my $resistance_data;
my $genes_data;
my $status;
my @files;
my @ft;
my $kraken_command;
my $kraken_out;
my $gz;
my $missing;
my @data;
# inclusion of /mnt/software libraries is messing up local workstation awk
my $ENV_save = $ENV{LD_LIBRARY_PATH};

sub print_usage {
  print "Usage: $0 <Run ID> [ -db ] [ -verbose ] [ -force|noforce ]\n";
  print "  With -db do the database import / update\n";
  print "  With -force then always update (Fastq table)\n";
  print "  With -verbose print out all the information collected\n";
}

if (!defined $ARGV[0] || !length($ARGV[0])) {
  &print_usage;
  exit;
}
# check for files
my $runID = $ARGV[0];
my $DBH = GERMS::dbconnect();
my $DB_PARSER = DateTime::Format::DBI->new($DBH);
my $TIP = GERMS::get_browser_tip($runID, $DBH, $USE_DB);
if ($verbose) {
  print "RunID: $runID\n";
  print "TIP: $TIP\n";
  if ($USE_DB) {
    print "Doing database operations (-db is on)\n";
  } else {
    print "No database operations (no -db flag)\n";
  }
}
# just check if things are complete
if ($TIP && !$force) {
  $missing = 0;
  $sql = "SELECT * FROM Fastq JOIN Tips ON Tips.TIP = Fastq.TIP WHERE Tips.TIP = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($TIP);
  while (@data = $sth->fetchrow_array()) {
    foreach $i (0..$#data) {
      if (!defined $data[$i] || $data[$i] eq "") {
        $missing++;
      }
    }
    last;
  }
  $missing = 1 if $sth->rows <= 0;
  if (!$missing) {
    if ($verbose) {
      print "Information already present for $runID ($TIP):\n";
      foreach $i (0..$#data) {
        print "  $data[$i]\n";
      }
    }
    exit;
  }
}

# Fastq
# use get_files.pl to get files
# Paired, NumReads, ReadLength, Kraken from there
# NumReads from the first file
# ReadLength from the first read
# For GenBank SRA
#   Technology is platform in the sra database
#   Machine is instrument_model in the sra database
# For GIS Libraries
#   Technology and Machine should come from Studies table
#
$fastq_data->{TIP} = $TIP;
@files = `$GET_FILES -checkonly $runID`;
foreach $i (reverse 0..$#files) {
  chomp $files[$i];
  if (!-f $files[$i] || !-s $files[$i]) {
    splice @files, $i, 1;
  }
}
if ($#files >= 1) {
  $fastq_data->{Paired} = "Yes";
  if ($#files > 1) {
    print STDERR "Warning: ", scalar(@files), " files found for $runID, calling this paired anyway\n";
  }
} elsif ($#files == 0) {
  $fastq_data->{Paired} = "No";
}
die "No files found (or all empty) for $runID\n" if !@files;
if (-f $files[0] && -s $files[0]) {
  if ($verbose) {
    print "Files:\n", join ("\n", @files), "\n";
  }
  @ft = GERMS::file_type($files[0]);
  if (lc($ft[1]) eq "gzip") {
    $ENV{LD_LIBRARY_PATH} = "";
    $fastq_data->{NumReads} = `zcat $files[0] | wc -l | awk '{printf "\%.0f", \$1/4}'`;
    $ENV{LD_LIBRARY_PATH} = $ENV_save;
    chomp $fastq_data->{NumReads};
    # use Compress::Zlib because getting stdout: Broken pipe errors
    $gz = Compress::Zlib::gzopen($files[0], "rb");
    $gz->gzreadline($fastq_data->{ReadLength});
    $gz->gzreadline($fastq_data->{ReadLength});
    $gz->gzclose;
    chomp $fastq_data->{ReadLength};
    $fastq_data->{ReadLength} = length($fastq_data->{ReadLength});
  } else {
    $ENV{LD_LIBRARY_PATH} = "";
    $fastq_data->{NumReads} = `wc -l $files[0] | awk '{printf "\%.0f", \$1/4}'`;
    $ENV{LD_LIBRARY_PATH} = $ENV_save;
    chomp $fastq_data->{NumReads};
    $fastq_data->{ReadLength} = `head -n 2 $files[0] | tail -n 1`;
    chomp $fastq_data->{ReadLength};
    $fastq_data->{ReadLength} = length($fastq_data->{ReadLength});
  }

  # check if there's already Kraken data
  $sql = "SELECT Kraken, KrakenPercent, Kraken2, Kraken2Percent FROM Fastq WHERE TIP = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($fastq_data->{TIP});
  while (@f = $sth->fetchrow_array()) {
    if (defined $f[0] && length($f[0])) {
      $fastq_data->{Kraken} = $f[0];
      print "Found existing Kraken data $f[0] for $runID ($TIP)\n" if $verbose;
    }
    if (defined $f[1] && length($f[1])) {
      $fastq_data->{KrakenPercent} = $f[1];
    }
    if (defined $f[2] && length($f[2])) {
      $fastq_data->{Kraken2} = $f[2];
    }
    if (defined $f[3] && length($f[3])) {
      $fastq_data->{Kraken2Percent} = $f[3];
    }
    last;
  }
  if (!defined $fastq_data->{Kraken} || $force) {
    $kraken_command = "$RUN_KRAKEN -q1 $files[0]";
    if (defined $files[1] && -f $files[1] && -s $files[1]) {
      $kraken_command .= " -q2 $files[1]";
    }
    print "Running $kraken_command\n" if $verbose;
    $kraken_out = `$kraken_command`;
    @f = split /\n/, $kraken_out;
    foreach $i (0..$#f) {
      chomp $f[$i];
      @g = split /\t/, $f[$i];
      if ($i == 0) {
        $fastq_data->{Kraken} = $g[0];
        $fastq_data->{KrakenPercent} = $g[1];
      } elsif ($i == 1) {
        $fastq_data->{Kraken2} = $g[0];
        $fastq_data->{Kraken2Percent} = $g[1];
      }
    }
  }

  # get Machine and Technology data
  if ($runID =~ /^[SED]RR\d+/) {
#    my $sradbh = DBI->connect("DBI:SQLite:dbname=$SRADB", "", "", \%ATTR);
#    my $srasql = "SELECT platform, instrument_model FROM sra WHERE run_accession = ?";
#    my $srasth = $sradbh->prepare($srasql);
#    $srasth->execute($runID);
#    while (@f = $srasth->fetchrow_array()) {
#      $fastq_data->{Technology} = $f[0];
#      $fastq_data->{Technology} =~ s/^\s+//;
#      $fastq_data->{Technology} =~ s/\s+$//;
#      $fastq_data->{Machine} = $f[1];
#      $fastq_data->{Machine} =~ s/^\s+//;
#      $fastq_data->{Machine} =~ s/\s+$//;
#    }
#    $sradbh->disconnect();
     $fastq_data->{Technology} = "";
     $fastq_data->{Machine} = "";
  } else {
    $sql = "SELECT Machine FROM Studies WHERE Run = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($runID);
    while (@f = $sth->fetchrow_array()) {
      if (defined $f[0] && length $f[0]) {
        $fastq_data->{Machine} = $f[0];
        if ($fastq_data->{Machine} =~ /PacBio/i) {
          $fastq_data->{Technology} = "PACBIO";
        } else {
          $fastq_data->{Technology} = "ILLUMINA";
        }
      }
      last;
    }
  }
  if (!defined $fastq_data->{Technology}) {
    $fastq_data->{Technology} = "ILLUMINA";
  }
  if (!defined $fastq_data->{Machine}) {
    $fastq_data->{Machine} = "";
  }

  $status = GERMS::do_db_withSourceFile($fastq_data, "Fastq", $force, $DBH, $USE_DB);
  if ($verbose) {
    foreach $i (sort keys %$fastq_data) {
      print "  $i: $fastq_data->{$i}\n";
    }
    print "Fastq: $status\n";
  }
} # Fastq
