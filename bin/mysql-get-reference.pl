#!/usr/bin/perl -w
#
# Take in a run number
# Print out species, reference sequence
#
use warnings;
use strict;
use File::Spec;
use File::Basename;
use File::Path;
use File::Fetch;
use GERMS;
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $quiet = 1;
my $force_species = "NONE";
my $force_reference = "NONE";
my $OUTBREAK_BASE = $ENV{'OUTBREAK_BASE'};
my $makedir = 1;
GetOptions (
  'quiet!' => \$quiet,
  'species=s' => \$force_species,
  'reference=s' => \$force_reference,
  'makedir!' => \$makedir
);

if (!defined $ARGV[0] || !length($ARGV[0]) || !defined $OUTBREAK_BASE || -f $OUTBREAK_BASE) {
  print "Usage: $0 <Run ID> [ -species <species> ] [ -reference <reference> ] [ -makedir|nomakedir ]\n";
  print "Will print out:\n";
  print "  Directory for rest of files\n";
  print "  Species (abbreviation) for SRST2\n";
  print "  Reference sequence for mapping\n";
  print "Can force Species, use \"Ecoli\" or \"Escherichia coli\" (quoted as needed)\n";
  print "Can force Reference, use Reference filename, Reference md5sum, Reference genome name (like \"EC958\" or \"LT2\")\n";
  print "With -makedir will make the directories and download and index necessary references in OUTBREAK_BASE\n";
  if (!defined $OUTBREAK_BASE || -f $OUTBREAK_BASE) {
    print "-- Error, OUTBREAK_BASE environment variable not set\n";
  }
  exit;
}

my @f;
my @data;
my $directory = "";
my $fullspecies = "";	# this is like "Escherichia coli"
my $species = "";	# this is like "Ecoli"
my $reference = "";
my $URL = "";
my $DBH = GERMS::dbconnect();
my $sql;
my $sth;
my $temp;
my $ff;

if (uc($force_species) eq "NONE") {
  $sql = "SELECT Kraken FROM Fastq JOIN Tips ON Fastq.TIP = Tips.TIP WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($ARGV[0]);
  while (@data = $sth->fetchrow_array()) {
    if (length $data[0]) {
      $fullspecies = $data[0];
      last;
    }
  }
} else {
  $fullspecies = $force_species;
}
if (length $fullspecies) {
  # convert "Escherichia coli" to "Ecoli"
  # should also keep "Ecoli" as "Ecoli"
  if ($fullspecies =~ /^[A-Z]\w+$/) {
    # this should be only if $force_species was like "Ecoli" already
    # try to rescue this from species we know about
    $species = $fullspecies;
    $sql = "SELECT DISTINCT Kraken FROM Fastq";
    $sth = $DBH->prepare($sql);
    $sth->execute();
    while (@data = $sth->fetchrow_array()) {
      next if !defined $data[0];
      @f = split /\s+/, $data[0];
      if (defined $f[1]) {
        $temp = substr($f[0], 0, 1) . $f[1];
      }
      if ($temp eq $species) {
        $fullspecies = $data[0];
        last;
      }
    }
  } else {
    @f = split /\s+/, $fullspecies;
    if (defined $f[1]) {
      $species = substr($f[0], 0, 1) . $f[1];
    }
  }
}

if (uc($force_reference) eq "NONE") {
  $sql = "SELECT ReferenceFile, URL FROM ReferenceGenomes WHERE DefaultReference = 1 AND Species = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($species);
  while (@data = $sth->fetchrow_array()) {
    if (length $data[0]) {
      $reference = $data[0];
      $URL = $data[1];
      last;
    }
  }
} else {
  # this should enable using just the base file name, md5sum,
  # or reference name as $force_reference
  $species = "";
  $sql = "SELECT Species, ReferenceFile, URL from ReferenceGenomes WHERE ReferenceFile LIKE '%?' OR ReferenceMD5 = ? OR ReferenceName = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute("%" . $force_reference, $force_reference, $force_reference);
  while (@data = $sth->fetchrow_array()) {
    $species = $data[0];
    $reference = $data[1];
    $URL = $data[2];
    last;
  }
  if ($species eq "" && -f $force_reference) {
    $reference = File::Spec->rel2abs($force_reference);
  }
}

# check for files, make dirs and download if needed
if (-f $reference) {
  $reference = File::Spec->rel2abs($reference);
  $directory = File::Basename::dirname($reference);
} elsif (-f "$OUTBREAK_BASE/$reference") {
  $reference = File::Spec->rel2abs("$OUTBREAK_BASE/$reference");
  $directory = File::Basename::dirname($reference);
} elsif ($reference ne "" && $reference ne "NONE") {
  $reference =~ s/^$OUTBREAK_BASE\///;
  $directory = "$OUTBREAK_BASE/" . File::Basename::dirname($reference);
  if ($makedir && defined $URL && length($URL)) {
    if (-d $directory) {
      File::Path::make_path($directory) || die "Cannot make $directory and -makedir is true\n";
    }
    $ff = File::Fetch->new( uri => $URL );
    $ff->fetch( to => $directory );
    if (!defined $ff->output_file || !-f "$directory/" . $ff->output_file) {
      die "Cannot download $URL (to $directory/" . $ff->output_file . ") and -makedir is true\n";
    }
    if ($ff->output_file =~ /\.gz$/) {
      system("gunzip $directory/" . $ff->output_file);
    }
    if (-f "$OUTBREAK_BASE/$reference") {
      $reference = File::Spec->rel2abs("$OUTBREAK_BASE/$reference");
      system("samtools faidx $reference");
      system("bwa index $reference");
    } else {
      die "Error, seemed to download $URL ok to $directory but can't find $OUTBREAK_BASE/$reference\n";
    }
  }
}
if (!-f $reference) {
  $reference = "NONE";
  $directory = "$OUTBREAK_BASE/Others/Noref";
  if (!-d $directory && $makedir) {
    File::Path::make_path($directory);
    if (!-d $directory) {
      die "Cannot make $directory and -makedir is true\n";
    }
  }
}

if (!defined $species || $species eq "") {
  $species = "NONE";
}
if (!defined $fullspecies || $fullspecies eq "") {
  $fullspecies = "NONE";
}

print "WORK=$directory\n";
print "REF=$reference\n";
print "SPECIES=$species\n";
print "FULLSPECIES=\"$fullspecies\"\n";

$sth->finish();
$DBH->disconnect;
