#!/usr/bin/perl -w
#
# Take in a run number
# Print out species, reference sequence
#
use warnings;
use strict;
use File::Spec;
use File::Basename;
use GERMS;
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

# this needs to get set
my $NOREF_DIR = ""

my $sid = 'germs_browser';
my $user;
my $pass;
my $host;
my $quiet = 1;
my $force_species = "NONE";
my $force_reference = "NONE";
GetOptions (
  'sid|db=s' => \$sid,
  'quiet!' => \$quiet,
  'user=s' => \$user,
  'pass=s' => \$pass,
  'host=s' => \$host,
  'species=s' => \$force_species,
  'reference=s' => \$force_reference
);

if (!defined $ARGV[0] || !length($ARGV[0])) {
  print "Usage: $0 <Run ID> [ -species <species> ] [ -reference <reference> ]\n";
  print "Will print out:\n";
  print "  Directory for rest of files\n";
  print "  Species (abbreviation) for SRST2\n";
  print "  Reference sequence for mapping\n";
  print "Can force Species, use \"Ecoli\" or \"Escherichia coli\" (quoted as needed)\n";
  print "Can force Reference, use Reference filename, Reference md5sum, Reference genome name (like \"EC958\" or \"LT2\")\n";
  exit;
}

my @f;
my @data;
my $directory = "";
my $fullspecies = "";	# this is like "Escherichia coli"
my $species = "";	# this is like "Ecoli"
my $reference = "";
my $dbh = GERMS::dbconnect($sid, $host, $user, $pass);
my $sql;
my $sth;
my $temp;

if (uc($force_species) eq "NONE") {
  $sql = "SELECT Kraken FROM Fastq JOIN Tips ON Fastq.TIP = Tips.TIP WHERE Run = ?";
  $sth = $dbh->prepare($sql);
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
    $sth = $dbh->prepare($sql);
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
  $sql = "SELECT ReferenceFile FROM ReferenceGenomes WHERE DefaultReference = 1 AND Species = ?";
  $sth = $dbh->prepare($sql);
  $sth->execute($species);
  while (@data = $sth->fetchrow_array()) {
  if (length $data[0]) {
      $reference = File::Spec->rel2abs($data[0]);
      last;
    }
  }
} else {
  if (-f $force_reference) {
    $reference = File::Spec->rel2abs($force_reference);
  } else {
    # this should enable using just the base file name, md5sum,
    # or reference name as $force_reference
    $sql = "SELECT Species, ReferenceFile from ReferenceGenomes WHERE ReferenceFile LIKE ? OR ReferenceMD5 = ? OR ReferenceName = ?";
    $sth = $dbh->prepare($sql);
    $sth->execute("%" . $force_reference, $force_reference, $force_reference);
    while (@data = $sth->fetchrow_array()) {
      if (length $data[0] && length $data[1]) {
        $species = $data[0];
        $reference = File::Spec->rel2abs($data[1]);
        last;
      }
    }
  }
}
if (length $reference && -f $reference) {
  $directory = File::Basename::dirname($reference);
} else {
  $reference = "NONE";
  $directory = $NOREF_DIR;
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
$dbh->disconnect;
