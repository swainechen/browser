#!/usr/bin/perl
#
# Insert tracking information
#
use warnings;
use strict;
use DBI;
use DateTime::Format::DBI;
use DateTime::Format::DateParse;
use GERMS;
use Getopt::Long;
&Getopt::Long::Configure ("pass_through");

my $verbose = 0;
my $runID = $ARGV[0];
my $operation = lc($ARGV[1]);
my $TIP;
my $date_now = DateTime->now;
my $STAGING = $ENV{"STAGING"};
my $DBH;
my $DB_PARSER;
my $TABLE = "Tracker";
my $sql;
my $sth;
my $log;
my @data;
my $count;
my $species;
my $tempdate;
my $command_output;
my $tempsuccess;

GetOptions (
  "verbose!" => \$verbose
);

if (!defined $runID || !length($runID) ||
    !defined $operation || ($operation ne "start" && $operation ne "finish")) {
  print "Usage: $0 <runID> start [ -verbose ]\n";
  print "       $0 <runID> finish [ -verbose ]\n";
  print "start will insert a row into the Tracking table\n";
  print "finish will check for files and update the Tracking table\n";
  exit;
}

$DBH = GERMS::dbconnect();
$DB_PARSER = DateTime::Format::DBI->new($DBH);

if ($operation eq "start") {
  # first check if anything is already there
  $sql = "SELECT Run, Started, Finished, Success FROM $TABLE WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID);
  while (@data = $sth->fetchrow_array()) {
    if (defined $data[3] && $data[3] eq "Y") {
      $tempdate = $data[1] if defined $data[1];
    }
  }
  $count = $sth->rows;
  if ($count) {
    print "$runID has $count entries already\n";
    if (defined $tempdate) {
      print "The run was successful when started on $tempdate\n";
    } else {
      print "But there is no finished date for any successful run\n";
    }
    print "Not inserting\n";
    exit(1);
  }
  # check for completion already
  $command_output = `browser-library-status.pl $runID > /dev/null 2>&1`;
  if ($? == 0) {
    $TIP = GERMS::get_browser_tip($runID, $DBH, 0);
    print "Appears that Run is already done, Tip $TIP\n";
    print "Not inserting\n";
    exit(1);
  }
  # should be safe to insert now
  $sql = "INSERT INTO $TABLE (Run, Started) VALUES (?, ?)";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID, $date_now) || die "Error: " . $sth->errstr . "\n";
  if ($verbose) {
    print "Inserted row for $runID, started at $date_now\n";
  }
}

if ($operation eq "finish") {
  # some sanity - there should be a row there already
  $sql = "SELECT Run, Finished, Species, Success FROM $TABLE WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID);
  undef $tempdate;
  $species = "";
  while (@data = $sth->fetchrow_array()) {
    $tempdate = $data[1];
    $species = $data[2] if defined $data[2];
    $tempsuccess = $data[3];
  }
  if ($sth->rows() == 0) {
    die "Operation finish but no existing row for $runID, exiting...\n";
  }
  if (defined $tempdate && length($tempdate)) {
    die "Already have a finished row for $runID, species $species, success $tempsuccess\n";
  }
  # last parameter is typically $USE_DB but we're not changing anything in Tips
  $TIP = GERMS::get_browser_tip($runID, $DBH, 0);
  $command_output = `browser-library-status.pl $runID > /dev/null 2>&1`;
  undef $species;
  if (defined $TIP) {
    # get species if available
    $sql = "SELECT Kraken FROM Fastq WHERE TIP = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($TIP);
    while (@data = $sth->fetchrow_array()) {
      $species = $data[0];
    }
  }
  if (defined $TIP && $? == 0) {
    # this looks successful - last get species data
    $sql = "UPDATE $TABLE SET Finished = ?, TIP = ?, Species = ?, Success = ? WHERE Run = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($date_now, $TIP, $species, "Y", $runID) || die "Error: " . $sth->errstr . "\n";
    if ($verbose) {
      print "Updated row for $runID, finished successfully at $date_now\n";
    }
  } else {
    # had some error, save the log
    undef $log;
    if (-f "$STAGING/$runID.log") {
      system "gzip $STAGING/$runID.log";
      open F, "$STAGING/$runID.log.gz";
      $log = <F>;
      close F;
    }
    $sql = "UPDATE $TABLE SET Finished = ?, TIP = ?, Species = ?, Success = ?, LogBlob = ? WHERE Run = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($date_now, $TIP, $species, "N", $log, $runID) || die "Error: " . $sth->errstr . "\n";
    if ($verbose) {
      print "Updated row for $runID, finished unsuccessfully at $date_now\n";
      if (defined $log) {
        print "Log at $STAGING/$runID.log saved\n";
      } else {
        print "Could not find $STAGING/$runID.log\n";
      }
    }
  }
}
