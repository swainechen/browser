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
my $operation = $ARGV[1];
if (defined $operation) { $operation = lc($operation); }
my $TIP;
my $date_now = DateTime->now;
my $INSTANCE_ID = `GET http://169.254.169.254/latest/meta-data/instance-id`;
chomp $INSTANCE_ID;
my $INSTANCE_TYPE = `GET http://169.254.169.254/latest/meta-data/instance-type`;
chomp $INSTANCE_TYPE;
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
my $output_log = "";

GetOptions (
  "output=s" => \$output_log,
  "verbose!" => \$verbose
);

if (!defined $runID || !length($runID) ||
    !defined $operation || ($operation ne "start" && $operation ne "finish" && $operation ne "clear" && $operation ne "getlog")) {
  print "Usage: $0 <runID> start [ -verbose ]\n";
  print "       $0 <runID> finish [ -verbose ]\n";
  print "       $0 <runID> clear [ -verbose ]\n";
  print "       $0 <runID> getlog [ -output <log output filename> ] [ -verbose ]\n";
  print "start will insert a row into the Tracking table\n";
  print "finish will check for files and update the Tracking table\n";
  print "clear will remove a row from the Tracking table - this will always ask for confirmation\n";
  print "getlog will get the gzipped log file if it exists and output to specified file (default is runID.log.gz in the current directory)\n";
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
    print "Inserting anyway...\n";
#    exit(1);
  }
  # should be safe to insert now
  $sql = "INSERT INTO $TABLE (Run, InstanceID, InstanceType, Started) VALUES (?, ?, ?, ?)";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID, $INSTANCE_ID, $INSTANCE_TYPE, $date_now) || die "Error: " . $sth->errstr . "\n";
  if ($verbose) {
    print "Inserted row for $runID, started at $date_now on $INSTANCE_ID\n";
  }
}

if ($operation eq "finish") {
  # some sanity - there should be a row there already
  $sql = "SELECT Run, Finished, Species, Success FROM $TABLE WHERE Run = ? AND InstanceID = ? AND InstanceType = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID, $INSTANCE_ID, $INSTANCE_TYPE);
  undef $tempdate;
  $species = "";
  while (@data = $sth->fetchrow_array()) {
    $tempdate = $data[1];
    $species = $data[2] if defined $data[2];
    $tempsuccess = $data[3];
  }
  if ($sth->rows() == 0) {
    die "Operation finish but no existing row for $runID on $INSTANCE_ID ($INSTANCE_TYPE), exiting...\n";
  }
  if (defined $tempdate && length($tempdate)) {
    die "Already have a finished row for $runID on $INSTANCE_ID ($INSTANCE_TYPE), species $species, success $tempsuccess\n";
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
    $sql = "UPDATE $TABLE SET Finished = ?, TIP = ?, Species = ?, Success = ? WHERE Run = ? AND InstanceID = ? AND InstanceType = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($date_now, $TIP, $species, "Y", $runID, $INSTANCE_ID, $INSTANCE_TYPE) || die "Error: " . $sth->errstr . "\n";
    if ($verbose) {
      print "Updated row for $runID on $INSTANCE_ID ($INSTANCE_TYPE), finished successfully at $date_now\n";
    }
  } else {
    # had some error, save the log
    undef $log;
    if (-f "$STAGING/$runID.log") {
      system "gzip $STAGING/$runID.log";
      open F, "$STAGING/$runID.log.gz";
      local $/ = undef;
      $log = <F>;
      close F;
    }
    $sql = "UPDATE $TABLE SET Finished = ?, TIP = ?, Species = ?, Success = ?, LogBlob = ? WHERE Run = ? AND InstanceID = ? AND InstanceType = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($date_now, $TIP, $species, "N", $log, $runID, $INSTANCE_ID, $INSTANCE_TYPE) || die "Error: " . $sth->errstr . "\n";
    if ($verbose) {
      print "Updated row for $runID on $INSTANCE_ID ($INSTANCE_TYPE), finished unsuccessfully at $date_now\n";
      if (defined $log) {
        print "Log at $STAGING/$runID.log saved\n";
      } else {
        print "Could not find $STAGING/$runID.log\n";
      }
    }
  }
}

if ($operation eq "clear") {
  $sql = "SELECT Run, InstanceID, InstanceType, Started, Success FROM $TABLE WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID);
  print "Found the following for $runID in $TABLE Table:\n";
  while (@data = $sth->fetchrow_array()) {
    $data[4] = "" if !defined $data[4];	# all else should be not null
    print join("\t", @data), "\n";
  }
  if ($sth->rows() == 0) {
    print "No rows found for $runID in $TABLE Table. Exiting...\n";
    exit;
  }
  print "Total: ", $sth->rows(), " rows will be deleted\n";
  print "Continue? (y/N): ";
  $command_output = <STDIN>;
  if ($command_output =~ /^y/i) {
    $sql = "DELETE FROM $TABLE WHERE Run = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($runID) || die "Error: " . $sth->errstr . "\n";
    $sql = "SELECT Run, InstanceID FROM $TABLE WHERE Run = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($runID);
    while (@data = $sth->fetchrow_array()) {
      next;
    }
    print $sth->rows(), " rows left for $runID in $TABLE Table\n";
  } else {
    print "Aborting, no changes to $TABLE Table...\n";
  }
}

if ($operation eq "getlog") {
  if ($output_log eq "") {
    $output_log = "$runID.log.gz";
  }
  $sql = "SELECT LogBlob FROM $TABLE WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($runID);
  @data = $sth->fetchrow_array();
  if (defined $data[0]) {
    if (-f $output_log) {
      print "$output_log already exists, refusing to overwrite...\n";
    } else {
      open F, ">$output_log";
      print F $data[0];
      close F;
    }
  } else {
    print "No Log stored in $TABLE Table for $runID\n";
  }
}
