#!/usr/bin/perl -w
#
# source file is a Seqsero_result.txt file
# check for changes and do updates
# database should be germs_browser
# table should be SENTERICA_Serotype
#
use lib '/home/slchen/bin';
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
my $runID = "";	# we should try to read this from the Seqsero_result.txt file

&Getopt::Long::Configure("pass_through");
GetOptions (
  'run=s' => \$runID,
  'db!' => \$USE_DB,
  'verbose!' => \$verbose,
  'force!' => \$force,
  'help' => \$show_help
);

my $data;
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
my $input;
my $Oag;
my $H1;
my $H2;
my $AgProfile;
my $Sdf;
my $Serotype;

sub print_usage {
  print "Usage: $0 <Seqsero_result.txt file> [ -run <RunID> ] [ -db ] [ -verbose ] [ -force ]\n";
  print "  Will import/update the data to the SENTERICA_Serotype table.\n";
  print "  With -db do database import / update (default just show what would be done)\n";
  print "  With -force then always update\n";
  print "  With -verbose print out all the information collected\n";
}

if ($show_help) {
  &print_usage;
  exit;
}
if (!defined $ARGV[0] || !-f $ARGV[0]) {
  &print_usage;
  exit;
}

open S, $ARGV[0];
while (<S>) {
  next if /^$/;
  next if /^#/;
  if (/^Input files:(.*)/) {
    $input = $1;
    $input =~ s/^\s+//;
    $input =~ s/\s+$//;
    if ($runID eq "") {
      $runID = $input;
      $runID =~ s/\.fna$//;
      $runID =~ s/\.assembly$//;
      if ($runID =~ /(^GC.*?\.\d+)/) {
        $runID = $1;
      }
    }
  }
  if (/^O antigen prediction:(.*)/) {
    $Oag = $1;
    $Oag =~ s/^\s+//;
    $Oag =~ s/\s+$//;
  }
  if (/^H1 antigen prediction\(fliC\):(.*)/) {
    $H1 = $1;
    $H1 =~ s/^\s+//;
    $H1 =~ s/\s+$//;
  }
  if (/^H2 antigen prediction\(fljB\):(.*)/) {
    $H2 = $1;
    $H2 =~ s/^\s+//;
    $H2 =~ s/\s+$//;
  }
  if (/^Predicted antigen profile:(.*)/) {
    $AgProfile = $1;
    $AgProfile =~ s/^\s+//;
    $AgProfile =~ s/\s+$//;
  }
  if (/^Sdf prediction:(.*)/) {
    $Sdf = $1;
    $Sdf =~ s/^\s+//;
    $Sdf =~ s/\s+$//;
  }
  if (/^Predicted serotype\(s\):(.*)/) {
    $Serotype = $1;
    $Serotype =~ s/^\s+//;
    $Serotype =~ s/\s+$//;
  }
}
if ($runID eq "") {
  die "No Run ID given and can't figure it out from $input\n";
}
my $DBH = GERMS::dbconnect("germs_browser");
my $DB_PARSER = DateTime::Format::DBI->new($DBH);
my $TIP = GERMS::get_browser_tip($runID, $DBH, $USE_DB);

$data->{TIP} = $TIP;
if (defined $Oag) {
  $data->{OAntigen} = $Oag;
}
if (defined $H1) {
  $data->{H1} = $H1;
}
if (defined $H2) {
  $data->{H2} = $H2;
}
if (defined $AgProfile) {
  $data->{AntigenicProfile} = $AgProfile;
}
if (defined $Sdf) {
  $data->{Sdf} = $Sdf;
}
if (defined $Serotype) {
  $data->{Serotype} = $Serotype;
}

if ($verbose) {
  print "RunID: $runID\n";
  print "TIP: $TIP\n";
  if ($USE_DB) {
    print "Doing database operations (-db is on)\n";
  } else {
    print "No database operations (no -db flag)\n";
  }
}

$status = GERMS::do_db($data, "SENTERICA_Serotype", $force, $DBH, $USE_DB);
if ($verbose) {
  foreach $i (sort keys %$data) {
    print "  $i: $data->{$i}\n";
  }
  print "Serotype: $status\n";
}
