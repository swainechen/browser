#!/usr/bin/perl -w
#
# take tab delimited text file
# expected "/home/slchen/rdrive/ID/ID3/Swaine/Chen lab sequencing libraries.txt"
# check for changes and do updates
# database should be germs_browser
# table should be Studies
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
use Getopt::Long;
use GERMS;

my $show_help = 0;
my $USE_DB = 0;
my $verbose = 0;
my $force = 0;	# update no matter what - ignore date stamps
my $libraries_file = "/home/slchen/rdrive/ID/ID3/Swaine/Chen lab sequencing libraries.txt";

&Getopt::Long::Configure("pass_through");
GetOptions (
  'db!' => \$USE_DB,
  'verbose!' => \$verbose,
  'force!' => \$force,
  'help' => \$show_help
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
my $studies_data;
my $status;

sub print_usage {
  print "Usage: $0 [ -db ] [ -verbose ] [ -force ]\n";
  print "  Will look for the following file:\n";
  print "  $libraries_file\n";
  print "  and import/update the data to the Studies table.\n";
  print "  With -db do database import / update (default just show what would be done)\n";
  print "  With -force then always update\n";
  print "  With -verbose print out all the information collected\n";
  if (!-f $libraries_file) {
    print STDERR "Can't find $libraries_file\n";
  }
}

if ($show_help || !-f $libraries_file) {
  &print_usage;
  exit;
}

my $DBH = GERMS::dbconnect("germs_browser");

# we are expecting these columns:
#   Run VARCHAR(255) NOT NULL,
#   MUX VARCHAR(63),
#   SampleName VARCHAR(255),
#   Machine VARCHAR(63),
#   Paired VARCHAR(4) NOT NULL,
#   ReadLength INT,
#   LibraryPrep VARCHAR(255),
#   SequencingType VARCHAR(63) NOT NULL,
#   ChenLabContact VARCHAR(63) NOT NULL,
#   OtherContact VARCHAR(255),
#   Study VARCHAR(255),
#   PublicAccession VARCHAR(255),
#   PrivacyCode VARCHAR(63),
#   Notes MEDIUMTEXT,
#   DateStamp VARCHAR(63) NOT NULL
#
$status = ();
my @header = ();
my $datestamp = DateTime->from_epoch(epoch => (stat($libraries_file))[9]);
open F, $libraries_file;
while (<F>) {
  chomp;
  s/\r$//;
  next if /^#/;
  next if /^$/;
  if (scalar @header == 0) {
    @header = split /\t/, $_;
    next;
  }
  $studies_data = ();
  $studies_data->{DateStamp} = $datestamp;
  @f = split /\t/, $_;
  foreach $i (0..$#header) {
    if (defined $f[$i]) {
      $studies_data->{$header[$i]} = $f[$i];
    } else {
      undef $studies_data->{$header[$i]};
    }
  }
  $status->{GERMS::do_db($studies_data, "Studies", $force, $DBH, $USE_DB)}++;
}
close F;
foreach $i ( qw(INSERT UPDATE NOTHING) ) {
  $status->{$i} = 0 if !defined $status->{$i};
}
if ($verbose) {
  print "For file $libraries_file\n";
  print "  DateStamp: $datestamp\n";
  print "  INSERT for $status->{INSERT} lines\n";
  print "  UPDATE for $status->{UPDATE} lines\n";
  print "  NOTHING for $status->{NOTHING} lines\n";
}
