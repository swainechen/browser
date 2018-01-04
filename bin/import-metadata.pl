#!/usr/bin/perl -w
#
# take tab delimited text file
# expected "<species>metadata.txt"
# check for changes and do updates
# database should be germs_browser
# table should be based on species name
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
use Getopt::Long;
use GERMS;

my $db_info = {
  database => "",
  host => "",
  username => "",
  password => ""
};
my $show_help = 0;
my $USE_DB = 0;
my $verbose = 0;
my $force = 0;	# update no matter what - ignore date stamps

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
my %ATTR;	# needed for SQLite connection - ?
my $species_data;
my $species_table;
my $status;

sub print_usage {
  print "Usage: $0 <Speciesmetadata.txt> [ -db ] [ -verbose ] [ -force ]\n";
  print "  Will import/update the data to the Species table based on input filename.\n";
  print "  With -db do database import / update (default just show what would be done)\n";
  print "  With -force then always update\n";
  print "  With -verbose print out all the information collected\n";
}

if ($show_help || !defined $ARGV[0] || !-f $ARGV[0]) {
  &print_usage;
  exit;
}
if (File::Basename::basename($ARGV[0]) =~ /^(\w+)[.-_]?metadata.txt/) {
  $species_table = $1;
} else {
  print STDERR "Can't parse filename ", File::Basename::basename($ARGV[0]), " for species name\n";
  print STDERR "Expect \"speciesmetadata.txt\" as filename\n";
  &print_usage;
  exit;
}

my $DBH = dbconnect($db_info->{database}, $db_info->{host}, $db_info->{username}, $db_info->{password});

# pick off the columns from the first line
$status = ();
my @header = ();
my $datestamp = DateTime->from_epoch(epoch => (stat($ARGV[0]))[9]);
open F, $ARGV[0];
while (<F>) {
  chomp;
  s/\r$//;
  next if /^#/;
  next if /^$/;
  if (scalar @header == 0) {
    @header = split /\t/, $_;
    next;
  }
  $species_data = ();
  $species_data->{DateStamp} = $datestamp;
  @f = split /\t/, $_;
  foreach $i (0..$#header) {
    if (defined $f[$i]) {
      $species_data->{$header[$i]} = $f[$i];
    } else {
      undef $species_data->{$header[$i]};
    }
  }
  $status->{GERMS::do_db($species_data, $species_table, $force, $DBH, $USE_DB)}++;
}
close F;
foreach $i ( qw(INSERT UPDATE NOTHING) ) {
  $status->{$i} = 0 if !defined $status->{$i};
}
if ($verbose) {
  print "For file ", File::Basename::basename($ARGV[0]), "\n";
  print "  DateStamp: $datestamp\n";
  print "  INSERT for $status->{INSERT} lines\n";
  print "  UPDATE for $status->{UPDATE} lines\n";
  print "  NOTHING for $status->{NOTHING} lines\n";
}
