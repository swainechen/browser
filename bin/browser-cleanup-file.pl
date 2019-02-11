#!/usr/bin/perl -w
#
# clean up various things
#
use warnings;
use strict;
use File::Spec;
use File::Basename;
use DateTime;
use DateTime::Format::DateParse;
use GERMS;
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $quiet = 0;
my $show_help = 0;
my $USE_DB = 0;
my $verbose = 1;
my $filetype = "";
my $all = 0;
my $force = 0;
GetOptions (
  'help' => \$show_help,
  'quiet' => \$quiet,
  'verbose!' => \$verbose,
  'type=s' => \$filetype,
  'all!' => \$all,
  'db!' => \$USE_DB,
  'force!' => \$force
);

if (!defined $ARGV[0] || !length($ARGV[0])) {
  $show_help = 1;
}
if ($show_help) {
  print "Usage: $0 <RunID> -type <filetype> [ -all ] [ -verbose ]\n";
  print "Usage: $0 <md5sum> [ -verbose ]\n";
  print "  Will clean up old versions of files for the given file type or file specified by the MD5sum\n";
  print "  For srst2 will clear out any information associated with that file in Genes, Resistance, MLST\n";
  exit;
}

my $DBH = GERMS::dbconnect("germs_browser");
my $table_list;
@{$table_list->{"srst2.gz"}} = qw(Genes Resistance MLST);
@{$table_list->{"tgz"}} = qw(Assembly);
@{$table_list->{"gcov.gz"}} = ();
@{$table_list->{"lofreq.gz"}} = ();
my $runID;
my $MD5;
my $TIP;
my $type;
my $sql;
my $sth;
my @data;
my $table;
my $DateStamp;
my $latest = "";
my $FileTimes;
my $allfiles;
my $fileinfo;
my $meta;
my @other;

# if we get an MD5 just do one file
if ($ARGV[0] =~ /[0-9a-fA-F]{32}/) {
  $MD5 = $ARGV[0];
  ($TIP, $type, $DateStamp) = tiptype_by_md5($MD5);
  if ($TIP == -1) {
    die "Can't find any TIP or Type in Files for $MD5\n";
  }
  if (!defined $runID || $runID eq "") {
    $sql = "SELECT Run FROM Tips where TIP = ?";
    $sth = $DBH->prepare($sql);
    $sth->execute($TIP);
    while (@data = $sth->fetchrow_array()) {
      $runID = $data[0];
      last;
    }
  }
  $fileinfo->{TYPE} = $type;
  $fileinfo->{MD5} = $MD5;
  $fileinfo->{DATESTAMP} = $DateStamp;
  $FileTimes = get_times($TIP, $type);
  if (defined $FileTimes->{LATEST} && defined $FileTimes->{$MD5} &&
      $FileTimes->{LATEST} eq $FileTimes->{$MD5}) {
    $latest = " LATEST";
  }
  $meta = get_metadata($TIP, $runID);
  if ($verbose) {
    print "For $runID (Tip $TIP; ", meta_string($meta), "), Filetype $type\n";
    print "  MD5 $MD5 ($DateStamp)$latest:\n";
    foreach $table (@{$table_list->{$type}}) {
      print "    Table $table: Would delete ", rows_by_md5($TIP, $table, $type, $MD5), " rows\n";
    }
  }
  &do_one_file($TIP, $fileinfo, $force);
  if ($verbose) {
    if (scalar(keys %$FileTimes) > 2) {
      # should have LATEST and the file we're using
      my $other = scalar(keys %$FileTimes) - 2;
      print "Other files of this type ($other):\n";
      foreach $other (sort { DateTime->compare(dt($FileTimes->{$a}), dt($FileTimes->{$b})) } keys %$FileTimes) {
        next if $other eq $MD5;
        next if $other eq "LATEST";
        if ($FileTimes->{$other} eq $FileTimes->{LATEST}) {
          print "  MD5 $other ($FileTimes->{$other}) LATEST\n";
        } else {
          print "  MD5 $other ($FileTimes->{$other})\n";
        }
      }
    }
  }
  exit;
}

# we should need a run at this point
$runID = $ARGV[0];
$TIP = GERMS::get_browser_tip($runID, $DBH, $USE_DB);
$allfiles = files_by_tip($TIP);

# check if we only wanted one file type
if (uc($filetype) =~ /^SRST/) {
  foreach $type (keys %$allfiles) {
    if ($type ne "srst2.gz") {
      delete $allfiles->{$type};
    }
  }
} elsif (uc($filetype) =~ /^TGZ/) {
  foreach $type (keys %$allfiles) {
    if ($type ne "tgz") {
      delete $allfiles->{$type};
    }
  }
}

$meta = get_metadata($TIP, $runID);
print "For $runID (Tip $TIP; ", meta_string($meta), "):\n" if $verbose;
foreach $type (sort keys %$allfiles) {
  my $num = scalar(keys %{$allfiles->{$type}});
  my $count = 0;
  if ($verbose) {
    print "Type: $type ($num total)\n";
    foreach $MD5 (sort { DateTime->compare(dt($allfiles->{$type}->{$a}->{DateStamp}), dt($allfiles->{$type}->{$b}->{DateStamp})) } keys %{$allfiles->{$type}}) {
      $count++;
      if ($count == $num) {
        print "  MD5 $MD5 ($allfiles->{$type}->{$MD5}->{DateStamp}) LATEST\n";
        print "    ", short_filename($allfiles->{$type}->{$MD5}->{Filename}), "\n";
      } else {
        print "  MD5 $MD5 ($allfiles->{$type}->{$MD5}->{DateStamp})\n";
        print "    ", short_filename($allfiles->{$type}->{$MD5}->{Filename}), "\n";
      }
    }
    print "\n";
  }
  $count = 0;
  foreach $MD5 (sort { DateTime->compare(dt($allfiles->{$type}->{$a}->{DateStamp}), dt($allfiles->{$type}->{$b}->{DateStamp})) } keys %{$allfiles->{$type}}) {
    $count++;
    if ($count != $num) {
      print "  MD5 $MD5 ($allfiles->{$type}->{$MD5}->{DateStamp})\n" if $verbose;
      $fileinfo->{TYPE} = $type;
      $fileinfo->{MD5} = $MD5;
      $fileinfo->{DATESTAMP} = $allfiles->{$type}->{$MD5};
      &do_one_file($TIP, $fileinfo, $force);
    }
  }
}
exit;

if (defined $TIP && defined $type && defined $MD5) {
  $fileinfo->{TYPE} = $type;
  $fileinfo->{MD5} = $MD5;
  $fileinfo->{DATESTAMP} = $DateStamp;
  &do_one_file($TIP, $fileinfo, $force);
}

sub do_one_file {
  my ($TIP, $fileinfo, $force) = @_;
  my $type = $fileinfo->{TYPE};
  my $MD5 = $fileinfo->{MD5};
  my $DateStamp = $fileinfo->{DATESTAMP};
  my $FileTimes = get_times($TIP, $type);
  if (defined $FileTimes->{LATEST} && defined $FileTimes->{$MD5} &&
      $FileTimes->{LATEST} eq $FileTimes->{$MD5}) {
    $latest = " LATEST";
  }
  if ($USE_DB) {
    foreach $table (@{$table_list->{$type}}) {
      &check_delete_other($table, $runID, $TIP, $type, $MD5, $force);
    }
    &check_delete_file($runID, $TIP, $type, $MD5, $force);
  } elsif ($verbose) {
    foreach $table (@{$table_list->{$type}}) {
      print "    Would delete ", rows_by_md5($TIP, $table, $type, $MD5), " rows from $table\n";
    }
    print "    Would delete 1 row from Files\n";
  }
}

sub tiptype_by_md5 {
  my ($MD5) = @_;
  my @data;
  my $sql = "SELECT TIP, Type, DateStamp FROM Files WHERE MD5 = ?";
  my $sth = $DBH->prepare($sql);
  $sth->execute($MD5);
  while (@data = $sth->fetchrow_array) {
    return ($data[0], $data[1], $data[2]);
  }
  return (-1, "", "");
}

sub md5_by_type {
  my ($TIP, $type) = @_;
  my @data;
  my $sql = "SELECT MD5, DateStamp FROM Files WHERE TIP = ? AND Type = ? ORDER by DateStamp";
  my $sth = $DBH->prepare($sql);
  $sth->execute($TIP, $type);
  while (@data = $sth->fetchrow_array) {
    return ($data[0], $data[1]);
  }
  return (-1, "");
}

sub rows_by_md5 {
  my ($TIP, $table, $type, $MD5) = @_;
  my @data;
  my $sql;
  my $sth;
  $sql = "SELECT COUNT(*) FROM $table WHERE TIP = ? AND SourceFileType = ? AND SourceFileMD5 = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($TIP, $type, $MD5);
  while (@data = $sth->fetchrow_array) {
    return $data[0];
  }
  return -1;
}

sub get_times {
  my ($TIP, $type) = @_;
  my $times;
  my @data;
  my $latest = "";
  my $latest_dt;
  my $sql = "SELECT MD5, DateStamp FROM Files WHERE TIP = ? AND Type = ?";
  my $sth = $DBH->prepare($sql);
  $sth->execute($TIP, $type);
  while (@data = $sth->fetchrow_array) {
    $times->{$data[0]} = $data[1];
    if ($latest eq "") {
      $latest = $data[1];
      $latest_dt = dt($data[1]);
    } elsif (DateTime->compare(dt($data[1]), $latest_dt) == 1) {
      $latest = $data[1];
      $latest_dt = dt($data[1]);
    }
  }
  if ($latest ne "") {
    $times->{LATEST} = $latest;
  }
  return $times;
}

sub files_by_tip {
  my ($TIP) = @_;
  my $sql = "SELECT Type, MD5, Filename, DateStamp FROM Files WHERE TIP = ?";
  my $files;
  my @data;
  my $sth = $DBH->prepare($sql);
  $sth->execute($TIP);
  while (@data = $sth->fetchrow_array) {
    $files->{$data[0]}->{$data[1]}->{DateStamp} = $data[3];
    $files->{$data[0]}->{$data[1]}->{Filename} = $data[2];
  }
  return $files;
}

sub dt {
  # take a string / DateStamp, return DateTime object
  my ($s) = @_;
  return(DateTime::Format::DateParse->parse_datetime($s));
}

sub check_delete_other {
  my ($table, $runID, $TIP, $type, $MD5, $force) = @_;
  my $base = "FROM $table WHERE TIP = $TIP AND SourceFileType = '$type' AND SourceFileMD5 = '$MD5'";
  my $check = "SELECT COUNT(*) $base";
  my $delete = "DELETE $base";
  my $response = "";
  my $num = 0;
  my @data;
  my $sth = $DBH->prepare($check);
  $sth->execute;
  while (@data = $sth->fetchrow_array) {
    $num = $data[0];
    last;
  }
  return(0) if !$num;
  if (!$force) {
    print STDERR "About to delete $num rows from $table for $runID ($TIP) $type ($MD5)\n";
    print STDERR "Proceed (y/N)? ";
    $response = <STDIN>;
    return(0) if $response !~ /^y/i;
  }
  $sth = $DBH->prepare($delete);
  $sth->execute();
  $num = $sth->rows;
  print STDERR "Deleted $num rows from $table\n";
  return($num);
}

sub check_delete_file {
  my ($runID, $TIP, $type, $MD5, $force) = @_;
  my $base = "FROM Files WHERE TIP = $TIP AND Type = '$type' AND MD5 = '$MD5'";
  my $check = "SELECT COUNT(*) $base";
  my $delete = "DELETE $base";
  my $response = "";
  my $num = 0;
  my @data;
  my $sth = $DBH->prepare($check);
  $sth->execute;
  while (@data = $sth->fetchrow_array) {
    $num = $data[0];
    last;
  }
  return(0) if !$num;
  if (!$force) {
    print STDERR "About to delete $num rows from Files for $runID ($TIP) $type ($MD5)\n";
    print STDERR "Proceed (y/N)? ";
    $response = <STDIN>;
    return(0) if $response !~ /^y/i;
  }
  $sth = $DBH->prepare($delete);
  $sth->execute();
  $num = $sth->rows;
  print STDERR "Deleted $num rows from Files\n";
  return($num);
}

sub get_metadata {
  my ($TIP, $runID) = @_;
  # get some metadata to give some context
  my $meta;
  my @data;
  my $sql;
  my $sth;

  $sql = "SELECT Kraken, KrakenPercent FROM Fastq WHERE TIP = $TIP";
  $sth = $DBH->prepare($sql);
  $sth->execute;
  while (@data = $sth->fetchrow_array) {
    $meta->{KRAKEN} = $data[0];
    $meta->{KRAKENPERCENT} = $data[1];
    last;
  }

  $sql = "SELECT SampleName, Study, PrivacyCode FROM Studies WHERE Run = '$runID'";
  $sth = $DBH->prepare($sql);
  $sth->execute;
  while (@data = $sth->fetchrow_array) {
    $meta->{SAMPLENAME} = $data[0];
    $meta->{STUDY} = $data[1];
    $meta->{PRIVACYCODE} = $data[2];
    last;
  }

  return($meta);
}

sub meta_string {
  my ($meta) = @_;
  # build a text string to show some metadata
  my @i = ();
  if (defined $meta->{SAMPLENAME}) {
    push @i, $meta->{SAMPLENAME};
  }
  if (defined $meta->{STUDY}) {
    push @i, $meta->{STUDY};
  }
  if (defined $meta->{PRIVACYCODE}) {
    push @i, $meta->{PRIVACYCODE};
  }
  if (defined $meta->{KRAKEN}) {
    push @i, "$meta->{KRAKEN} ($meta->{KRAKENPERCENT}%)";
  }
  if (scalar @i) {
    return (join ("; ", @i));
  } else {
    return ("");
  }
}

sub short_filename {
  my ($filename) = @_;
  if ($filename =~ /^s3/) {
    $filename =~ s/^.*Outbreaks\///;
  } else {
    $filename =~ s/^.*Outbreaks\///;
  }
  return($filename);
}
