#!/usr/bin/perl -w
#
# try to figure things out
# GIS is [DWR][A-Z]B###(#) usually
# URL should go to GenBank SRA
#
# use GERMS_DATA environment variable
#
use File::Basename;
use File::Spec;
use Cwd;
use FindBin;
use Getopt::Long;
use DBI;
&Getopt::Long::Configure("pass_through");

# needs configuration
my $S3_BASE = "";
my $S3_PROFILE = "";

# some globals
my $SRADB = $FindBin::Bin . "/../lib/ncbi/SRAmetadb.sqlite";
my $DBH;
my %ATTR;
my $TEMPDIR = "/tmp";
my $FQD_DOWNLOAD = 0;

my $dir = "";
my $fqd = "fastq-dump";
my $combine = 1;	# whether to combine multiple files
my $delimiter = "\n";	# for output
my $delete = -1;	# whether to delete files after combining
			# leave this at -1 to figure out if user set something
my $wait = 1;		# run job with -sync y so it doesn't return until files
			# are there already
my $debug = 0;
my @out;
my $out;
my $tries = 5;
my $i;		# generic counter
my $study;
my $sample;
my $url;
my $current;
my @f;
my $checkonly = 0;	# just look for files, no download
my $showhelp = 0;
my $internal = 0;	# use this to redirect combining to an ionode job
my $command;		# temp to hold the command to be run
my $base_command = File::Spec->rel2abs($0);
my $ENA_DOWNLOAD = `which ena-fast-download.py`;
chomp $ENA_DOWNLOAD;

GetOptions (
  'debug!' => \$debug,
  'combine!' => \$combine,
  'delete!' => \$delete,
  'delimiter=s' => \$delimiter,
  'tries=i' => \$tries,
  'wait!' => \$wait,
  'checkonly' => \$checkonly,
  'internal!' => \$internal,
  'fqd_download!' => \$FQD_DOWNLOAD,
  'help' => \$showhelp
);

if (defined $ENV{GERMS_DATA}) {
  $dir = $ENV{GERMS_DATA};
} else {
  die "Please set GERMS_DATA environment variable\n";
}
if ($showhelp || !defined $ARGV[0] || !length($ARGV[0])) {
  print "Usage: $0 <ID> [ -delimiter <char> ] [ -debug ] [ -combine|nocombine ] [ -delete|nodelete ] [ -nowait ] [ -hrt <hrt string> ] [ -mem <mem_free string> ] [ -checkonly ] [ -tries <int> ] [ -fqd_download|nofqd_download ]\n";
  print "ID is either some GIS sequencing library ID or a Genbank SRA URL.\n";
  print "This script relies on the GERMS_DATA environment variable being set\n";
  print "Default is to combine - for example HiSeq4K comes in 4 fastq files otherwise\n";
  print "Default is to delete files that are combined, leaving only one set of ID-combined_R[12].fastq.gz files\n";
  print "Default delimiter is newline, this is the separator for multiple filenames in output. This can also be a string, one good one is to use -delimiter \" -q2 \" in a command line.\n";
  print "-checkonly option only looks for files and won't download anything. Will return error code 1 if no files found (i.e. not yet downloaded).\n";
  print "Default tries is 5 (number of retries for downloading from Genbank)\n";
  print "-fqd_download uses fastq-dump to download sra files from Genbank (default is not to use fastq-dump)\n";
  exit;
}

if ($internal) {
  # this should only be for combine
  if ($combine) {
    @out = combine($ARGV[0], $dir, $delete, $internal);
    print join ("\n", @out), "\n";
    exit;
  }
}

$current = Cwd::getcwd();
chdir $dir;
@out = q12($ARGV[0], $dir);
if (scalar @out) {
  if (!$checkonly && $combine && scalar @out > 2) {
    # don't delete unless user said otherwise
    if ($delete == -1) { $delete = 0; }
    # this will be not internal for combine by definition
    @out = combine($ARGV[0], $dir, $delete, 0);
  }
  print join ($delimiter, @out), "\n";
  exit;
}

if ($ARGV[0] =~ /[SED]RR\d+/ && -x $ENA_DOWNLOAD) {
  # we can handle this with ena-fast-download
  $run = $ARGV[0];
  @out = sra_files($run, $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
    exit;
  }
  $command = "mkdir -p $dir/$run && $ENA_DOWNLOAD $run --output-directory $dir/$run --quiet > /dev/null 2>&1";
  if ($debug) {
    print STDERR "Trying to get files from ENA. Running command:\n";
    print STDERR "$command\n";
  }
  if ($checkonly) {
    exit(1);
  } else {
    system ($command);
  }
  if ($? != 0) {
    print STDERR "Error from $ENA_DOWNLOAD. Falling back to Genbank.\n";
  } else {
    @out = sra_files($run, $dir);
    if (scalar @out) {
      print join ($delimiter, @out), "\n";
      exit;
    }
  }
}
if ($ARGV[0] =~ /[SED]R[APSXR]\d+/) {
  $url = "";
  $db_query = sra_query($ARGV[0]);
  if (!defined $db_query || scalar keys %$db_query == 0) {
    die "Can't find data for $ARGV[0].\n";
  }
  if (scalar keys %$db_query > 1) {
    print STDERR "Ambiguous identifier $ARGV[0] - potential matches:\n";
    foreach $i (keys %$db_query) {
      print "$i\t$db_query->{$i}\t", make_url($i, $db_query->{$i});
    }
    die;
  }
  foreach $i (keys %$db_query) {
    $run = $i;
    $url = make_url($i, $db_query->{$i});
    last;
  }
  if (!length $url) {
    die "Couldn't figure out what $ARGV[0] is.\n";
  }
  print STDERR "$url\n" if $debug;
  @out = sra_files($run, $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
    exit;
  }
  # first check if it's on S3
  $command = "aws --profile $S3_PROFILE s3 ls $S3_BASE/$ARGV[0]/ > /dev/null";
  system ($command);
  if ($? != 0) {
    if (-x $ENA_DOWNLOAD) {
      $command = "mkdir -p $dir/$run && $ENA_DOWNLOAD $run --output-directory $dir/$run --quiet > /dev/null 2>&1";
      if ($debug) {
        print STDERR "Trying to get files from ENA. Running command:\n";
        print STDERR "$command\n";
      }
      if ($checkonly) {
        exit(1);
      } else {
        system ($command);
      }
      if ($? != 0) {
        print STDERR "Error from $ENA_DOWNLOAD. Falling back to Genbank.\n";
      } else {
        @out = sra_files($run, $dir);
        if (scalar @out) {
          print join ($delimiter, @out), "\n";
          exit;
        }
      }
    }
    if ($FQD_DOWNLOAD) {
      # let sra utils take care of the urls - can remove make_url procedure and reliance on it now
      $command = "$fqd --split-files --gzip -O $dir/$run $run > /dev/null 2>&1 ; rm -f ~/ncbi/public/sra/$run.sra";
    } else {
      $command = "wget $url --tries=$tries --quiet -O $TEMPDIR/$run.sra && $fqd --split-files --gzip -O $dir/$run $TEMPDIR/$run.sra && rm -f $TEMPDIR/$run.sra > /dev/null 2>&1";
    }
    if ($debug) {
      print STDERR "Trying to get files from GenBank SRA. Running command:\n";
      print STDERR "$command\n";
    }
  } else {
    $command = "mkdir -p $dir/$ARGV[0]; cd $dir/$ARGV[0]; aws --profile $S3_PROFILE s3 sync $S3_BASE/$ARGV[0]/ . > /dev/null 2>&1";
    if ($debug) {
      print STDERR "Trying to get files from AWS S3. Running command:\n";
      print STDERR "$command\n";
    }
  }
  if ($checkonly) {
    exit 1;
  } else {
    system ($command);
  }
  if ($? != 0) {
    die "Some error ($?) with SRA job...\n";
  }
  @out = sra_files($run, $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
  }
} elsif ($ARGV[0] =~ /[DMWR][A-Z][A-Z]\d\d\d+/) {
  $command = "mkdir -p $dir/$ARGV[0]; cd $dir/$ARGV[0]; aws --profile $S3_PROFILE s3 sync $S3_BASE/$ARGV[0]/ . > /dev/null 2>&1";
  if ($debug) {
    print STDERR "Trying to get GIS style files from S3. Running command:\n";
    print STDERR "$command\n";
  }
  if ($checkonly) {
    exit 1;
  } else {
    system($command);
  }
  if ($? != 0) {
    die "Some error ($?) with S3 job...\n";
  }
  @out = q12($ARGV[0], $dir);
  if (scalar(@out) == 0) {
    if ($debug) {
      print STDERR "No error in command, but no files yet, params $ARGV[0] and $dir. Waiting up to 30 seconds...\n";
    }
    $i = 0;
    while (scalar(@out) == 0) {
      sleep 5;
      $i += 5;
      @out = q12($ARGV[0], $dir);
      last if $i >= 30;
    }
  }
  if ($combine && scalar @out > 2) {
    # delete source files unless user said otherwise
    if ($delete == -1) { $delete = 1; }
    @out = combine($ARGV[0], $dir, $delete);
  }
  print join ($delimiter, @out), "\n";
} else {
  @out = sra_files($ARGV[0], $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
    exit;
  }
  # check if it's on S3
  $command = "aws --profile $S3_PROFILE s3 ls $S3_BASE/$ARGV[0]/ > /dev/null";
  system ($command);
  if ($? == 0) {
    $command = "mkdir -p $dir/$ARGV[0]; cd $dir/$ARGV[0]; aws --profile $S3_PROFILE s3 sync $S3_BASE/$ARGV[0]/ . > /dev/null 2>&1";
    if ($checkonly) {
      exit 1;
    } else {
      system ($command);
    }
    if ($? != 0) {
      die "Some error ($?) with AWS S3 copy...\n";
    }
  }
  chdir $current;
  @out = sra_files($ARGV[0], $dir);
  if (scalar @out) {
    print join ($delimiter, @out), "\n";
  } else {
    die "Couldn't figure out what $ARGV[0] is.\n";
  }
}
chdir $current;

sub sra_files {
  my ($run, $dir) = @_;
  my $f;
  my @r = ();
  my @r2;
  my $size;
  my $maxsize = 0;
  my $filename;
  if (!-d "$dir/$run") {
    return ();
  }
  opendir(D, "$dir/$run");
  while ($f = readdir D) {
    if ($f =~ /($run)_\d\.fastq\.gz/) {
      $filename = "$dir/$run/$f";
      push @r, $filename;
      $size->{"$filename"} = (stat("$filename"))[7];
      $maxsize = $size->{"$filename"} if $maxsize < $size->{"$filename"};
    }
  }
  @r = sort { $a cmp $b } (@r);
  if (scalar @r == 2) {
    foreach $f (@r) {
      if ($size->{$f} / $maxsize > 0.5) {
        # for some miseq they seem to be trimmed already
        # filesizes can be quite different
        push @r2, $f;
      }
    }
  } else {
    foreach $f (@r) {
      if ($size->{$f} / $maxsize > 0.7) {
        # for some sra data sets there are indexes or other things
        # and the right files appear to be _1 and _3 or something else
        # try to take the biggest files here and hope that works
        push @r2, $f;
      }
    }
  }
  return @r2;
}

sub sra_query {
  # plan to deprecate in favor of ENA - no extraction needed
  my ($search) = @_;
  my $q;
  my $sql;
  if (!defined $DBH) {
    $DBH = DBI->connect("DBI:SQLite:dbname=$SRADB", "", "", \%ATTR);
  }
  if ($search =~ /^.RR\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE run_accession = ?";
  } elsif ($search =~ /^.RX\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE experiment_accession = ?";
  } elsif ($search =~ /^.RS\d/) {
    $sql = "SELECT run_accession, experiment_accession FROM sra WHERE sample_accession = ?";
  }
  my $sth = $DBH->prepare($sql);
  $sth->execute($search);
  my $r;
  my @r;
  while (@r = $sth->fetchrow_array) {
    $r->{$r[0]} = $r[1];
  }
  return $r;
}

sub q12 {
  my ($lib, $dir) = @_;
  my $i = 0;
  my $j = 0;
  my @return = ();
  my @filelist = ();
  my $file;
  my @combined = ();
  my $out_dir = "";
  return @return if !defined $lib || !length $lib || !-d $dir;
  if (-d "$dir/$lib.h5") {
    opendir D, "$dir/$lib.h5";
    $out_dir = "$lib.h5";
  } elsif (-d "$dir/$lib") {
    opendir D, "$dir/$lib";
    $out_dir = "$lib";
  }
  if ($out_dir ne "") {
    while ($file = readdir(D)) {
      push @filelist, $file;
    }
    @filelist = sort { $a cmp $b } @filelist;
    foreach $file (@filelist) {
      if (
          # standard HiSeq
          $file =~ /_R[12]\.fastq$/     || $file =~ /_R[12]\.fastq.gz$/ ||
          $file =~ /_R[12]_\d+\.fastq$/ || $file =~ /_R[12]_\d+\.fastq.gz$/ ||
          # NextSeq at GIS looks like WEB1880-CGATGT_S2_L001_R1_001_NS001-PE-R00488.fastq.gz (July 2018)
          # HiSeq at GIS looks like WEB3019-AGCAGTAC_S94_L007_R1_001_HS008-PE-R00095.fastq.gz (July 2018)
          # or WEB2675_HS007-PE-R00300-TGCTACAT_S1_L001_R1_001.fastq.gz
          # iSeq is WEB2675-TGCTACAT_S1_L001_R1_001_IS001-PE-R00002.fastq.gz
          $file =~ /_R[12]_\d+_\w+\d+-[PS]E-R\d+.fastq$/ ||
          $file =~ /_R[12]_\d+_\w+\d+-[PS]E-R\d+.fastq.gz$/ ||
          # can also see WEB4055_HS008-PE-R00136_L004_R1_unaligned_001.fastq.gz
          $file =~ /_R[12]_unaligned_\d+.fastq.gz$/ ||

          # SRA from Genbank
          $file =~ /_[12]\.fq$/      || $file =~ /_[12]\.fq.gz$/
         ) {
        push @return, "$dir/$out_dir/$file";
        push @combined, "$dir/$out_dir/$file" if $file =~ /-combined_R[12].fastq.gz$/;
      }
    }
  }
  if (scalar @combined == 2) {
    return @combined;
  } else {
    return (@return);
  }
}

sub combine {
  # need to know whether to delete files
  my ($lib, $dir, $delete, $internal) = @_;
  my @files = q12($lib, $dir);
  my $r1in = "";
  my $r2in = "";
  my $r1out;
  my $r2out;
  my $i;
  my @return;
  my $command;
  my $delete_param = "-nodelete";
  my $cluster_out;
  my $cluster_err;

  if (!$internal) {
    if ($debug) {
      $cluster_out = "$dir/$ARGV[0]-get_files.o";
      $cluster_err = "$dir/$ARGV[0]-get_files.e";
    } else {
      $cluster_out = "/dev/null";
      $cluster_err = "/dev/null";
    }
    if ($delete) { $delete_param = "-delete"; }
    # have to manage environment and absolute path to our initial program
    $command = "$base_command $ARGV[0] -internal -combine $delete_param";
    if ($debug) {
      print STDERR "$command\n";
    }
    system("$command > /dev/null 2>&1");
    @return = q12($lib, $dir);
    if (scalar(@return) == 2) {	# should be 2 only at this point, and they should be combined
      return (@return);
    }
  } else {
    if ($files[0] =~ /$ARGV[0]\.h5/) {
      $r1out = "$dir/$ARGV[0].h5/$ARGV[0]-combined_R1.fastq.gz";
      $r2out = "$dir/$ARGV[0].h5/$ARGV[0]-combined_R2.fastq.gz";
    } else {
      $r1out = "$dir/$ARGV[0]/$ARGV[0]-combined_R1.fastq.gz";
      $r2out = "$dir/$ARGV[0]/$ARGV[0]-combined_R2.fastq.gz";
    }
    @return = ($r1out, $r2out);
    open R1O, "| gzip > $r1out";
    open R2O, "| gzip > $r2out";
    foreach $i (0..$#files) {
      if ($files[$i] =~ /_R1\.fastq$/     ||
          $files[$i] =~ /_R1_\d+\.fastq$/ ||
          $files[$i] =~ /_R1_\d+_\w+\d+-[PS]E-R\d+.fastq$/ ||
          $files[$i] =~ /_R1_unaligned_\d+.fastq$/ ||
          $files[$i] =~ /_1\.fq$/
         ) {
        open R1I, $files[$i];
        while (<R1I>) {
          print R1O;
        }
        close R1I;
        unlink $files[$i] if $delete;
      } elsif ($files[$i] =~ /_R1\.fastq.gz$/     ||
               $files[$i] =~ /_R1_\d+\.fastq.gz$/ ||
               $files[$i] =~ /_R1_\d+_\w+\d+-[PS]E-R\d+.fastq.gz$/ ||
               $files[$i] =~ /_R1_unaligned_\d+.fastq.gz$/ ||
               $files[$i] =~ /_1\.fq\.gz$/
         ) {
        open R1I, "zcat $files[$i] |";
        while (<R1I>) {
          print R1O;
        }
        close R1I;
        unlink $files[$i] if $delete;
      } elsif ($files[$i] =~ /_R2\.fastq$/ ||
          $files[$i] =~ /_R2_\d+\.fastq$/  ||
          $files[$i] =~ /_R2_\d+_\w+\d+-[PS]E-R\d+.fastq$/ ||
          $files[$i] =~ /_R2_unaligned_\d+.fastq$/ ||
          $files[$i] =~ /_2\.fq$/
         ) {
        open R2I, $files[$i];
        while (<R2I>) {
          print R2O;
        }
        close R2I;
        unlink $files[$i] if $delete;
      } elsif ($files[$i] =~ /_R2\.fastq.gz$/     ||
               $files[$i] =~ /_R2_\d+\.fastq.gz$/ ||
               $files[$i] =~ /_R2_\d+_\w+\d+-[PS]E-R\d+.fastq.gz$/ ||
               $files[$i] =~ /_R2_unaligned_\d+.fastq.gz$/ ||
               $files[$i] =~ /_2\.fq\.gz$/
         ) {
        open R2I, "zcat $files[$i] |";
        while (<R2I>) {
          print R2O;
        }
        close R2I;
        unlink $files[$i] if $delete;
      }
    }
    close R1O;
    close R2O;
    if (-f $r1out && -f $r2out) {
      return @return;
    } else {
      die "Some problem, no $r1out or $r2out in combine...\n";
    }
  }
}

sub make_url {
  my ($run, $exp) = @_;
  return undef if !defined $run || !length $run || !defined $exp || !length $exp;
#  my $url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/";
#  $url .= substr($exp, 0, 3) . "/" . substr($exp, 0, 6) . "/" . $exp . "/";
#  $url .= $run . "/" . $run . ".sra";
# this changed sometime, edit Nov 2017 - now by run
  my $url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/";  $url .= substr($run, 0, 3) . "/" . substr($run, 0, 6) . "/" . $run . "/";
  $url .= $run . ".sra";
  return($url);
}
