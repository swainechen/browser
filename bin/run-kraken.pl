#!/usr/bin/perl -w #
# just a kraken run
# taken from GERMS-wgs.pl initially
#
use File::Spec;
use File::Temp;
use Cwd;
use GERMS;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

#my $KRAKEN_BIN = `which kraken`; chomp $KRAKEN_BIN;
#my $KRAKEN_REPORT = `which kraken-report`; chomp $KRAKEN_REPORT;
#my $KRAKEN_DB = "/usr/local/lib/kraken/minikraken_20171019_8GB";
my $KRAKEN2_BIN = `which kraken2`; chomp $KRAKEN2_BIN;
my $KRAKEN2_DB = "/usr/local/lib/Kraken2/minikraken2_v2_8GB_201904_UPDATE";
my $SEQTK = `which seqtk`; chomp $SEQTK;
my $q1 = "";
my $q2 = "";
my $show_help = 0;
my $VERBOSE = 0;
my $num_reads = 100000;
my $seed = 11;
my $current_dir;
my $classification;

GetOptions (
  'q1=s' => \$q1,
  'q2=s' => \$q2,
  'n=i' => \$num_reads,
  'seed=i' => \$seed,
  'help' => \$show_help,
  'verbose' => \$VERBOSE
);

if ($show_help) {
  &print_usage;
  exit;
}
if ($q1 eq "") {
  &print_usage;
  print "Error: No q1 file specified\n";
  exit;
}
if (!-f $q1) {
  &print_usage;
  print "Error: can't find q1 file $q1\n";
  exit;
}
if ($q2 ne "" && !-f $q2) {
  &print_usage;
  print "Error: can't find q2 file $q2\n";
  exit;
}

$q1 = File::Spec->rel2abs($q1);
if ($q2 ne "") {
  $q2 = File::Spec->rel2abs($q2);
}

$current_dir = getcwd;
$tempdir = File::Temp::tempdir( CLEANUP => 1 );
chdir($tempdir);

$q1 = GERMS::downsample($q1, $SEQTK, $num_reads, $seed, $tempdir);
if ($q2 ne "") {
  $q2 = GERMS::downsample($q2, $SEQTK, $num_reads, $seed, $tempdir);
}

if ($q1 ne "" && -f $q1) {
  $classification = run_kraken($q1, $q2);
  if (defined $classification->{Species}) {
    print join ("\t", $classification->{Species}->{Classification}, $classification->{Species}->{Percent}), "\n";
    if (defined $classification->{Species2}) {
      print join ("\t", $classification->{Species2}->{Classification}, $classification->{Species2}->{Percent}), "\n";
    }
  } elsif (defined $classification->{Genus}) {
    print join ("\t", $classification->{Genus}->{Classification}, $classification->{Genus}->{Percent}), "\n";
    if (defined $classification->{Genus2}) {
      print join ("\t", $classification->{Genus2}->{Classification}, $classification->{Genus2}->{Percent}), "\n";
    }
  }
}

chdir($current_dir);

sub print_usage {
  print "Usage: $0 -q1 <fastq1 file> [ -q2 <fastq2 file> ] [ -n <num reads> ] [ -seed <seed> ] [ -verbose ]\n";
  print "  Will run Kraken on the fastq files\n";
  print "  Default is to process $num_reads reads\n";
  print "  Seed sets the seed for seqtk to do downsampling for # of reads (default $seed)\n";
  print "  Default is to just give a final classification\n";
  print "  With -verbose provide all Kraken output\n";
  print "  Note: only q1 is required, q2 is optional\n";
  print "        q1 and q2 files should be the same type (fasta/fastq, gzipped or not)\n";
}

sub run_kraken {
  # these files should be downsampled already
  my ($q1, $q2) = @_;
  my $expected_file = "kraken2-out.report";
  my $command;
  my @f;
  my @g;
  my $i;
  my $output;
  my $parse = (); # keyed as level (G or S), line #, then Percent/Classification
  my $return;
  my @ft = GERMS::file_type($q1);

  $command = "$KRAKEN2_BIN --db $KRAKEN2_DB --report $expected_file";
  if ($q1 eq $q2 || $q2 eq "") {
    $command .= " $q1";
  } else {
    $command .= " $q1 $q2";
  }
  $command .= " 2>&1 > $expected_file";
  # for the kraken run, we're using the trick of 'command 2>&1 > output'
  # this puts stdout into the output file, and stderr shows up on stdout
  # so we can capture stderr with backticks
  # this is because we want the classifications into the output file
  # and the summary comes on stderr which we want to parse
  $output = `$command`;
  $output =~ s/\r/\n/g;
  if ($VERBOSE) {
    print $output;
    @f = split /\n/, $output;
    foreach $i (@f) {
      if ($i =~ /(\d+) sequences .* processed/) {
        print "Kraken2: $1 total sequences processed\n";
      } elsif ($i =~ /^\s*\d+ sequences classified/) {
        $i =~ s/^\s+//;
        print "Kraken2: $i\n";
      } elsif ($i =~ /^\s*\d+ sequences unclassified/) {
        $i =~ s/^\s+//;
        print "Kraken2: $i\n";
      }
    }
  }
  if (-f $expected_file) {
    open REPORT, $expected_file;
    @f = <REPORT>;
    close REPORT;
    foreach $i (0..$#f) {
      chomp $f[$i];
      $f[$i] =~ s/^\s+//;
      @g = split /\t/, $f[$i];
      if ($g[3] eq "G" || $g[3] eq "S") {
        $g[5] =~ s/^\s+//;
        $parse->{$g[3]}->{$i}->{Percent} = $g[0];
        $parse->{$g[3]}->{$i}->{Classification} = $g[5];
      }
    }
    unlink($expected_file);
    if (defined $parse->{S}) {
      @g = reverse sort 
             { $parse->{S}->{$a}->{Percent} <=> $parse->{S}->{$b}->{Percent} }
             keys %{$parse->{S}};
      if (scalar(@g)) {
        print "Kraken2 classification line: $f[$g[0]]\n" if $VERBOSE;
        $return->{Species}->{Classification} = $parse->{S}->{$g[0]}->{Classification};
        $return->{Species}->{Percent} = $parse->{S}->{$g[0]}->{Percent};
        if ($#g >= 1) {
          print "Kraken2 classification line: $f[$g[1]]\n" if $VERBOSE;
          $return->{Species2}->{Classification} = $parse->{S}->{$g[1]}->{Classification};
          $return->{Species2}->{Percent} = $parse->{S}->{$g[1]}->{Percent};
        }
      }
    }
    if (defined $parse->{G}) {
      @g = reverse sort
             { $parse->{G}->{$a}->{Percent} <=> $parse->{G}->{$b}->{Percent} }
             keys %{$parse->{G}};
      if (scalar(@g)) {
        print "Kraken2 classification line: $f[$g[0]]\n" if $VERBOSE;
        $return->{Genus}->{Classification} = $parse->{G}->{$g[0]}->{Classification};
        $return->{Genus}->{Percent} = $parse->{G}->{$g[0]}->{Percent};
        if ($#g >= 1) {
          print "Kraken2 classification line: $f[$g[1]]\n" if $VERBOSE;
          $return->{Genus2}->{Classification} = $parse->{G}->{$g[1]}->{Classification};
          $return->{Genus2}->{Percent} = $parse->{G}->{$g[1]}->{Percent};
        }
      }
    }
    return($return);
  } else {
    print "Couldn't find $expected_file file after initial Kraken2 run, skipping...\n" if $VERBOSE;
    return(undef);
  }
  return(undef);
}
