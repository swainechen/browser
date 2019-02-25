#!/usr/bin/perl -w
#
# get info on a single library by run ID or tip number
#
use warnings;
use strict;
use File::Spec;
use File::Basename;
use GERMS;
use DBI;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $quiet = 0;	# if quiet, then return just an indication of completeness
my $show_help = 0;
GetOptions (
  'help' => \$show_help,
  'quiet' => \$quiet
);

if ($show_help || !defined $ARGV[0] || !length $ARGV[0]) {
  print "Usage: $0 <RunID> [ -quiet ]\n";
  print "Will print out information on this library\n";
  print "Can also give a TIP number but be careful about confusing with Run IDs that might be integers\n";
  print "With -quiet will just give an indication of completeness, explanation on STDERR:\n";
  print "  Fully complete: No output, exit code 0\n";
  print "  No TIP: TIP, exit code 1\n";
  print "  No fastq data: FASTQ, exit code 2\n";
  print "  No treefiles (depends on reference): TREEFILES, exit code 4\n";
  print "  No srst2: SRST2, exit code 8\n";
  print "  No assembly: ASSEMBLY, exit code 16\n";
  exit;
}

my $exit = 0;
my $USE_DB = 0;
my $DBH = GERMS::dbconnect();
my $sql;
my $sth;
my $tip;
my $run;
my ($i, $j, $k);
my @f;
my @data;
my $data;
my $output;
my $temp;
my $r;
my $keycol;
my $fullspecies;
my $shortspecies;
my $secondaryspecies;

# get basic info
if ($ARGV[0] =~ /^\d+$/) {
  # assume this is a TIP number
  $tip = $ARGV[0];
  $sql = "SELECT Run FROM Tips where TIP = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($tip);
  while (@data = $sth->fetchrow_array()) {
    $run = $data[0];
  }
} else {
  $run = $ARGV[0];
  $tip = GERMS::get_browser_tip($run, $DBH, $USE_DB);
}
if (!defined $run || !defined $tip || !length $run || $tip <= 0) {
  print STDERR "TIP\n";
  $exit += 1;
  exit($exit);
}

#$DBH->{FetchHashKeyName} = 'NAME_uc';

# Fastq
$r = get_all_from_table("Fastq", "TIP", $tip, "Technology");
$output->{FASTQ} = $r;
if (defined $output->{FASTQ}->{ILLUMINA}) {
  if (!defined $output->{FASTQ}->{ILLUMINA}->{Kraken} ||
      $output->{FASTQ}->{ILLUMINA}->{Kraken} eq "") {
    $exit += 2;
  }
} else {
  $exit += 2;
}

# Files - this has more keys (TIP, Type, MD5)
$sql = "SELECT * FROM Files WHERE TIP = ?";
$sth = $DBH->prepare($sql);
$sth->execute($tip);
while ($data = $sth->fetchrow_hashref) {
  push @{$output->{FILES}->{$data->{Type}}}, $data;
}
$sth->finish;

# SRST2 - Resistance and Genes
my $error_srst2 = 0;
$sql = "SELECT COUNT(*), SOURCE, SourceFileMD5 FROM Genes WHERE TIP = ? GROUP BY SOURCE, SourceFileMD5";
$sth = $DBH->prepare($sql);
$sth->execute($tip);
while (@data = $sth->fetchrow_array()) {
  $output->{SRST2}->{Genes}->{$data[1]}->{$data[2]} = $data[0];
}
if ($sth->rows != 1) {
  $error_srst2 = 1;
}
$sql = "SELECT COUNT(*), SOURCE, SourceFileMD5 FROM Resistance WHERE TIP = ? GROUP BY SOURCE, SourceFileMD5";
$sth = $DBH->prepare($sql);
$sth->execute($tip);
while (@data = $sth->fetchrow_array()) {
  $output->{SRST2}->{Resistance}->{$data[1]}->{$data[2]} = $data[0];
}
# sometimes we just don't have resistance, so don't count an error on this
#if ($sth->rows != 1) {
#  $error_srst2 = 1;
#}
$sql = "SELECT MLSTDatabase, MLST, SOURCE, SourceFileMD5 FROM MLST WHERE TIP = ? GROUP BY SOURCE, SourceFileMD5";
$sth = $DBH->prepare($sql);
$sth->execute($tip);
while (@data = $sth->fetchrow_array()) {
  $output->{SRST2}->{MLST}->{$data[2]}->{$data[3]} = "$data[1] ($data[0])";
}
if ($sth->rows != 1) {
  $error_srst2 = 1;
}
$exit += 8 if $error_srst2;

# Assembly
$r = get_all_from_table("Assembly", "TIP", $tip, "assembly");
$output->{ASSEMBLY} = $r;
if (scalar keys %{$output->{ASSEMBLY}} == 0) {
  $exit += 16;
} else {
  foreach $i (keys %{$output->{ASSEMBLY}}) {
    next if $i eq "assembly";
    if (!defined $output->{ASSEMBLY}->{$i}->{N50}) {
      $exit += 16;
    }
  }
}

# metadata
$r = get_all_from_table("Studies", "Run", $run, "Run");
$output->{META}->{Studies} = $r->{$run};
if (defined $output->{FASTQ}) {
  foreach $i (keys %{$output->{FASTQ}}) {
    if (defined $output->{FASTQ}->{$i}->{Kraken}) {
      $fullspecies = $output->{FASTQ}->{$i}->{Kraken};
      @f = split /\s+/, $fullspecies;
      if ($#f >= 1) {
        $shortspecies = uc(substr($f[0], 0, 1) . $f[1]);
      } else {
        $shortspecies = $fullspecies;
      }
      $secondaryspecies = $shortspecies;
      if ($shortspecies eq "SAGALACTIAE") {
        $secondaryspecies = "GBS";
      } elsif ($shortspecies eq "SPYOGENES") {
        $secondaryspecies = "GAS";
      }
      $sql = "SHOW TABLES";
      $sth = $DBH->prepare($sql);
      $sth->execute();
      while (@data = $sth->fetchrow_array()) {
        if (defined $data[0] && ( uc($data[0]) eq $shortspecies ||
                                  uc($data[0]) =~ /${shortspecies}_/ ||
                                  uc($data[0]) eq $secondaryspecies ||
                                  uc($data[0]) =~ /${secondaryspecies}_/ )) {
          $j = GERMS::get_columns($data[0], $DBH);
          if (defined $j->{run}) {
            $keycol = "run";
          } elsif (defined $j->{Run}) {
            $keycol = "Run";
          }
          $r = get_all_from_table($data[0], $keycol, $run, $keycol);
          $output->{META}->{$data[0]} = $r->{$run};
        }
      }
    }
  }
}

# Output
if (!$quiet) {
  print "Run: $run\n";
  print "Tip: $tip\n";
  print "FASTQ:\n";
  $r = $output->{FASTQ};
  foreach $i (keys %{$r}) {
    next if $i eq "Technology";
    print "  Technology: $i\n";
    print "  $r->{$i}->{NumReads} " if defined $r->{$i}->{NumReads};
    if (uc($r->{$i}->{Paired}) eq 'YES') {
      print "paired reads ";
    } else {
      print "unpaired reads ";
    }
    print "of length $r->{$i}->{ReadLength}\n" if defined $r->{$i}->{ReadLength};
    if (defined $r->{$i}->{Kraken}) {
      print "  Classified as $r->{$i}->{Kraken} ($r->{$i}->{KrakenPercent}\%)\n";
    }
    if (defined $r->{$i}->{Kraken2}) {
      print "  Secondary classification $r->{$i}->{Kraken2} ($r->{$i}->{Kraken2Percent}\%)\n";
    }
  }
  print "Files:\n";
  if (defined $output->{FILES}) {
    foreach $i (keys %{$output->{FILES}}) {
      $r = $output->{FILES}->{$i};
      foreach $j (0..$#$r) {
        if (-f $r->[$j]->{Filename}) {
          print "  $i Present: ";
        } else {
          print "  $i Not Present: ";
        }
        print "$r->[$j]->{Filename} ($r->[$j]->{MD5}; $r->[$j]->{DateStamp})\n";
      }
    }
  }
  print "SRST2:\n";
  if (defined $output->{SRST2}) {
    if (defined $output->{SRST2}->{MLST}) {
      $r = $output->{SRST2}->{MLST};
      foreach $j (sort keys %{$r}) {
        foreach $k (sort keys %{$r->{$j}}) {
          print "  MLST $r->{$j}->{$k} from $j ($k)\n";
        }
      }
    }
    foreach $i ("Genes", "Resistance") {
      $r = $output->{SRST2}->{$i};
      # $j is SOURCE
      # $k is SourceFileMD5
      foreach $j (sort keys %{$r}) {
        foreach $k (sort keys %{$r->{$j}}) {
          print "  $r->{$j}->{$k} calls for $i from $j ($k)\n";
        }
      }
    }
  }
  print "Assembly:\n";
  if (defined $output->{ASSEMBLY}) {
    $r = $output->{ASSEMBLY};
    foreach $i (sort keys %$r) {
      print "  $i assembly with $r->{$i}->{assembly_reads} reads, Total $r->{$i}->{total_length}, N50 $r->{$i}->{N50} ($r->{$i}->{SourceFileMD5})\n";
    }
  }
  print "Metadata:\n";
  if (defined $output->{META}) {
    if (defined $output->{META}->{Studies}) {
      print "  Studies table:\n";
      $r = $output->{META}->{Studies};
      foreach $i (sort keys %{$r}) {
        $r->{$i} = "" if !defined $r->{$i};
        print "    $i: $r->{$i}\n";
      }
    }
    foreach $i (keys %{$output->{META}}) {
      next if ($i eq "Studies");
      print "  $i table:\n";
      $r = $output->{META}->{$i};
      foreach $j (sort keys %$r) {
        $r->{$j} = "" if !defined $r->{$j};
        print "    $j: $r->{$j}\n";
      }
    }
  }
}
my $explain = $exit;
if ($explain >= 16) {
  $explain -= 16;
  print STDERR "Some problem with Assembly\n";
}
if ($explain >= 8) {
  $explain -= 8;
  print STDERR "Some problem with SRST2\n";
}
if ($explain >= 4) {
  $explain -= 4;
  print STDERR "Some problem with Treefiles\n";
}
if ($explain >= 2) {
  $explain -= 2;
  print STDERR "Some problem with Fastq\n";
}
if ($explain >= 1) {
  $explain -= 1;
  print STDERR "Some problem with Tip\n";
}
exit($exit);

sub get_all_from_table {
  my ($table, $keycolumn, $key, $fetchcolumn) = @_;
  my $sql = "SELECT * FROM $table WHERE $keycolumn = ?";
  my $sth = $DBH->prepare($sql);
  $sth->execute($key);
  my $data = $sth->fetchall_hashref($fetchcolumn);
  if ($sth->err || $sth->rows < 1) {
    return undef;
  } else {
    return($data);
  }
}
