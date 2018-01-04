#!/usr/bin/perl -w
#
# filter out certain regions for calling SNPs
# exclude phages, transposons, insertion sequences, rRNAs
# only use chromosomal sequences
# use gff file
#
use LWP::Simple;
my @blacklist = qw(
  transposase
  insertion.sequence
  ribosomal.RNA
);

my $exclude = ();
my $chromosomes = ();
while (<>) {
  chomp;
  next if /^$/; 
  next if /^#/;
  @f = split /\t/, $_;
  if ($f[8] =~ /;genome=(.*?);/) {
    $type->{$f[0]} = $1;
    $start->{$f[0]} = $f[3];
    $end->{$f[0]} = $f[4];
    next;
  }
  foreach $search (@blacklist) {
    if ($f[8] =~ /$search/i) {
      # track the reason
      @g = split /;/, $f[8];
      @r = ();
      foreach $j (0..$#g) {
        if ($g[$j] =~ /$search/i) {
          push @r, $g[$j];
        }
      }
      if (scalar @r) {
        $reason = join (";", @r);
      } else {
        $reason = "";
      }
      $exclude->{$f[0]}->{$f[3]}->{$f[4]} = $reason;
    }
  }
}

print "track name=exclude type=bedDetail description=\"Automatically generated list of excluded regions\"\n";
foreach $i (sort {$a cmp $b} keys %$type) {
  if ($type->{$i} eq "chromosome") {
    do_phast($i, $exclude);
    next if !defined $exclude->{$i};
    foreach $j (sort {$a <=> $b} keys %{$exclude->{$i}}) {
      foreach $k (sort {$a <=> $b} keys %{$exclude->{$i}->{$j}}) {
        print join ("\t", $i, $j, $k, $exclude->{$i}->{$j}->{$k}), "\n";
      }
    }
  } else {
    print join ("\t", $i, $start->{$i}, $end->{$i}, $type->{$i}), "\n";
  }
}

sub do_phast {
  my ($acc, $hash) = @_;
  # PHAST seems to not use the version number at the end
  my $url_acc = $acc;
  $url_acc =~ s/\.\d+$//;
#  my $PHAST_URL = 'http://phast.wishartlab.com/cgi-bin/change_to_html.cgi?num=__ACC__';
  my $PHAST_URL = 'http://phast.wishartlab.com/tmp/__ACC__/summary.txt';
  $PHAST_URL =~ s/__ACC__/$url_acc/;
  my $content = get($PHAST_URL);
  if (!length($content)) {
    print STDERR "Tried to get $PHAST_URL but had no content\n";
    print STDERR "Enter $url_acc into the form at http://phast.wishartlab.com to generate this\n";
    return;
  }
  my @content = split /\n/, $content;
  my $i;
  my @f;
  my $start;
  my $end;
  my $header = 1;
  foreach $i (0..$#content) {
    if ($content[$i] =~ /REGION\s*REGION_LENGTH/) {
      $header = 0;
      next;
    }
    if (!$header) {
      $content[$i] =~ s/\s+/\t/g;
      @f = split /\t/, $content[$i];
      next if !defined $f[5];
      if ($f[5] =~ /(\d+)-(\d+)/) {
        $start = $1;
        $end = $2;
        if ($f[3] =~ /intact/ || $f[3] =~ /incomplete/ || $f[3] =~ /questionable/) {
          $hash->{$acc}->{$start}->{$end} = "PHAST:$f[3],$f[14]";
        }
      }
    }
  }
}
