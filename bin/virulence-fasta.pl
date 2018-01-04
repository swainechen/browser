#!/usr/bin/perl -w
#
# take srst2.gz file
# parse out resistance genes fasta sequences
#
if (!defined $ARGV[0] || !-f $ARGV[0]) {
  print "Usage: $0 <srst2.gz file>\n";
  print "Parse out only those fasta sequences from ARGAnnot hits\n";
  print "Outputs to STDOUT.\n";
  exit;
}

if ($ARGV[0] =~ /\.gz$/) {
  open F, "zcat $ARGV[0] |";
} else {
  open F, $ARGV[0];
}

my $found_fasta = 0;
my @f;
my $hit;
while (<F>) {
  next if /^#/;
  next if /^$/;
  chomp;
  if (/^>/) {
    $found_fasta = 0;
    @f = split /__/, $_;
    if (defined $hit->{$f[1]} && $hit->{$f[1]}) {
      $found_fasta = 1;
      print $_, "\n";
    }
  } elsif (/\t/) {	# early part with tabular output
    @f = split /\t/, $_;
    if (defined $f[13] && $f[13] =~ /DB:ARGannot/) {
      $hit->{$f[2]} = 1;
    }
  } elsif ($found_fasta) {	# ought to be fasta sequence part
    print "$_\n";
  }
}
