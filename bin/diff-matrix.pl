#!/usr/bin/perl -w
#
# take vcf file
# output full matrix of SNP, INDEL, or all differences
#
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $mode = "SNP";	# SNP, INDEL, or BOTH
my $min_cov = 10;
my $dotisindel = 0;	# if no coverage (./.) then default is just ignore
			# if $dotisindel = 1, then count this as an indel

GetOptions (
  "mode=s" => \$mode,
  "min_coverage=i" => \$min_cov,
  "dotisindel!" => \$dotisindel
);

my $matrix = ();
my $i;
my $j;
my @f;
my @samples;
$mode = uc($mode);

while (<>) {
  next if /^##/;
  chomp;
  if (/^#CHROM/) {
    @f = split /\t/, $_;
    @samples = @f[9..$#f];
    foreach $i (0..$#f-10) {
      foreach $j ($i+1..$#f-9) {
        $matrix->[$i]->[$j] = 0;
      }
    }
    next;
  }
  @f = split /\t/, $_;
  $ref = $f[3];
  @alt = split /,/, $f[4];
  @alleles = ($ref, @alt);
  # for multiallelic entries, can mix both SNPs and INDELs
  # so make sure use bcftools merge -m all to combine them
  # so we have all the info on each line
  # will have to process every line
#  if ($ref eq "." || length($ref) > 1) {
#    $class = "INDEL";
#  } else {
#    $class = "SNP";
#    foreach $i (0..$#alt) {
#      if ($alt[$i] eq "." || length($alt[$i]) > 1) {
#        $class = "INDEL";
#      }
#    }
#  }
#  if ($mode eq "BOTH" || $mode eq $class) {
  foreach $i (9..$#f) {
    @g = split /:/, $f[$i];
    if ($g[1] =~ /^\d+$/ && $g[1] >= $min_cov) {
      @g = split /\//, $f[$i];
      $f[$i] = $g[0];
    } else {
      $f[$i] = ".";
    }
  }
  foreach $i (0..$#f-10) {
    foreach $j ($i+1..$#f-9) {
      if ($f[$i+9] ne $f[$j+9]) {
        $class = "";
        if ($f[$i+9] ne "." && $f[$j+9] ne ".") {
          # we have some difference, and both are real alleles
          # check the lenghts of the alelles
          # also see if $ref or $alt is a "."
          if (length($alleles[$f[$i+9]]) == length($alleles[$f[$j+9]])) {
            # this should never happen that these alleles are the same here
            next if $alleles[$f[$i+9]] eq $alleles[$f[$j+9]];
            if ($alleles[$f[$i+9]] eq "." || $alleles[$f[$j+9]] eq ".") {
              $class = "INDEL";
            } else {
              $class = "SNP";
            }
          } else {
            $class = "INDEL";
          }
        } elsif ($dotisindel && ($f[$i+9] ne "." || $f[$j+9] ne ".")) {
          # This else clause counts no coverage as an indel
          # need one of them to be a real allele
          $class = "INDEL";
        }
        if ($class ne "" && ($class eq $mode || $mode eq "BOTH")) {
          $matrix->[$i]->[$j]++;
        }
      }
    }
  }
}

print "# Mode: $mode\n";
print join ("\t", "Sample", @samples), "\n";
foreach $i (0..$#samples) {
  @out = ($samples[$i]);
  foreach $j (0..$#samples) {
    if ($j < $i) {
      push @out, $matrix->[$j]->[$i];
    } elsif ($j > $i) {
      push @out, $matrix->[$i]->[$j];
    } else {
      push @out, 0;
    }
  }
  print join ("\t", @out), "\n";
}
