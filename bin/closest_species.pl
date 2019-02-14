#!/usr/bin/perl
#
use warnings;
use GERMS;
use Orgmap;
use File::Temp;
use Getopt::Long;
Getopt::Long::Configure("pass_through");

my $debug = 0;
my $num_test = 200;
my $initial = "";
GetOptions (
  'initial=s' => \$initial,
  'num_test=i' => \$num_test,
  'debug' => \$debug
);

if (!defined $ARGV[0] || !-f $ARGV[0]) {
  print "Usage: $0 [ -initial <string> ] [ -num_test <int> ] <fasta file>\n";
  print "Will give a list, one per line, of .fna files that can be used collectively\n";
  print "as a species reference for the input fasta file (for novel sequence calling,\n";
  print "etc.  The file first on the list is the one that is closest to the input file.\n";
  print "-initial specifies an initial species ID (such as from Kraken or an orgcode) so we can skip the initial blast against the reduced species database\n";
  print "-num_test is how many of the longest sequences to test. Intended to be something like N50 number or similar. Defaults to 200.\n";
  exit;
}

my $reduced = "/mnt/genomeDB/ncbi/genomes/Bacteria/bacteria-genus.fna";
my $tempdir;
if ($debug) {
  $tempdir = File::Temp::tempdir( CLEANUP => 0 );
  print STDERR "Intermediate files in $tempdir\n";
} else {
  $tempdir = File::Temp::tempdir( CLEANUP => 1 );
}
my $tempfile = "__closest_species.tmp";

if (!-f $reduced) {
  print "Can't find reduced fna file $reduced.  Please regenerate with make-reduced-fna.pl\n";
  die;
}

# the program itself...
print join ("\n", closest_species($ARGV[0], $num_test, $initial)), "\n";

##########################
# subroutines below here #
##########################
sub closest_species {
  my ($assembly, $num_test, $organism) = @_;
  my $best_orgcode = "";
  my @reference = (); 
  my @f;

  if ($organism ne "") {
    @reference = fna_files($organism);
    if (scalar @reference == 0) {
      @f = split /\s+/, $organism;
      if ($#f > 0) {
        @reference = fna_files(join (" ", @f[0,1]));
      } else {
        @reference = fna_files($f[0]);
      }
    }
  }
  if (scalar @reference == 0) {
    print STDERR "Couldn't find any reference sequences for initial organism $initial\n" if $debug && $organism ne "";
    $organism = "";
    print STDERR "Starting blast against reduced fna file\n" if $debug;
    ($organism, $best_orgcode) = blast_assembly($assembly, $num_test, $reduced);
    print STDERR "Return from reduced fna blast\n" if $debug;
    print STDERR "organism: $organism\nbest_orgcode: $best_orgcode\n" if $debug;
  }

  if ($organism ne "" && $best_orgcode ne "") {
    $organism =~ s/^\s*//;
    $organism =~ s/\s.*$//;	# hope genus is separated by space from species
    $db = create_db(fna_files($organism));	# make blast databases
    ($organism, $best_orgcode) = blast_assembly($assembly, $num_test, $db);
    @reference = fna_files($organism);	# don't make blast databases
    &read_orgmap($best_orgcode);
    foreach $i (0..$#reference) {
      if ($reference[$i] eq $Orgmap::fnafile) {
        splice (@reference, $i, 1);
        unshift (@reference, $Orgmap::fnafile);
      }
    }
  }

  return (@reference);
}

sub blast_assembly {
  my ($assembly, $num_test, $db) = @_;
  open ASSEMBLY, $assembly;
  my @a = <ASSEMBLY>;
  close ASSEMBLY;
  my $i = 0;
  my $contigs = GERMS::fasta2hash(@a);
  open BLAST, "| blastn -query - -db $db -outfmt '6 qseqid sseqid qlen length pident' -max_target_seqs 1 > $tempdir/$tempfile";
  foreach $key (reverse sort { length($contigs->{$a}) <=> length($contigs->{$b}) } keys %$contigs) {
    print BLAST ">$key\n$contigs->{$key}\n";
    $i++;
    last if $i >= $num_test;
  }
  close BLAST;

  open BLAST, "$tempdir/$tempfile";
  my $blast;
  my $votes;
  while ($blast = <BLAST>) {
    chomp $blast;
    @f = split /\t/, $blast;
    if (defined $votes->{$f[1]}) {
#      $votes->{$f[1]} += $f[4]/100 * $f[3]/$f[2];
      $votes->{$f[1]} += $f[4]/100 * $f[3];
    } else {
#      $votes->{$f[1]} = $f[4]/100 * $f[3]/$f[2];
      $votes->{$f[1]} = $f[4]/100 * $f[3];
    }
  }
  close BLAST;
  unlink ("$tempdir/$tempfile") if !$debug;

  my $best_org;
  foreach $key (reverse sort {$votes->{$a} <=> $votes->{$b}} keys %$votes) {
    $best_org = $key;
    last;
  }

  print "Best org: $best_org\n" if $debug;
  if (length($best_org)) {
    ($best_org, $orgcode) = get_species($best_org);
  } else {
    ($best_org, $orgcode) = ("", "");
  }
  return ($best_org, $orgcode);
}

sub get_species {
  my ($line) = @_;
  my @f = split /\|/, $line;
  my $acc;
  my $genus = "";
  my $species = "";
  if ($f[0] eq "gi" && $f[2] eq "ref") {
    $acc = $f[3];
    $acc =~ s/\.\d+$//;
  }
  print "Accession: $acc\n" if $debug;
  read_orgmap($acc);

  if (!defined $Orgmap::fnafile || !-f $Orgmap::fnafile) {
    return ("");
  } else {
    return ($Orgmap::orgname, $orgcode);
  }
}

sub fna_files {
  # we get an organism designation, maybe genus or genus/species or orgcode
  # return a list of all the .fna files in an array
  my ($organism) = @_;
  my @f;
  my @files = ();
  open ORG, "$Orgmap::LIBPATH/org-map";
  while ($line = <ORG>) {
    next if $line =~ /^$/;
    next if $line =~ /^#/;
    @f = split /\t/, $line;
    if ($f[0] eq $organism || $f[5] =~ /$organism/) {
      &read_orgmap($f[0]);
      push @files, $Orgmap::fnafile;
    }
  }
  close ORG;
  return @files;
}

sub create_db {
  # we get a list of fasta files
  # put them all together and make a blast database
  # return the name of the database
  my (@files) = @_;
  my $file;
  my $line;
  my $db = "$tempdir/__closest_species.fna";
  open DB, ">$db";
  foreach $file (@files) {
    if (-f $file) {
      open FNA, $file;
      while ($line = <FNA>) { print DB $line; }
      close FNA;
    }
  }
  close DB;
  system "makeblastdb -in $db -title $db -dbtype nucl 2>&1 > /dev/null";
  return ($db);
}
