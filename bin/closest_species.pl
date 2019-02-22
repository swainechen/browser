#!/usr/bin/perl
#
use warnings;
use GERMS;
use Orgmap;
use File::Temp;
use File::Fetch;
use File::Basename;
use FindBin;
use Getopt::Long;
Getopt::Long::Configure("pass_through");

# this file comes from http://www.ncbi.nlm.nih.gov/genome/browse/#
# currently (14 Feb 2019) the fields are (with examples, quotes removed):
# 0 - Organism Name: Corynebacterium jeikeium K411
# 1 - Organism Groups: Bacteria;Terrabacteria group;Actinobacteria
# 2 - Strain: K411 = NCTC 11915
# 3 - BioSample: SAMEA3283089
# 4 - BioProject: PRJNA13967
# 5 - Assembly: GCA_000006605.1
# 6 - Level: Complete
# 7 - Size(Mb): 2.47682
# 8 - GC%: 61.3561
# 9 - Replicons: chromosome:NC_007164.1/CR931997.1; plasmid pKW4:NC_003080.1/AF401314.1
# 10 - WGS:
# 11 - Scaffolds: 2
# 12 - CDS: 2088
# 13 - Release Date: 2005-06-27T00:00:00Z
# 14 - GenBank FTP: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/605/GCA_000006605.1_ASM660v1
# 15 - RefSeq FTP: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/605/GCF_000006605.1_ASM660v1
my $PROKS_CSV = $FindBin::Bin . "/../lib/ncbi/prokaryotes.csv";
# number of genomes to get on the fly from the above file if needed
my $NEW_DOWNLOAD = 10;
my $DOWNLOAD = "";
my $debug = 0;
my $num_test = 200;
my $initial = "";
GetOptions (
  'initial=s' => \$initial,
  'num_test=i' => \$num_test,
  'download=s' => \$DOWNLOAD,	# where to leave downloaded genomes
  'debug' => \$debug
);

if (!defined $ARGV[0] || !-f $ARGV[0]) {
  print "Usage: $0 [ -initial <string> ] [ -num_test <int> ] [ -download <path> ] <fasta file>\n";
  print "Will give a list, one per line, of .fna files that can be used collectively\n";
  print "as a species reference for the input fasta file (for novel sequence calling,\n";
  print "etc.  The file first on the list is the one that is closest to the input file.\n";
  print "-initial specifies an initial species ID (such as from Kraken or an orgcode) so we can skip the initial blast against the reduced species database\n";
  print "-num_test is how many of the longest sequences to test. Intended to be something like N50 number or similar. Defaults to 200.\n";
  print "-download specifies a download directory for new genomes. If specified, these genomes will be left there and not cleaned up after program execution.\n";
  exit;
}

if ($initial eq "") {
  # we can't really do anything at this point
  exit;
}
my $tempdir;
if ($debug) {
  $tempdir = File::Temp::tempdir( CLEANUP => 0 );
  print STDERR "Intermediate files in $tempdir\n";
} else {
  $tempdir = File::Temp::tempdir( CLEANUP => 1 );
}
my $tempfile = "__closest_species.tmp";

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
    exit;
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
  # The database of fna files used to be at ftp://ftp.ncbi.nih.gov/genomes/Bacteria/
  # Now seems to be at ftp://ftp.ncbi.nih.gov/genomes/refseq/bacteria but follows the GCF/GCA designation now
  # In any case, it's huge...only have a few downloaded
  # If we don't have anything downloaded here then go to prokaryotes.csv file to fetch some on the fly
  my ($organism) = @_;
  my @f;
  my @g;
  my @files = ();
  my $ff;
  my $url;
  my $tempfile;
  my $output_dir;
  if ($DOWNLOAD eq "") {
    $output_dir = $tempdir;
  } else {
    $output_dir = $DOWNLOAD;
  }
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
  if (!scalar @files && -f $PROKS_CSV) {
    open F, $PROKS_CSV;
    while (<F>) {
      chomp;
      next if /^#/;
      @f = split /,/, $_;
      foreach $i (0..$#f) {
        $f[$i] =~ s/^"//;
        $f[$i] =~ s/"$//;
      }
      next if $f[6] ne "Complete";
      if ($f[0] =~ /^$organism/) {
        $tempfile = "";
        $url = $f[15];
        $url .= "/" . File::Basename::basename($url) . "_genomic.fna.gz";
        $ff = File::Fetch->new( uri => $url );
        $ff->fetch( to => $output_dir );
        if (-f $ff->output_file) {
          $tempfile = File::Spec->rel2abs($ff->output_file);
        } elsif (-f "$output_dir/" . $ff->output_file) {
          $tempfile = File::Spec->rel2abs("$output_dir/" . $ff->output_file);
        }
        if ($tempfile ne "") {
          if ($tempfile =~ /\.gz$/) {
            system "gunzip $tempfile";
            $tempfile =~ s/\.gz$//;
            push @files, $tempfile;
            if (scalar(@files) >= $NEW_DOWNLOAD) {
              return @files;
            }
          }
        }
      }
    }
  }
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
