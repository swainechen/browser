#!/usr/bin/perl -w
#
# run srst2 for generic bacteria - defaults to E. coli
# for MLST just keep the MLST line
# others take the full gene outputs
#
use warnings;
use strict;
use File::Basename;
use File::Temp;
use File::Spec;
use Cwd;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $name = "";
my $q1 = "";
my $q2 = "";
my @other_fasta = ();
my $mlst_fasta = "";
my $mlst_def = "";
my $mlst_delimiter = "_";
my $debug = 0;
my $samtools = "/usr/local/src/samtools-0.1.18/samtools";
my $bowtie2 = "/usr/local/src/bowtie2-2.2.9/bowtie2";
my $bowtie2_build = "/usr/local/src/bowtie2-2.2.9/bowtie2-build";
my $mlst_data = &set_mlst_data;
my $species = "";
my $showmeta = 0;	# just a help option to show metadata for a species
my $t;			# just a temp variable
my @to_show;

GetOptions (
  'name=s' => \$name,
  'q1=s' => \$q1,
  'q2=s' => \$q2,
  'mlst_db=s' => \$mlst_fasta,
  'mlst_def=s' => \$mlst_def,
  'mlst_delimiter=s' => \$mlst_delimiter,
  'species=s' => \$species,
  'other=s' => \@other_fasta,
  'showmeta' => \$showmeta,
  'debug!' => \$debug
);

if ($showmeta) {
  if (defined $species && defined($mlst_data->{$species})) {
    @to_show = ($species);
  } else {
    @to_show = sort keys %$mlst_data;
  }
  foreach $name (@to_show) {
    print "# Metadata for $name\n";
    foreach $t (sort keys %{$mlst_data->{$name}}) {
      if (ref($mlst_data->{$name}->{$t}) eq "ARRAY") {
        print join ("\t", $t, join (",", @{$mlst_data->{$name}->{$t}})), "\n";
      } else {
        print join ("\t", $t, $mlst_data->{$name}->{$t}), "\n";
      }
    }
  }
  exit;
}

if ($species ne "" && !defined($mlst_data->{$species})) {
  print STDERR "Invalid -species argument found.\n";
  print STDERR "Skipping MLST, only doing other files for E. coli\n";
  $species = "NONE";
  @other_fasta = @{$mlst_data->{'Ecoli'}->{OTHER}};
}

if (!length($name) || !length($q1) || !-f $q1) {
  if (!length($name)) {
    print STDERR "No -name argument found.\n";
  }
  if (!length($q1) || !-f $q1) {
    print STDERR "No valid -q1 argument found.\n";
  }
  print "Usage:\n";
  print "  $0 -q1 <fastq1> [ -q2 <fastq2> ] -name <sample name> -other <SRST2 database> -species <species>\n";
  print "  $0 -q1 <fastq1> [ -q2 <fastq2> ] -name <sample name> -other <SRST2 database> -mlst_db <mlst fasta> -mlst_def <mlst definitions> -mlst_delimiter <mlst delimiter>\n";
  print "Runs SRST2 for MLST and for any other databases\n";
  print "Files specified for fastq1 and fastq2 can be gzipped or not\n";
  print "Known values for species (show with -showmeta flag):\n";
  foreach $name (sort keys %$mlst_data) {
    print "  $name\n";
  }
  print "Other database is a fasta database formatted for SRST2 use\n";
  print "\n";
  print "You probably need to do something like\n";
  print "  source activate srst2\n";
  print "  $0 -q1 <q1> -q2 <q2> -name <name> -species <species> -other <other>\n";
  print "  source deactivate\n";
  exit;
}

if (length($species) && defined $mlst_data->{$species}) {
  $mlst_fasta = $mlst_data->{$species}->{FASTA};
  $mlst_def = $mlst_data->{$species}->{DEFINITIONS};
  $mlst_delimiter = $mlst_data->{$species}->{DELIMITER};
  if (!scalar @other_fasta && defined $mlst_data->{$species}->{OTHER}) {
    @other_fasta = @{$mlst_data->{$species}->{OTHER}};
  }
}


if (!File::Spec->file_name_is_absolute($q1)) {
  $q1 = File::Spec->rel2abs($q1);
}
if (defined($q2) && -f $q2) {
  if (!File::Spec->file_name_is_absolute($q2)) {
    $q2 = File::Spec->rel2abs($q2);
  }
} else {
  $q2 = "";
}

my $currentdir = getcwd;
my $tempdir;
my $filearg;
if ($debug) {
  $tempdir = File::Temp::tempdir( CLEANUP => 0 );
} else {
  $tempdir = File::Temp::tempdir( CLEANUP => 1 );
}
$ENV{'SRST2_SAMTOOLS'} = $samtools;
$ENV{'SRST2_BOWTIE2'} = $bowtie2;
$ENV{'SRST2_BOWTIE2_BUILD'} = $bowtie2_build;
chdir($tempdir);

# these filenames will follow the SRST2 naming system
my ($clean_q1, $clean_q2);
if ($q1 =~ /\.gz$/) {
  $clean_q1 = $name . "_1.fastq.gz";
} else {
  $clean_q1 = $name . "_1.fastq";
}
if (length $q2) {
  if ($q2 =~ /\.gz$/) {
    $clean_q2 = $name . "_2.fastq.gz";
  } else {
    $clean_q2 = $name . "_2.fastq";
  }
  $filearg = "--input_pe $clean_q1 $clean_q2";
} else {
  $filearg = "--input_se $clean_q1";
}

# we are using a local srst2 due to library differences - and we need to mangle the path for it so it gets the right bowtie2 also

`ln -s $q1 $tempdir/$clean_q1`;
if (length $q2) {
  `ln -s $q2 $tempdir/$clean_q2`;
}
if ($debug) {
  print STDERR "Path: $ENV{PATH}\n";
  print STDERR "bowtie2 ", `which bowtie2`, "\n";
  print STDERR "q1 $tempdir/$clean_q1\n";
  if (length $q2) {
    print STDERR "q2 $tempdir/$clean_q2\n";
  }
}

my $output;
my $mlst_final = "";
if ($species ne "NONE") {
  $output = `srst2 $filearg --log --output $name-mlst --mlst_db $mlst_fasta --mlst_definitions $mlst_def --mlst_delimiter "$mlst_delimiter" 2>&1`;
  
  if ($debug && defined $output) {
    print STDERR "srst2 $filearg --log --output $name-mlst --mlst_db $mlst_fasta --mlst_definitions $mlst_def --mlst_delimiter \"$mlst_delimiter\" 2>&1\n";
    print STDERR $output;
    }
  $mlst_final = `tail -n +2 -q $name-mlst__*__results.txt | sed -e 's/^/# /'`;
}

my $fullgenes_final = "";
my $fullgenes_fasta = "";
my $i;
my $result;
foreach $i (0..$#other_fasta) {
  $output = `srst2 $filearg --log --output $name-$i --gene_db $other_fasta[$i] --report_all_consensus 2>&1`;
  if ($debug && defined $output) {
    print STDERR "Command: srst2 $filearg --log --output $name-$i --gene_db $other_fasta[$i] --report_all_consensus 2>&1\n";
    print STDERR $output;
  }
  $result = "$name-$i" . "__fullgenes__*__results.txt";
  $fullgenes_final .= `tail -n +2 -q $result`;
  $result = "$name-$i.all_consensus_alleles.fasta";
  $fullgenes_fasta .= `cat $result`;
}

chdir($currentdir);
if (!length $q2) {
  # for some reason SRST2 puts a _1 if there's only one fastq file...
  $mlst_final =~ s/^# ${name}_1\t/# ${name}\t/;
  my @t = split /\n/, $fullgenes_final;
  foreach $i (0..$#t) {
    $t[$i] =~ s/^${name}_1\t/${name}\t/;
  }
  $fullgenes_final = join ("\n", @t);
  chomp $fullgenes_final;
  $fullgenes_final .= "\n";
}
print $mlst_final, $fullgenes_final, $fullgenes_fasta;

sub set_mlst_data {
  my $meta;
  my $base = "/usr/local/lib/SRST2";
  my $mlst_base = "$base/MLST";
  my $i;
  my $file;
  my @fasta;

  $meta->{"Abaumannii"}->{FASTA} = "$mlst_base/Acinetobacter_baumannii#1.fasta";
  $meta->{"Abaumannii"}->{DEFINITIONS} = "$mlst_base/abaumannii.txt";
  $meta->{"Abaumannii"}->{DELIMITER} = "_";

  $meta->{"Ecoli"}->{FASTA} = "$mlst_base/Escherichia_coli#1.fasta";
  $meta->{"Ecoli"}->{DEFINITIONS} = "$mlst_base/ecoli.txt";
  $meta->{"Ecoli"}->{DELIMITER} = "_";

  $meta->{"GBS"}->{FASTA} = "$mlst_base/Streptococcus_agalactiae.fasta";
  $meta->{"GBS"}->{DEFINITIONS} = "$mlst_base/sagalactiae.txt";
  $meta->{"GBS"}->{DELIMITER} = "_";

  $meta->{"Sagalactiae"}->{FASTA} = "$mlst_base/Streptococcus_agalactiae.fasta";
  $meta->{"Sagalactiae"}->{DEFINITIONS} = "$mlst_base/sagalactiae.txt";
  $meta->{"Sagalactiae"}->{DELIMITER} = "_";

  $meta->{"Kpneumoniae"}->{FASTA} = "$mlst_base/Klebsiella_pneumoniae.fasta";
  $meta->{"Kpneumoniae"}->{DEFINITIONS} = "$mlst_base/kpneumoniae.txt";
  $meta->{"Kpneumoniae"}->{DELIMITER} = "_";

  $meta->{"Ecloacae"}->{FASTA} = "$mlst_base/Enterobacter_cloacae.fasta";
  $meta->{"Ecloacae"}->{DEFINITIONS} = "$mlst_base/ecloacae.txt";
  $meta->{"Ecloacae"}->{DELIMITER} = "_";

  $meta->{"Ccoli"}->{FASTA} = "$mlst_base/Campylobacter_jejuni.fasta";
  $meta->{"Ccoli"}->{DEFINITIONS} = "$mlst_base/campylobacter.txt";
  $meta->{"Ccoli"}->{DELIMITER} = "_";

  $meta->{"Cjejuni"}->{FASTA} = "$mlst_base/Campylobacter_jejuni.fasta";
  $meta->{"Cjejuni"}->{DEFINITIONS} = "$mlst_base/campylobacter.txt";
  $meta->{"Cjejuni"}->{DELIMITER} = "_";

  $meta->{"Spneumoniae"}->{FASTA} = "$mlst_base/Streptococcus_pneumoniae.fasta";
  $meta->{"Spneumoniae"}->{DEFINITIONS} = "$mlst_base/spneumoniae.txt";
  $meta->{"Spneumoniae"}->{DELIMITER} = "_";

  $meta->{"Lmonocytogenes"}->{FASTA} = "$mlst_base/Listeria_monocytogenes.fasta";
  $meta->{"Lmonocytogenes"}->{DEFINITIONS} = "$mlst_base/lmonocytogenes.txt";
  $meta->{"Lmonocytogenes"}->{DELIMITER} = "_";

  $meta->{"GAS"}->{FASTA} = "$mlst_base/Streptococcus_pyogenes.fasta";
  $meta->{"GAS"}->{DEFINITIONS} = "$mlst_base/spyogenes.txt";
  $meta->{"GAS"}->{DELIMITER} = "_";

  $meta->{"Spyogenes"}->{FASTA} = "$mlst_base/Streptococcus_pyogenes.fasta";
  $meta->{"Spyogenes"}->{DEFINITIONS} = "$mlst_base/spyogenes.txt";
  $meta->{"Spyogenes"}->{DELIMITER} = "_";

  $meta->{"Senterica"}->{FASTA} = "$mlst_base/Salmonella_enterica.fasta";
  $meta->{"Senterica"}->{DEFINITIONS} = "$mlst_base/senterica.txt";
  $meta->{"Senterica"}->{DELIMITER} = "_";

  $meta->{"Efaecalis"}->{FASTA} = "$mlst_base/Enterococcus_faecalis.fasta";
  $meta->{"Efaecalis"}->{DEFINITIONS} = "$mlst_base/efaecalis.txt";
  $meta->{"Efaecalis"}->{DELIMITER} = "_";

  $meta->{"Paeruginosa"}->{FASTA} = "$mlst_base/Pseudomonas_aeruginosa.fasta";
  $meta->{"Paeruginosa"}->{DEFINITIONS} = "$mlst_base/paeruginosa.txt";
  $meta->{"Paeruginosa"}->{DELIMITER} = "_";

  $meta->{"Saureus"}->{FASTA} = "$mlst_base/Staphylococcus_aureus.fasta";
  $meta->{"Saureus"}->{DEFINITIONS} = "$mlst_base/saureus.txt";
  $meta->{"Saureus"}->{DELIMITER} = "_";

  foreach $i (keys %$meta) {
    if (!-f $meta->{$i}->{FASTA} || !-f $meta->{$i}->{DEFINITIONS}) {
      delete($meta->{$i});
    }
  }

  foreach $i (keys %$meta) {
    if (-d "$base/$i") {
      @fasta = ();
      opendir O, "$base/$i";
      while ($file = readdir O) {
        if ($file =~ /.*-combined-.*\.fasta$/) {
          push @fasta, $file;
        }
      }
      @fasta = sort {$a cmp $b} @fasta;
      if (scalar @fasta) {
        $meta->{$i}->{OTHER} = [ "$base/$i/$fasta[$#fasta]" ];
      }
    }
  }

  return $meta;
}
