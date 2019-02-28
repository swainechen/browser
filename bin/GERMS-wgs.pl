#!/usr/bin/perl
#
# wrapper to run the whole GERMS WGS pipeline
# - assembly (SOAPdenovo or velvetoptimizer)
# - scaffolding (SOAPdenovo)
# - reference selection (BLAST)
# - SNP calling (MUMmer)
# - SNP annotation (GERMS::ns, taken from Swaine - needs upgrade)
# - indel calling ()
# - new sequences (select_novel from Niranjan, based on MUMmer(?))
# - gene prediction (prodigal)
#
use warnings;
use strict;
use GERMS;
use Getopt::Long;
use LWP::Simple;
use File::Temp;
use File::Basename;
use File::Spec;
use File::Copy;
use File::Path;
use Bio::SeqIO;
use Cwd;
use JSON;

# programs we need
my $local_bin_dir = "/usr/local/bin";
my $bin_dir = "/usr/bin";
my $programs;
my $original_commandline = join (" ", $0, @ARGV);
my @arguments = @ARGV;
$programs->{SOAPdenovo} = "$bin_dir/soapdenovo-63mer";
$programs->{GapCloser} = "$local_bin_dir/GapCloser";
$programs->{OPERA} = "$local_bin_dir/opera";
$programs->{OPERA_preprocess} = "$local_bin_dir/preprocess_reads.pl";
$programs->{nucmer} = "$bin_dir/nucmer";
$programs->{showcoords} = "$bin_dir/show-coords";
$programs->{showsnps} = "$bin_dir/show-snps";
$programs->{mummerplot} = "$bin_dir/mummerplot";
$programs->{prodigal} = "$bin_dir/prodigal";
$programs->{selectnovel} = "$local_bin_dir/select_novel.pl";
#$programs->{qseq2fastq} = "$bin_dir/illumina2fastq";
$programs->{fastqc} = "$bin_dir/fastqc";
$programs->{kraken} = "$bin_dir/kraken";
$programs->{krakenreport} = "$bin_dir/kraken-report";
$programs->{closest_species} = "$local_bin_dir/closest_species.pl";
#$programs->{elm_client} = "/mnt/software/lib/elmClient.cml.jar";
$programs->{VelvetShuffle} = "$local_bin_dir/shuffleSequences_fastq.pl";
$programs->{VelvetOptimiser} = "$local_bin_dir/VelvetOptimiser.pl";
$programs->{velvetg} = "$local_bin_dir/velvetg";
$programs->{FinIS} = "$local_bin_dir/finis";
$programs->{spades} = "$local_bin_dir/spades.py";
$programs->{prokka} = "$local_bin_dir/prokka";
$programs->{seqtk} = "$local_bin_dir/seqtk";
foreach my $ex (keys %$programs) {
  if (!-f $programs->{$ex}) {
    print "Can't find required program/library $ex, expected at $programs->{$ex}.  Exiting...\n";
    die;
  }
}
my $mux_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__MUX__/solexaRun/json";
my $library_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__LIB__/expand/json";

# parameters
my $tempdir;
my $final_sample_name = "final_output";
my $libname;
my $sampleID = "";
my $q1;
my $q2;
my $cleanup = 1;
my $temp_base = "/tmp";
my $output_dir = "/tmp";
my $cluster = 0;
my $old = 0;	# 0 means CASAVA 182 (new); 1 means older version
my @not_wanted = (
  "Arc",
  "ContigIndex",
  "edge",
  "gapFilled.fill",
  "gapSeq",
  "kmerFreq",
  "links",
  "markOnEdge",
  "path",
  "peGrads",
  "preArc",
  "preGraphBasic",
  "readInGap",
  "readOnContig",
  "scaf_gap",
  "updated.edge",
  "vertex",
  "read-on-contig",
  "scafSeq.amb",
  "scafSeq.ann",
  "scafSeq.pac",
  "scafSeq.bwt",
  "scafSeq.sa",
  "contigs.txt.amb",
  "contigs.txt.ann",
  "contigs.txt.pac",
  "contigs.txt.bwt",
  "contigs.txt.sa",
  "contigs.txt.1.ebwt",
  "contigs.txt.2.ebwt",
  "contigs.txt.3.ebwt",
  "contigs.txt.4.ebwt",
  "contigs.txt.rev.1.ebwt",
  "contigs.txt.rev.2.ebwt",
  "spades.fasta.amb",
  "spades.fasta.ann",
  "spades.fasta.bwt",
  "spades.fasta.pac",
  "spades.fasta.sa",
  "mummerplot.filter",
  "mummerplot.fplot",
  "mummerplot.rplot",
  "mummerplot.gp"
);

# These option hashes are going to be global...
my $default_options = {};
$default_options->{VELVETOPTIMIZER}->{s} = 0;	# 0 is auto
$default_options->{VELVETOPTIMIZER}->{e} = 99;	# fallback if auto fails
$default_options->{SOAP}->{max_rd_len} = 188;
$default_options->{SOAP}->{avg_ins} = 300;
$default_options->{SOAP}->{reverse_seq} = 0;
$default_options->{SOAP}->{asm_flags} = 3;
$default_options->{SOAP}->{rd_len_cutoff} = 100;
$default_options->{SOAP}->{rank} = 1;
$default_options->{SOAP}->{k} = 41;
$default_options->{SOAP}->{n_cpu} = 16;
$default_options->{SOAP}->{overlap} = 31;
$default_options->{FASTQC}->{dofastqc} = 0;
$default_options->{threads} = 16;
$default_options->{ASSEMBLER} = "velvetoptimizer";
$default_options->{SCAFFOLDER} = "opera";
# As of OPERA_LG_v2.0.6 there are some syntax errors with bwa as mapper
$default_options->{OPERA}->{mapper} = "bowtie";
$default_options->{contig_min} = 500;
$default_options->{ANNOTATOR} = "prokka";	# versus "old"
$default_options->{MAXREADS} = 5000000;
$default_options->{MEMORYINGB} = 16;
$default_options->{SEED} = 11;
$default_options->{KRAKEN}->{reads} = 100000;
$default_options->{KRAKEN}->{DB} = "/mnt/genomeDB/misc/softwareDB/kraken/minikraken_20141208";
$default_options->{KRAKEN}->{classification} = "NONE";

my $user_options = {};

# variables
my $coords;
my $rank;
my $novel;
my @reference;
my $reference;
my $delta;
my $snps;
my $genes;
my $contigs;
my $ordering;
my $assembly;
my @results;
my $index_dir;
my $current_dir;
my $lane;
my $pttfile;
my $ptt_ref;
my $best_reference_sequence;
my $mux_data;
my $mux;
my $plex;
my $global_threads = 0;
my $kraken_classification = "";
my $preferred_reference = "";
my $tempparent;

GetOptions (
  'sample=s' => \$sampleID,
  'tempdir=s' => \$temp_base,
  'q1=s' => \$q1,
  'q2=s' => \$q2,
  'name=s' => \$final_sample_name,
  't=i' => \$global_threads,
  'insert=i' => \$user_options->{SOAP}->{avg_ins},
  'soap=s%' => \$user_options->{SOAP},
  'velvetoptimizer=s%' => \$user_options->{VELVETOPTIMIZER},
  'fastqc=s%' => \$user_options->{FASTQC},
  'qc!' => \$user_options->{FASTQC}->{dofastqc},
  'contig_min=i' => \$user_options->{contig_min},
  'index_dir=s' => \$index_dir,
  'lane=i' => \$lane,
  'mux=s' => \$mux,
  'libname=s' => \$libname,
  'cleanup!' => \$cleanup,
  'cluster!' => \$cluster,
  'output_dir=s' => \$output_dir,
  'old' => \$old,
  'reference=s' => \$user_options->{REFERENCE},
  'assembler=s' => \$user_options->{ASSEMBLER},
  'scaffolder=s' => \$user_options->{SCAFFOLDER},
  'annotator=s' => \$user_options->{ANNOTATOR},
  'maxreads=i' => \$user_options->{MAXREADS},
  'memory=i' => \$user_options->{MEMORYINGB},
  'seed=i' => \$user_options->{SEED},
  'krakenreads=i' => \$user_options->{KRAKEN}->{reads},
  'species=s' => \$user_options->{KRAKEN}->{classification}
);

# postprocess some options

if ($global_threads) {
  $user_options->{threads} = $global_threads;
}

# make sure we have absolute paths
if (!File::Spec->file_name_is_absolute($temp_base)) {
  $temp_base = File::Spec->rel2abs($temp_base);
}
if (!File::Spec->file_name_is_absolute($output_dir)) {
  $output_dir = File::Spec->rel2abs($output_dir);
}
if (!-d $output_dir) {
  File::Path::make_path($output_dir) || die "Can't find or make output directory $output_dir\n";
}

# figure out how we're getting the files and info
# if we have MUX that's all we need and we go get the info and run it all
# if we see the -cluster switch we just get info and print out commands to run everything
# these commands specify: -name -index_dir -lane -insert -nocluster -output_dir
# to run things manually, we can also use this info: -q1 -q2 -name -insert -output_dir
# if running manually another common option would be -reference

if (defined $mux && length $mux) {
  if ($old) {
    $mux_data = GERMS::info_from_mux($mux);
  } else {
    $mux_data = GERMS::info_from_mux_182($mux);
  }
} elsif (defined $index_dir && -d $index_dir &&
         defined $lane &&
         defined $final_sample_name &&
         length $final_sample_name) {
  undef $mux_data;
} elsif (!defined $q1 || !-f $q1 || (defined $q2 && !-f $q2) || !length($final_sample_name)) {
    print "Can't find files or have no sample name, exiting...\n";
    &print_usage;
    exit;
}


####################
# Do all the stuff #
####################

# easier with absolute filenames and paths, we'll check $q1 and $q2 later
$current_dir = getcwd;
$output_dir = File::Spec->rel2abs($output_dir);

if (defined $mux_data->{RUNID}) {
  foreach $plex (keys %$mux_data) {
    next if $plex !~ /^\d+$/;

    if((defined $libname) && (uc($libname) ne $mux_data->{$plex}->{LIBRARY})){
	next;
    }
    # set up all the information we need
#    $final_sample_name = join ("_", $mux_data->{$plex}->{LIBRARY}, $mux_data->{$plex}->{BARCODE}, $mux_data->{$plex}->{DESCRIPTION}, $mux_data->{$plex}->{INDEXDIR});
    $final_sample_name = join ("_", $mux_data->{$plex}->{LIBRARY}, $mux_data->{$plex}->{BARCODE}, $mux_data->{$plex}->{DESCRIPTION});
    $final_sample_name = clean_final_name($final_sample_name);

    $index_dir = $mux_data->{$plex}->{FULLDIR};
    $lane = $mux_data->{LANEID};
    $default_options->{SOAP}->{avg_ins} = $mux_data->{INSERTLENGTH};


    if ($cluster) {
      foreach my $i (reverse (0..$#arguments)) {
        if ($arguments[$i] eq "-mux" || $arguments[$i] eq "-libname") {
          splice @arguments, $i, 2;
        }
      }
      if ($old) {
        print "'bsub $0 -old -name $final_sample_name -index_dir $index_dir -lane $lane -insert $default_options->{SOAP}->{avg_ins} -nocluster ", join (" ", @arguments), "'\n";
      } else {
        print "bsub '$0 -name $final_sample_name -index_dir $index_dir -lane $lane -insert $default_options->{SOAP}->{avg_ins} -nocluster ", join (" ", @arguments), "'\n";
      }
    } else {
      # set up files and directories we need
      $user_options->{tempparent} = File::Temp::tempdir('WGS-XXXXXX', DIR => $temp_base, CLEANUP => $cleanup );
      $user_options->{tempdir} = "$user_options->{tempparent}/$final_sample_name";
      $tempdir = set_option("tempdir");
      mkdir($tempdir);
      open LOG, ">$tempdir/$final_sample_name-log.txt";
      open SHORTLOG, ">$tempdir/$final_sample_name-summary.txt";
      &shortlog("Your original command line was $original_commandline\n");

      if ($old) {
        ($q1, $q2) = GERMS::get_data_bydir($index_dir, $lane);
      } else {
        ($q1, $q2) = GERMS::get_data_bydir_182($index_dir, $tempdir, $final_sample_name, "nogz");
      }

      # Since Opera can't handle *.gz files, we need to unzip before further processing.
      if (($q1 =~ /\.gz$/) && ($q2 =~ /\.gz$/)) {
	  print `gunzip $q1`;
	  print `gunzip $q2`;
	  #remove .gz from the file extension : from here unzipped files are used! 
	  $q1 =~ s/\.[^.]*$//;
	  $q2 =~ s/\.[^.]*$//;
      }

      &shortlog(do_single_library($final_sample_name, $q1, $q2));
      &post_process($final_sample_name);
    }
  }
} else {
  $final_sample_name = clean_final_name($final_sample_name);

  # set up files and directories we need
  $user_options->{tempparent} = File::Temp::tempdir('WGS-XXXXXX', DIR => $temp_base, CLEANUP => $cleanup );
  $user_options->{tempdir} = "$user_options->{tempparent}/$final_sample_name";
  $tempdir = set_option("tempdir");
  mkdir($tempdir);
  open LOG, ">$tempdir/$final_sample_name-log.txt";
  open SHORTLOG, ">$tempdir/$final_sample_name-summary.txt";
  &shortlog("Your original command line was $original_commandline\n");

  # we should have defined index_dir, lane, final_sample_name from command line
  # also insert_length, maybe default or from command line
  if (defined $index_dir && defined $lane) {
    if ($old) {
      ($q1, $q2) = GERMS::get_data_bydir($index_dir, $lane);
    } else {
      ($q1, $q2) = GERMS::get_data_bydir_182($index_dir,$tempdir,$final_sample_name, "nogz");
    }
  }

  &shortlog(do_single_library($final_sample_name, $q1, $q2));
  &post_process($final_sample_name);
}
# return to our original place so we can clean up temp directories if needed
chdir($current_dir);


#########################
# Subroutines down here #
#########################

sub clean_final_name {
  my ($name) = @_;
  if (defined $name && length $name) {
    $name =~ s/\s+/_/g;
    $name =~ s/\//_/g;
    $name =~ s/\\/_/g;
    $name =~ s/\?/_/g;
    $name =~ s/\!/_/g;
    $name =~ s/\(/_/g;
    $name =~ s/\)/_/g;
    $name =~ s/'/_/g;
    $name =~ s/"/_/g;
  } else {
    $name = "";
  }
  return ($name);
}

sub config_SOAP {
  my ($q1, $q2) = @_;
  # put this in a separate SOAP directory
  my $tempdir = set_option("tempdir");
  my $config_file = "$final_sample_name.SOAP.cfg";
  my $config_string = "";
  if (!-f "$tempdir/SOAP/$config_file") {
    $config_string .= "max_rd_len=" . set_option("SOAP", "max_rd_len") . "\n";
    $config_string .= "[LIB]\n";
    $config_string .= "avg_ins=" . set_option("SOAP", "avg_ins") . "\n";
    $config_string .= "reverse_seq=" . set_option("SOAP", "reverse_seq") . "\n";
    $config_string .= "asm_flags=" . set_option("SOAP", "asm_flags") . "\n";
    $config_string .= "rd_len_cutoff=" . set_option("SOAP", "rd_len_cutoff") . "\n";
    $config_string .= "rank=" . set_option("SOAP", "rank") . "\n";
    $config_string .= "q1=$q1\n";
    if ($q1 ne $q2) {
      $config_string .= "q2=$q2\n";
    }
    open SOAP_CONFIG, ">$tempdir/SOAP/$config_file";
    print SOAP_CONFIG $config_string;
    close SOAP_CONFIG;
  }
  return ("$tempdir/SOAP/$config_file");
}

sub assemble_SOAPDENOVO {
  # we need fastq sequences
  my ($q1, $q2) = @_;
  # make the config file
  my $tempdir = set_option("tempdir");
  mkdir("$tempdir/SOAP");
  chdir("$tempdir/SOAP");
  my $expected_file = "$tempdir/$final_sample_name.scafSeq";
  my $SOAP_output_contigs = "$tempdir/SOAP/$final_sample_name.scafSeq";
  my $config_file = config_SOAP($q1, $q2);
  my $command = "$programs->{SOAPdenovo} all";
  $command .= " -K " . set_option("SOAP", "k");
  $command .= " -R -d 1 -p " . set_option("SOAP", "n_cpu");
  $command .= " -s $config_file -o $tempdir/SOAP/$final_sample_name";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  chdir($tempdir);
  if (-f $SOAP_output_contigs) {
    File::Copy::copy($SOAP_output_contigs, $expected_file);
    return ($expected_file);
  } else {
    &log("Couldn't find $SOAP_output_contigs after ".(caller(0))[3].", SOAPdenovo run");
    return (0);
  }
}

sub assemble_SPAdes {
  # we need fastq sequences
  my ($q1, $q2) = @_;
  my $tempdir = set_option("tempdir");
  my $expected_file = "$tempdir/$final_sample_name.spades.fasta";
  my $spades_scaffolds = "$tempdir/spades_working/scaffolds.fasta";
  my $command = "$programs->{spades}";
  my @out;
  my $in;
  my $head;
  my $s;
  if ($q1 eq $q2) {
    # single-end
    $command .= " --s1 $q1";
  } else {
    $command .= " --pe1-1 $q1 --pe1-2 $q2";
  }
  $command .= " -t " . set_option("threads");
  $command .= " -m " . set_option("MEMORYINGB");
  $command .= " --cov-cutoff auto --careful -o $tempdir/spades_working";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  chdir($tempdir);
  if (-f $spades_scaffolds) {
    # we have to do our own length filtering for SPAdes assembly
    open SPADES_SCAF, $spades_scaffolds;
    open SPADES_FILT, ">$expected_file";
    @out = ();
    $head = "";
    $s = "";
    while ($in = <SPADES_SCAF>) {
      chomp $in;
      if ($in =~ /^>/) {
        if ($head ne "" && length($s) > set_option("contig_min")) {
          print SPADES_FILT "$head\n";
          print SPADES_FILT "$s\n";
        }
        $head = $in;
        $s = "";
      } else {
        $s .= $in;
      }
    }
    if ($head ne "" && length($s) > set_option("contig_min")) {
      print SPADES_FILT "$head\n";
      print SPADES_FILT "$s\n";
    }
    close SPADES_FILT;
    close SPADES_SCAF;
    return($expected_file);
  } else {
    &log("Couldn't find $spades_scaffolds after ".(caller(0))[3].", SPAdes run");
    return (0);
  }
}

sub assemble_velvet {
  # use velvet optimiser
  # options should only be start and end kmer to use
  my ($q1, $q2) = @_;
  my $tempdir = set_option("tempdir");
  my $kstart = set_option("VELVETOPTIMIZER", "s");
  my $kend = set_option("VELVETOPTIMIZER", "e");
  my $readlength = `head -n 4 $q1 | grep -A1 '^\@' | tail -n 1 | wc | awk '{printf "\%.0f", \$3 - 1}'`;
  my $best_middle = 0;
  my $velvetoptimizer_output_dir = "$tempdir/velvetoptimiser";
  my $velvetlog = "";
  my $vlog;
  my $velvetg_options;
  my $velveth_k;
  my $single_fastq;
  my $file_spec;
  my $shuffle;
  my $expected_file;
  my $command;
  my $velvet_output_contigs;
  my $output;
  # we need to make a single file with all the reads for velvet
  if ($q1 eq $q2) {
    $single_fastq = 1;
    $file_spec = "'-short -fastq $q1'";
  } else {
    $single_fastq = 0;
    $shuffle = "$tempdir/$final_sample_name.shuffle-fastq.txt";
    $expected_file = $shuffle;
    $command = "$programs->{VelvetShuffle} $q1 $q2 $shuffle";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (!-f $expected_file) {
      &log("Couldn't find $expected_file after ".(caller(0))[3].", velvet shuffle run");
      return (0);
    }
    $file_spec = "'-shortPaired -fastq $shuffle'";
  }
  # VelvetOptimizer with multiple threads tries to use OpenMP which isn't available on the SMPs.  Force threads to be 1 (-t option) for now
  if ($kstart == 0) {
    # try do to an optimization
    # just guess the kmer should be around 70% of the read length - and odd
    # try to adjust up or down appropriately
    # do 10 kmer values at a time
    # max kmer should be 255, but let's say reasonable starting point no matter
    # what the read length is 99, and a reasonable ending point is 31
    # when we find the right one we just set kstart and kend to the optimal
    # and run the original version
    $kstart = int(($readlength * 0.7) / 2) * 2 + 1 - 10;
    $kstart = 99 if $kstart > 99;
    &shortlog("Beginning auto optimization of velvet kmer");
    while (!$best_middle) {
      $kend = $kstart + 20;
      $kend = 255 if $kend > 255;
      $command = "$programs->{VelvetOptimiser} -s $kstart -e $kend -d $tempdir/velvetoptimiser -f $file_spec -t 1 --o '-min_contig_lgth " . set_option("contig_min") . "'";
      &shortlog("Auto optimization of velvet kmer: $command");
      $output = `$command 2>&1`;
      &log($output);
      opendir DIR, $velvetoptimizer_output_dir;
      while ($velvetlog = readdir(DIR)) {
        last if $velvetlog =~ /Logfile.txt$/;
      }
      closedir DIR;
      if (-f "$velvetoptimizer_output_dir/$velvetlog") {
        open VLOG, "$velvetoptimizer_output_dir/$velvetlog";
        while ($vlog = <VLOG>) {
          # the final options should be the last one listed
          if ($vlog =~ /^Velvetg parameter string: (.*)$/) {
            $velvetg_options = $1;
          }
        }
        close VLOG;
        if ($velvetg_options =~ /auto_data_(\d+)/) {
          $velveth_k = $1;
        }
        # stop if we hit our limits or are done
        if ($velveth_k == 31 || $velveth_k == 255 ||
            ($kstart < $velveth_k && $velveth_k < $kend) ) {
          $best_middle = 1;
          $kstart = $velveth_k;
          $kend = $velveth_k;
          &shortlog("Best kmer: $velveth_k");
        } else {
          if ($velveth_k == $kstart) {
            $kstart -= 18;	# don't move the full 20 so we overlap by 2,
				# in case the optimal is the edge of one of our
				# intervals
            $kstart = 31 if $kstart < 31;
          } elsif ($velveth_k == $kend) {
            $kstart += 18;
            $kstart = 255 if $kstart > 255;	# though this shouldn't happen
          }
        }
      } else {
        # we failed at this velvetoptimiser run for some reason...
        $best_middle = 1;
        $kend = set_option("VELVETOPTIMIZER", "e");
        $kstart = $kend - 20;
        &log("Couldn't auto-optimize velvet assembly, reverting to range $kstart to $kend for kmer");
      }
      # clean things up - space is one of the reasons to do this optimization
      # VelvetOptimiser should already clean up the auto_data dirs
      if (-d $velvetoptimizer_output_dir) {
        File::Path::remove_tree($velvetoptimizer_output_dir);
      }
    }
  }

  # the final VelvetOptimiser run
  $command = "$programs->{VelvetOptimiser} -s $kstart -e $kend -d $tempdir/velvetoptimiser -f $file_spec -t 1 --o '-min_contig_lgth " . set_option("contig_min") . "'";
  $velvet_output_contigs = "$velvetoptimizer_output_dir/contigs.fa";
  $expected_file = "$tempdir/$final_sample_name.contigs.txt";
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  # VelvetOptimizer runs velvetg with the -clean yes switch.  This deletes the
  # LastGraph file that's needed for FinIS finishing...
  # the last velvetg command is in the log file in that directory
  # velvetg on the final directory
  $velvetg_options = "-clean no";
  $velveth_k = 0;
  my $velvetg_exp_cov = 0;
  my $velvetg_cov_cutoff = 0;
  opendir DIR, $velvetoptimizer_output_dir;
  while ($velvetlog = readdir(DIR)) {
    last if $velvetlog =~ /Logfile.txt$/;
  }
  closedir DIR;
  if (-f "$velvetoptimizer_output_dir/$velvetlog") {
    open VLOG, "$velvetoptimizer_output_dir/$velvetlog";
    while ($vlog = <VLOG>) {
      # the final options should be the last one listed
      if ($vlog =~ /^Velvetg parameter string: (.*)$/) {
        $velvetg_options = $1;
      }
    }
    close VLOG;
    if ($velvetg_options =~ /auto_data_(\d+)/) {
      $velveth_k = $1;
    }
    if ($velvetg_options =~ /-exp_cov\s+(\d+)/) {
      $velvetg_exp_cov = $1;
    }
    if ($velvetg_options =~ /-cov_cutoff\s+(\d+\.?\d*)/) {
      $velvetg_cov_cutoff = $1;
    }
    $velvetg_options =~ s/-clean yes/-clean no/;
    $velvetg_options =~ s/auto_data_\d+//;
  }

  $command = "$programs->{velvetg} $velvetoptimizer_output_dir $velvetg_options";
  &shortlog("velveth kmer: $velveth_k");
  &shortlog("velvetg exp_cov: $velvetg_exp_cov");
  &shortlog("velvetg cov_cutoff: $velvetg_cov_cutoff");
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  
  # this gets put into $tempdir/velvetoptimiser/
  if (-f $velvet_output_contigs) {
    File::Copy::copy("$velvet_output_contigs", $expected_file);
    if (!$single_fastq && defined $shuffle && -f $shuffle) {
      unlink($shuffle);
    }
    return ($expected_file);
  } else {
    &log("Couldn't find $velvet_output_contigs after ".(caller(0))[3].", velvet optimiser run");
    if (!$single_fastq && defined $shuffle && -f $shuffle) {
      unlink($shuffle);
    }
    return (0);
  }
}

sub scaffold_GAPCLOSER {
  my ($q1, $q2, $contigs) = @_;
  my $tempdir = set_option("tempdir");
  my $expected_file = "$tempdir/$final_sample_name.gapFilled";
  my $config_file = config_SOAP($q1, $q2);
  my $command = "$programs->{GapCloser} -b $config_file -a $contigs";
  $command .= " -p " . set_option("SOAP", "overlap");
  $command .= " -t " . set_option("SOAP", "n_cpu");
  $command .= " -o $expected_file";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  if (-f $expected_file) {
    return ($expected_file);
  } else {
    &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
    return (0);
  }
}

sub scaffold_OPERA {
  my ($q1, $q2, $contigs) = @_;
  my $tempdir = set_option("tempdir");
  my $opera_output_dir = "$tempdir/OPERA";
  my $mapper = set_option("OPERA", "mapper");
  my $expected_file = "$tempdir/$final_sample_name.gapFilled";
  my $map_file = "$tempdir/$final_sample_name.read-on-contig";
  if ($q1 eq $q2) {
    File::Copy::copy($contigs, $expected_file);
    &shortlog("q1 and q2 files the same.  Skipping scaffolding.");
  } else {
    # syntax changed at least as of OPERA_LG_v2.0.6
    my $command = "$programs->{OPERA_preprocess} --contig $contigs --illumina-read1 $q1 --illumina-read2 $q2 --out $map_file --map-tool $mapper";
    &shortlog($command);
    my $output = `$command 2>&1`;
    # OPERA's preprocess leaves this read.sai file
    unlink("$tempdir/read.sai");
    &log($output);
    if (!-f $map_file) {
      &log("Couldn't find mapping file $map_file after ".(caller(0))[3]." run");
      return(0);
    }
    $command = "$programs->{OPERA} $contigs $map_file $opera_output_dir";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (-f "$opera_output_dir/scaffoldSeq.fasta") {
      File::Copy::copy("$opera_output_dir/scaffoldSeq.fasta", $expected_file);
    } else {
      &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
      return(0);
    }
  }

  return($expected_file);
}

sub finish_FinIS {
  my ($velvet_folder, $scaffolds, $output_dir) = @_;
  my $expected_file;
  # as of 0.3 the number of threads is a required command line parameter
  my $command = "$programs->{FinIS} $velvet_folder $scaffolds $output_dir " . set_option("threads");
  my $output;
  my $finis_file = "$output_dir/scaffolds.filled.fasta";
  # check for what we need
  if ((-f "$velvet_folder/contigs.fa" && -f "$velvet_folder/LastGraph") ||
      (-f "$velvet_folder/$final_sample_name.scafSeq")
     ) {
    mkdir($output_dir);
    $expected_file = "$tempdir/$final_sample_name.gapFilled.FinIS";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (-f $finis_file) {
      File::Copy::copy($finis_file, $expected_file);
    } else {
      &log("Couldn't find $finis_file after ".(caller(0))[3]." run");
      return(0);
    }
    return($expected_file);
  } else {
    &shortlog("Couldn't find contigs.fa or LastGraph file in $velvet_folder");
    return 0;
  }
}

sub species_kraken {
  my ($q1, $q2) = @_;
  my $expected_file = "kraken-out.tmp";
  my $command = "$programs->{kraken}";
  my $kraken_ds1_filename;
  my $kraken_ds2_filename;
  my $kraken_dscommand;
  my @f;
  my @g;
  my $i;
  my $max;
  my $backup_max;
  my $output;
  my $return;
  my $return_index;
  my $backup_return;
  my $backup_return_index;
  my $n = set_option("KRAKEN", "reads");

  # just return if this was set on the command line
  $return = set_option("KRAKEN", "classification");
  if (defined $return && uc($return) ne "NONE" && $return ne "") {
    return $return;
  }

  $command .= " --threads " . set_option("threads");
  $command .= " --db " . set_option("KRAKEN", "DB");
  $command .= " --fastq-input";
  $kraken_ds1_filename = $q1 . "-" . set_option("MAXREADS") . "_" . set_option("SEED") . ".fastq";
  $kraken_ds2_filename = $q2 . "-" . set_option("MAXREADS") . "_" . set_option("SEED") . ".fastq";
  $kraken_dscommand = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q1 " . $n . " > $kraken_ds1_filename";
  &shortlog($kraken_dscommand);
  $output = `$kraken_dscommand 2>&1`;
  if ($q1 eq $q2) {
    $command .= " $kraken_ds1_filename 2>&1 > $expected_file";
  } else {
    $kraken_dscommand = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q2 " . $n . " > $kraken_ds2_filename";
    &shortlog($kraken_dscommand);
    $output = `$kraken_dscommand 2>&1`;
    $command .= " --paired $kraken_ds1_filename $kraken_ds2_filename 2>&1 > $expected_file";
  }
  &shortlog($command);
  # for the kraken run, we're using the trick of 'command 2>&1 > output'
  # this puts stdout into the output file, and stderr shows up on stdout
  # so we can capture stderr with backticks
  # this is because we want the classifications into the output file
  # and the summary comes on stderr which we want to parse
  $output = `$command`;
  &log($output);
  $output =~ s/\r/\n/g;
  @f = split /\n/, $output;
  foreach $i (@f) {
    if ($i =~ /(\d+) sequences .* processed/) {
      &shortlog("Kraken: $1 total sequences processed");
    } elsif ($i =~ /^\s*\d+ sequences classified/) {
      $i =~ s/^\s+//;
      &shortlog("Kraken: $i");
    } elsif ($i =~ /^\s*\d+ sequences unclassified/) {
      $i =~ s/^\s+//;
      &shortlog("Kraken: $i");
    }
  }
  if (-f $expected_file) {
    $command = "$programs->{krakenreport} --db " . set_option("KRAKEN", "DB") . " $expected_file";
    &shortlog($command);
    $output = `$command`;
    &log($output);
    @f = split /\n/, $output;
    $max = 0;
    $backup_max = 0;
    $return = "";
    $return_index = -1;
    $backup_return = "";
    $backup_return_index = -1;
    foreach $i (0..$#f) {
      $f[$i] =~ s/^\s+//;
      @g = split /\t/, $f[$i];
      if ($g[3] eq "G") {
        if ($g[0] > $backup_max) {
          $backup_max = $g[0];
          $backup_return = $g[5];
          $backup_return =~ s/^\s+//;
          $backup_return_index = $i;
        }
      }
      if ($g[3] eq "S") {
        if ($g[0] > $max) {
          $max = $g[0];
          $return = $g[5];
          $return =~ s/^\s+//;
          $return_index = $i;
        }
      }
    }
    &log("Cleaning up kraken intermediate files $expected_file and downsampled inputs...");
    unlink($expected_file);
    unlink($kraken_ds1_filename);
    unlink($kraken_ds2_filename);
    if ($return_index >= 0) {
      &shortlog("Kraken classification line: $f[$return_index]");
      return($return);
    } elsif ($backup_return_index >= 0) {
      &shortlog("Kraken classification line: $f[$backup_return_index]");
      return($backup_return);
    }
  } else {
    &shortlog("Couldn't find $expected_file file after initial Kraken run, skipping...");
    return("");
  }
  return("");
}

sub select_reference {
  my ($assembly, $num_to_test, $preferred_reference, $initial) = @_;
  my $tempdir = set_option("tempdir");
  $tempdir = "$tempdir/reference_genomes";
  mkdir($tempdir);
  if (!defined $num_to_test || $num_to_test <= 0) {
    $num_to_test = 200;
  }
  my $command = "$programs->{closest_species} -num_test $num_to_test $assembly -download $tempdir";
  if ($initial ne "") {
    $command .= " -initial $initial";
  }
  my $output = `$command`;
  my @reference = split /\n/, $output;
  foreach my $i (reverse 0..$#reference) {
    if (!-f $reference[$i]) {
      splice @reference, $i, 1;
    }
  }
  # if there was a reference specified, we'll use that as the "best"
  # move it around if it's already in the list, if not then just stick it
  # at the front
  if (defined $preferred_reference && length $preferred_reference && -f $preferred_reference) {
    unshift @reference, $preferred_reference;
    if (scalar @reference > 1) {
      foreach my $i (1..$#reference) {
        if ($preferred_reference eq $reference[$i]) {
          splice @reference, $i, 1;
          last;
        }
      }
    }
  }
  return (@reference);
}

sub compare_NUCMER {
  my ($reference, $assembly, $prefix) = @_;
  my $tempdir = set_option("tempdir");
  my @return = ();
  my $ordering;
  my $expected_file = "$prefix.delta";
  my $mummer_plot_output = $prefix . ".mummerplot";
  my $command = "$programs->{nucmer} --maxmatch $reference $assembly -p $prefix";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  if (!-f $expected_file) {
    &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
    return (0);
  }
  push @return, $expected_file;

  $command = "$programs->{mummerplot} $expected_file -p $mummer_plot_output -l --postscript";
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  $ordering = ();
  open GP, "$mummer_plot_output.gp";
  while (<GP>) {
    chomp;
    if (/^set ytics/) {
      while (<GP>) {
        chomp;
        if (/^\s*"(.*?)\"\s\d+/) {
          push @{$ordering}, $1;
        }
        last if /^\)/;
      }
    }
  }
  close GP;
  # clean up mummerplot stuff
  unlink("$mummer_plot_output.filter");
  unlink("$mummer_plot_output.fplot");
  unlink("$mummer_plot_output.rplot");
  unlink("$mummer_plot_output.gp");
#
# we don't need to die here if we can't find the mummerplot file? - no expected file was set
#
#  if (!-f $expected_file) {
#    &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
#    return (0);
#  }

  # process a little bit - coords
  $expected_file = "$prefix.coords";
  $command = "$programs->{showcoords} -HTclq $prefix.delta > $expected_file";
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
   if (!-f $expected_file) {
     &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
    return (0);
  }
  push @return, $expected_file;

  # process a little bit - snps
  $expected_file = "$prefix.snps";
  $command = "$programs->{showsnps} -CT $prefix.delta > $expected_file";
  &shortlog($command);
  $output = `$command 2>&1`;
  &log($output);
  if (!-f $expected_file) {
    &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
    return (0);
  }

  push @return, $expected_file;
  
  return (@return, $ordering);
}

sub annotate_prokka {
  my ($assembly, $prefix, $refspecies) = @_;
  my $genus;
  my $species;
  my @wanted = qw(faa ffn fna gbk gff);
  my $outdir = "$tempdir/prokka";
  my $species_flag;
  my $command;
  my $output;
  my $i;
  my $ext;
  my $error;
  undef $genus;
  undef $species;
  ($genus, $species) = split / /, $refspecies;
  if (defined $genus && length $genus > 0) {
    if (defined $species && length $species > 0) {
      $species_flag = "--species $species";
    }
    $command = "$programs->{prokka} --outdir $outdir --prefix $prefix --genus $genus $species_flag --quiet --cpus " . set_option("threads") . " --strain $prefix --kingdom Bacteria --gcode 0 --addgenes --locustag $prefix --mincontiglen " . set_option("contig_min") . " --usegenus $assembly";
  } else {
    # if no genus, use neither genus nor species
    $command = "$programs->{prokka} --outdir $outdir --prefix $prefix --quiet --cpus " . set_option("threads") . " --strain $prefix --kingdom Bacteria --gcode 0 --addgenes --locustag $prefix --mincontiglen " . set_option("contig_min") . " $assembly";
  }
  &shortlog($command);

  # overall output log
  $output = `$command 2>&1`;
  # should be no output here actually - get log from the log file
  if (-f "$outdir/$prefix.log") {
    $output = "";
    open PROKKAF, "$outdir/$prefix.log";
    while ($i = <PROKKAF>) {
      $output .= $i;
    }
  }
  &log($output);

  # summary of the run
  if (-f "$outdir/$prefix.txt") {
    $output = "\n==== Prokka Summary ====\n";
    open PROKKAF, "$outdir/$prefix.txt";
    while ($i = <PROKKAF>) {
      $output .= $i;
    }
  }
  &shortlog($output);

  $error = 0;
  foreach $ext (@wanted) {
    if (-f "$outdir/$prefix.$ext") {
      File::Copy::copy("$outdir/$prefix.$ext", "$tempdir/$prefix.$ext");
    } else {
      &log("Couldn't find $prefix.$ext after ".(caller(0))[3]." run");
      $error++;
    }
  }
  if ($cleanup) {
    File::Path::remove_tree($outdir);
  }
  return($error);
}

sub gene_prodigal {
  my ($assembly, $prefix) = @_;
  my $prodigal_temp;
  my $seq;
  my $out;
  my $sequence;
  my $annotation;
  my @data;
  my @f;
  my $name;
  my $locus;
  my $gbk;
  my $hdr;
  my $s;
  my $found_origin;
  my $expected_file = "$prefix.prodigal";
  my $command = "$programs->{prodigal} -i $assembly -m -o $expected_file -a $prefix.faa -d $prefix.ffn";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);

  if (-f $expected_file) {
    # combine the sequence and annotation into a .gbk file that can actually
    # be imported into programs like VectorNTI
    $gbk = $expected_file;
    $expected_file = "$prefix.gbk";
    $prodigal_temp = "$tempdir/prodigal.tempfile";
    $seq = Bio::SeqIO->new(-file=>$assembly);
    $out = Bio::SeqIO->new(-file=>">$prodigal_temp",-format=>"genbank");
    while ($s = $seq->next_seq) {
      $out->write_seq($s)
    }
    $sequence = ();
    $annotation = ();

    open PRODIGAL, $prodigal_temp;
    @data = ();
    $found_origin = 0;
    while (<PRODIGAL>) {
      if (/^LOCUS/) {
        @f = split /\s+/, $_;
        $name = $f[1];
        @data = ();
        $found_origin = 0;
        $locus->{$name} = $_;
        next;
      }
      $found_origin = 1 if /^ORIGIN/;
      next if !$found_origin;
      push @data, $_;
      if (/^\/\//) {
        $sequence->{$name} = join ("", @data);
        $name = "";
        next;
      }
    }
    close PRODIGAL;
    unlink $prodigal_temp;

    open PRODIGAL, $gbk;
    while (<PRODIGAL>) {
      if (/^DEFINITION/) {
        /seqhdr="(.*?)"/;
        $hdr = $1;
        @f = split /\s+/, $hdr;
        $name = $f[0];
        @data = ();
      }
      if (/^\/\//) {
        $annotation->{$name} = join ("", @data);
        $name = "";
        next;
      }
      push @data, $_;
    }
    close PRODIGAL;
    
    open PRODIGAL, ">$expected_file";
    foreach $name (keys %$locus) {
      print PRODIGAL $locus->{$name};
      if (defined $annotation->{$name}) {
        print PRODIGAL $annotation->{$name};
      } else {
        print PRODIGAL "FEATURES             Location/Qualifiers\n";
      }
      print PRODIGAL $sequence->{$name};
    }
    close PRODIGAL;

    return ($expected_file);
  } else {
    &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
    return (0);
  }
}

sub novel_NUCMER {
  my $assembly = pop;
  my @reference = @_;
  my $contig_min = set_option("contig_min");
  my $tempdir = set_option("tempdir");
  my ($delta, $coords, $snps);
  my $output_base;
  my $command;
  my $all_refs = "$tempdir/all_references.fna";
  my $in;

  open ALL, ">$all_refs";
  foreach my $reference (@reference) {
    open REF, $reference;
    while ($in = <REF>) { print ALL $in; }
    close REF;
  }
  close ALL;
  $output_base = "$tempdir/$final_sample_name-species_comparison";
  ($delta, $coords, $snps) = compare_NUCMER($all_refs, $assembly, $output_base);
  $command = "$programs->{selectnovel} $coords $assembly $output_base";
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
  if ($cleanup) {
    unlink ("$output_base.delta");
    unlink ("$output_base.coords");
    unlink ("$output_base.snps");
    unlink ("$all_refs");
  }
}

sub snp_info {
  my ($snpfile, $fnafile) = @_;
  my $pttfile;
  my $ptt_ref;
  my $sequence;
  my $data;
  my $snpline;
  my $position;
  my $mutation;
  my $line;
  my @f;
  my @out;
  my @info_output;
  my $output_file = $snpfile;
#  $output_file .= ".annotated";

  if ($fnafile =~ /fna$/) {
    # this should be a Genbank genome, so we should have a ptt file
    $pttfile = $fnafile;
    $pttfile =~ s/fna$/ptt/;
  }
  ($ptt_ref, $sequence) = GERMS::ns_setup($pttfile, $fnafile);
  @info_output = ();
  if (defined $ptt_ref && defined $sequence) {
    open SNP, $snpfile;
    while ($line = <SNP>) {
      chomp $line;
      @f = split /\t/, $line;
      if ($#f < 7 || $f[0] !~ /^\d+$/) {
        push @info_output, "$line\n";
        next;
      }
      @out = @f;
      $position = $f[0];
      $mutation = $f[2];
      # this doesn't handle consecutive mutations right now
      # It also needs a nearestorf component for intergenics
      $data = GERMS::ns($position, $mutation, $ptt_ref, $sequence);
      if ($data->{INTERGENIC}) {
        if ($f[1] eq '.') {
          push @out, "Insertion";
        } elsif ($mutation eq '.') {
          push @out, "Deletion";
        } else {
          push @out, "Substitution";
        }
        push @out, "Intergenic", "", "";
      } else {
        if ($f[1] eq '.') {
          push @out, "Insertion";
        } elsif ($mutation eq '.') {
          push @out, "Deletion";
        } else {
          push @out, "Substitution";
          if ($data->{SYNONYMOUS}) {
            push @out, "Synonymous";
            push @out, "";
          } else {
            push @out, "Nonsynonymous";
            push @out, $data->{ORIGINALAA} . $data->{AAPOSITION} . $data->{NEWAA};
          }
        }
        push @out, $data->{GID}, $data->{SYSTEMATIC}, $ptt_ref->{$data->{GID}}->{TRIVIAL}, $ptt_ref->{$data->{GID}}->{ANNOTATION};
      }
      push @info_output, (join ("\t", @out) . "\n");
    }
    close SNP;
    if (scalar @info_output) {
      unlink($output_file);
      open SNP, ">$output_file";
      print SNP @info_output;
      close SNP;
    }
    return ($output_file);
  } else {
    return "";
  }
}

sub join_assembly {
  my ($assembly, $ordering, $prefix) = @_;
  my @pf;
  my $pfseqs;
  my @pfout;
  my %used;
  my $key;
  my $expected_file = "$assembly-joined";
  &log("Joining assembly for gene calling\n");
  open PF, $assembly;
  open PS, ">$expected_file";
  print PS ">", substr($prefix,0,20), "\n";
  @pf = <PF>;
  # get rid of spaces for other info like length, just keep contig name
  foreach my $i (0..$#pf) {
    if ($pf[$i] =~ /^>/) {
      my @f = split /\s+/, $pf[$i];
      $pf[$i] = $f[0];
    }
  }
  $pfseqs = GERMS::fasta2hash(@pf);
  @pfout = ();
  %used = ();
  undef $key;
  foreach $key (@{$ordering}) {
    if ($key =~ s/^\*//) {
      next if !defined($pfseqs->{$key});
      push @pfout, GERMS::revcomp($pfseqs->{$key});
    } else {
      next if !defined($pfseqs->{$key});
      push @pfout, $pfseqs->{$key};
    }
    $used{$key} = 1;
  }
  foreach $key (keys %$pfseqs) {
   push @pfout, $pfseqs->{$key} if !defined($used{$key}) || !$used{$key};
  }
  print PS join ("N" x 100, @pfout);
  close PS;
  close PF;
  return ($expected_file) if -f $expected_file;
}

sub qc_FASTQC {
  my (@files) = @_;
  my $tempdir = set_option("tempdir");
  my $f;
  my $command = "$programs->{fastqc} -o $tempdir --noextract";
  $command .= " -t " . set_option("threads");
  $command .= " " . join (" ", @files);
  &shortlog($command);
  my $output = `$command 2>&1`;
  &log($output);
}

sub shortlog {
  my ($message) = @_;
  # this will also do a regular log of everything as well
  chomp $message;
  &log($message);
  print SHORTLOG "$message\n";
}

sub log {
  my ($message) = @_;
  print LOG "# ", basename($0), " ## ", scalar(localtime), "\n";
  chomp $message;
  print LOG "# ", basename($0), " ## $message\n";
}

sub do_single_library {
  my ($final_sample_name, $q1, $q2) = @_;

  my $tempparent = set_option("tempparent");
  my $tempdir = set_option("tempdir");
  my $contigs;
  my $assembly;
  my @reference;
  my ($delta, $coords, $snps);
  my $novel;
  my $genes;
  my $contig_stats;
  my $scaffold_stats;
  my $current_dir;
  my $final_information = "\n==== Final Summary ====\n";
  my $tempinfo;
  my $single_fastq = 0;

  if (defined $q2) {
    print STDERR "STARTED $final_sample_name using files $q1 and $q2 in $tempdir by $ENV{USER} on ", `hostname`;
    &shortlog("STARTED $final_sample_name using files $q1 and $q2 in $tempdir by $ENV{USER} on " . `hostname`);
    &shortlog("Files:\n" . `md5sum $0 $q1 $q2`);
  } else {
    print STDERR "STARTED $final_sample_name using file $q1 in $tempdir by $ENV{USER} on ", `hostname`;
    $q2 = $q1;	# this is how we'll deal with single end internally
    $single_fastq = 1;
    &shortlog("STARTED $final_sample_name using file $q1 in $tempdir by $ENV{USER} on ", `hostname`);
    &shortlog("Files:\n" . `md5sum $0 $q1`);
  }
  print STDERR "You can follow the progress by running 'tail -f $tempdir/$final_sample_name-log.txt'; you may need to be logged in to ", `hostname`;
  &shortlog("You can follow the progress using \'tail -f $tempdir/$final_sample_name-log.txt\; you may need to be logged in to " . `hostname`);

  if (!defined $q1 || !-f $q1 || !defined $q2 || !-f $q2) {
    &log("Can't find file $q1 and/or $q2 for $final_sample_name, skipping...\n");
    close LOG;
    close SHORTLOG;
    return;
  }

  # partly this gives us absolute paths for $q1 and $q2 so we can move around
  ($q1, $q2) = GERMS::check_files($q1, $q2);
  $current_dir = getcwd;
  chdir($tempdir);
  my $cleanup_fastq1 = 0;
  my $cleanup_fastq2 = 0;
  if ($q1 =~ /\.gz$/) {
    my $origq = $q1;
    $q1 =~ s/\.gz$//;
    $q1 = $tempdir . "/" . basename($q1);
    my $unzip = `zcat $origq > $q1`;
    &shortlog("Unzipping gzipped $origq into temp directory ($tempdir)");
    $cleanup_fastq1 = 1;
  }
  if ($single_fastq) {
    $q2 = $q1;
  } elsif ($q2 =~ /\.gz$/) {
    my $origq = $q2;
    $q2 =~ s/\.gz$//;
    $q2 = $tempdir . "/" . basename($q2);
    my $unzip = `zcat $origq > $q2`;
    &shortlog("Unzipping gzipped $origq into temp directory ($tempdir)");
    $cleanup_fastq2 = 1;
  }

  # collect info
  $tempinfo = `wc -l $q1 | awk '{printf "\%.0f", \$1/4}'`;
  chomp $tempinfo;
  if ($single_fastq) {
    $final_information .= "Number of single-end reads: $tempinfo\n";
  } else {
    $final_information .= "Number of paired-end reads: $tempinfo\n";
  }
  $kraken_classification = species_kraken($q1, $q2);
  $final_information .= "Kraken classification: $kraken_classification\n";
  if (set_option("MAXREADS") > 0 && set_option("MAXREADS") < $tempinfo) {
    &shortlog("Downsampling to " . set_option("MAXREADS") . " reads with seed " . set_option("SEED") . "\n");
    $final_information .= "Reads used for assembly: " . set_option("MAXREADS") . "\n";
    my $ds_filename = $q1 . "-" . set_option("MAXREADS") . "_" . set_option("SEED");
    my $command = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q1 " . set_option("MAXREADS") . " > $ds_filename";
    &shortlog($command);
    my $output = `$command`;
    if ($cleanup_fastq1) {
      unlink("$tempdir/".basename($q1));
    }
    $q1 = $ds_filename;
    if ($single_fastq) {
      $q2 = $ds_filename;
    } else {
      $ds_filename = $q2 . "-" . set_option("MAXREADS") . "_" . set_option("SEED");
      $command = $programs->{seqtk} . " sample -s " . set_option("SEED") . " $q2 " . set_option("MAXREADS") . " > $ds_filename";
      &shortlog($command);
      $output = `$command`;
      if ($cleanup_fastq2) {
        unlink("$tempdir/".basename($q2));
      }
      $q2 = $ds_filename;
    }
  }
  # this should be read length +1 because of newline
  $tempinfo = `head -n 4 $q1 | grep -A1 '^\@' | tail -n 1 | wc | awk '{printf "\%.0f", \$3 - 1}'`;
  chomp $tempinfo;
  $final_information .= "Read length (based on first read): $tempinfo\n";

  #
  # FastQC
  #
  if (set_option("FASTQC", "dofastqc")) {
    if ($single_fastq) {
      &qc_FASTQC($q1);
    } else {
      &qc_FASTQC($q1, $q2);
    }
  }

  #
  # assembly
  #
  if (uc(set_option("ASSEMBLER")) eq "VELVETOPTIMIZER" ||
      uc(set_option("ASSEMBLER")) eq "VELVETOPTIMISER") {
     $contigs = assemble_velvet($q1, $q2);
    if (!-f $contigs) {
      $contigs = assemble_SOAPDENOVO($q1, $q2);
      $final_information .= "Assembly: SOAPDENOVO\n";
    } else {
      $final_information .= "Assembly: VelvetOptimizer\n";
    }
  } elsif (uc(set_option("ASSEMBLER")) eq "SPADES") {
    $contigs = assemble_SPAdes($q1, $q2);
    if (!-f $contigs) {
      $contigs = assemble_SOAPDENOVO($q1, $q2);
      $final_information .= "Assembly: SOAPDENOVO\n";
    } else {
      $final_information .= "Assembly: SPAdes\n";
    }
  } else {
    $contigs = assemble_SOAPDENOVO($q1, $q2);
    if (!-f $contigs) {
      $contigs = assemble_velvet($q1, $q2);
      $final_information .= "Assembly: VelvetOptimizer\n";
    } else {
      $final_information .= "Assembly: SOAPDENOVO\n";
    }
  }
  $contig_stats = GERMS::assembly_stats($contigs, set_option("contig_min"));
  &shortlog("\n----- Contig stats for initial assembly -----\n$contig_stats->{text}\n-----");
  if ($contig_stats->{total_length} == 0) {
    &log("Couldn't do initial assembly for $final_sample_name, aborting...");
    return;
  }

  #
  # scaffolding and finishing
  #
  if ($q1 ne $q2) {
    &log("Running scaffolding\n");
    if (uc(set_option("SCAFFOLDER")) eq "OPERA") {
      $assembly = scaffold_OPERA($q1, $q2, $contigs);
      $final_information .= "Scaffolding: OPERA\n";
    } else {
      $assembly = scaffold_GAPCLOSER($q1, $q2, $contigs);
      $final_information .= "Scaffolding: GAPCLOSER\n";
    }
    #
    # finishing if we used velvet and OPERA
    #
    if (!-f $assembly) {
      $assembly = $contigs;
      &shortlog("Some problem with scaffolder; continuing with contigs.");
    } elsif (uc(set_option("SCAFFOLDER")) eq "OPERA" && uc(set_option("ASSEMBLER")) ne "SPADES") {
      &log("Running Finishing\n");
      my $finished = "";
      if (uc(set_option("ASSEMBLER")) eq "VELVETOPTIMIZER" || uc(set_option("ASSEMBLER")) eq "VELVETOPTIMISER") {
        $finished = finish_FinIS("$tempdir/velvetoptimiser", "$tempdir/OPERA/scaffolds.scaf", "$tempdir/FinIS");
      } elsif (uc(set_option("ASSEMBLER")) eq "SOAP") {
        $finished = finish_FinIS("$tempdir/SOAP", "$tempdir/OPERA/scaffolds.scaf", "$tempdir/FinIS");
      }
      if (-f $finished) {
        $assembly = $finished;
        $final_information .= "Finishing: FinIS\n";
      } else {
        $final_information .= "Finishing: Failed for some reason\n";
      }
    } else {
      $final_information .= "Finishing: None\n";
    }
  } else {
    &log("Single end.  No scaffolding or finishing.\n");
    $assembly = $contigs;
  }

  $scaffold_stats = GERMS::assembly_stats($assembly, set_option("contig_min"));
  &shortlog("\n---- Final assembly stats -----\n$scaffold_stats->{text}\n-----");
  if ($scaffold_stats->{total_length} == 0) {
    &log("Got zero total length during scaffolding for $final_sample_name, aborting...");
    return;
  }

  &shortlog("Renaming final assembly to $final_sample_name.assembly\n");
  rename($assembly, "$tempdir/$final_sample_name.assembly");
  $assembly = "$tempdir/$final_sample_name.assembly";
  if ($cleanup) {
    &shortlog("Cleaning up intermediate assembly files - contigs.txt, gapFilled, etc.\n");
    unlink("$tempdir/$final_sample_name.contigs.txt");
    unlink("$tempdir/$final_sample_name.gapFilled");
    unlink("$tempdir/$final_sample_name.spades.fasta");
  }
  if ($cleanup_fastq1) {
    unlink("$tempdir/".basename($q1));
  }
  if ($cleanup_fastq2) {
    unlink("$tempdir/".basename($q2));
  }
  $final_information .= $scaffold_stats->{text};

  #
  # choose reference sequence, SNP calling, novel sequences
  #
  if (defined $scaffold_stats->{N50_number} && $scaffold_stats->{N50_number} > 0) {
    @reference = select_reference($assembly, $scaffold_stats->{N50_number}, set_option("REFERENCE"), $kraken_classification);
  } elsif (defined $contig_stats->{N50_number} && $contig_stats->{N50_number} > 0) {
    @reference = select_reference($assembly, $contig_stats->{N50_number}, set_option("REFERENCE"), $kraken_classification);
  } else {
    # 200 is the closest_species default for # of sequences to look at
    @reference = select_reference($assembly, 200, set_option("REFERENCE"), $kraken_classification);
  }
  if (scalar(@reference)) {
    &shortlog(join ("\n", "References:", @reference));
    # $ordering is a reference to an array for mummer based contig ordering
    &log("Calling NUCMER with $reference[0], $assembly, $tempdir/$final_sample_name-".basename($reference[0])."\n");
    ($delta, $coords, $snps, $ordering) = compare_NUCMER($reference[0], $assembly, "$tempdir/$final_sample_name-".basename($reference[0]));
    &log("Getting snp info with $snps, $reference[0]\n");
    &snp_info($snps, $reference[0]);
    &log("Getting novel sequences with NUCMER\n");
    $novel = novel_NUCMER(@reference, $assembly);
  } else {
    &shortlog("Can not find a reference genome, skipping running SNP detection  ...\n");
  }
  $final_information .= "Reference sequence: $reference[0]\n";

  #
  # Annotation
  #
  # old way - run prodigal and do a blast-based annotation to best reference
  # this unfortunately doesn't have all the files we might want, though they
  # can be generated - ex. .gff
  # new way - prokka - this by default gives all the ncbi refseq files
  # including .fna, .faa, .ffn, .gff, .gbk
  # pick genus / species from closest_species
  #
  my $joined_assembly = join_assembly($assembly, $ordering, $final_sample_name);
  if (uc(set_option("ANNOTATOR")) eq "PROKKA") {
    &annotate_prokka($joined_assembly, $final_sample_name, $kraken_classification);
  } else {
    #
    # gene prediction
    #
    # at first we call on just the final assembly - so can trace back to contigs
    # don't do this until we join everything together - this is more for annotation
    &gene_prodigal($joined_assembly, "$tempdir/$final_sample_name");

    #
    # do a blast-based annotation against the best reference
    #
    if (scalar(@reference)) {
      &shortlog("Annotating $final_sample_name.gbk / $final_sample_name.faa using $reference[0]\n");
      &GERMS::annotate_by_ref("$final_sample_name.faa", "$final_sample_name.gbk", $reference[0]);
    }
  }
  unlink($joined_assembly);

  chdir($current_dir);
  return($final_information);
}

sub post_process {
  my ($final_sample_name) = @_;
  my $tempdir = set_option("tempdir");
  my $tempparent = set_option("tempparent");
  if (!$cleanup) {
    print STDERR "NOTE -- Intermediate files have been left in:\n$tempdir\n";
    &shortlog("NOTE -- Intermediate files have been left in:\n$tempdir\n");
  } else {
    # fastq files
    unlink("$tempdir/$final_sample_name"."_R1.fastq");
    unlink("$tempdir/$final_sample_name"."_R2.fastq");
    foreach my $ext (@not_wanted) {
      unlink("$tempdir/$final_sample_name.$ext");
      unlink("$tempdir/$final_sample_name-species_comparison.$ext");
    }
    for my $i (set_option("VELVETOPTIMIZER", "s") .. set_option("VELVETOPTIMIZER", "e")) {
      if (-d "$tempdir/auto_data_$i") {
        # these are velvetoptimizer directories, they get left if for some
        # reason velvetoptimizer could do a final assembly
        File::Path::remove_tree("$tempdir/auto_data_$i");
      }
    }
    if (-d "$tempdir/OPERA") {
      File::Path::remove_tree("$tempdir/OPERA");
    }
    if (-d "$tempdir/velvetoptimiser") {
      File::Path::remove_tree("$tempdir/velvetoptimiser");
    }
    if (-d "$tempdir/SOAP") {
      File::Path::remove_tree("$tempdir/SOAP");
    }
    if (-d "$tempdir/FinIS") {
      File::Path::remove_tree("$tempdir/FinIS");
    }
    if (-d "$tempdir/spades_working") {
      File::Path::remove_tree("$tempdir/spades_working");
    }
    if (-d "$tempdir/reference_genomes") {
      # this would have been made by select_reference
      File::Path::remove_tree("$tempdir/reference_genomes");
    }
  }
  &log("tar cvzf $output_dir/$final_sample_name.tgz $final_sample_name --exclude=*_fastq.txt");
  close LOG;
  close SHORTLOG;
  chdir($tempparent);
  `tar cvzf $output_dir/$final_sample_name.tgz $final_sample_name --exclude=*_fastq.txt`;
  if (-f "$output_dir/$final_sample_name.tgz" && -s "$output_dir/$final_sample_name.tgz") {
    print STDERR "FINISHED Final output is in $output_dir/$final_sample_name.tgz\n";
    # we don't shortlog here because the file is closed
  } else {
    if (-f "$output_dir/$final_sample_name.tgz") {
      print STDERR "Error - somehow didn't create final output file $output_dir/$final_sample_name.tgz\n";
      # we don't shortlog here because the file is closed
    } else {
      print STDERR "Error - final output file $output_dir/$final_sample_name.tgz has length 0!\n";
      # we don't shortlog here because the file is closed
    }
  }
  return ("$output_dir/$final_sample_name.tgz");
}

sub set_option {
  my (@key) = @_;
  my $return = "";
  my $done = 0;
  my $ref = $user_options;
  foreach my $i (0..$#key) {
    if (defined $ref->{$key[$i]}) {
      $ref = $ref->{$key[$i]};
      if ($i == $#key) {
        $done = 1;
        $return = $ref;
      }
    }
  }
  if (!$done) {
    $ref = $default_options;
    foreach my $i (0..$#key) {
      if (defined $ref->{$key[$i]}) {
        $ref = $ref->{$key[$i]};
        if ($i == $#key) {
          $done = 1;
          $return = $ref;
        }
      }
    }
  }
  return $return;
}

sub print_usage {
  my $name = basename($0);
  print <<__USAGE__;
Usage example:
---------------
   GERMS-wgs.pl -q1 <fastq file> -q2 <fastq file> -name <output file name> -output_dir <output directory> -tempdir <temp directory>
   OR
   GERMS-wgs.pl -mux <MUXID> [ -libname <LIBNAME> ] [ -cluster|nocluster ]
   

This will give you a file named <output file name>.tgz in the directory <output directory>.  If you choose velvetoptimizer and opera, it will run FinIS for you also.

Common options:
---------------
  -output_dir <dir>                  --  where output files go - you probably
                                           always want to specify this
  -tempdir                           --  temporary directory to use -
                                           default $temp_base
  -assembler <velvetoptimizer|spades|soap>
                                     --  which assembler to use - default $default_options->{ASSEMBLER}
  -scaffolder <opera|soap>           --  which scaffolder - default $default_options->{SCAFFOLDER}
  -annotator <prokka|old>            --  how to annotate - default $default_options->{ANNOTATOR}
  -reference <fna file>              --  force a certain reference sequence -
                                           this assumes NCBI refseq files are
                                           available - ptt, faa, ffn
  -t <threads>                       --  number of CPUs to use - default $default_options->{threads}
  -memory <int>                      --  memory to specify (in GB) - default $default_options->{MEMORYINGB}
  -contig_min <int>                  --  contig minimum cutoff - default $default_options->{contig_min}
  -cleanup|nocleanup                 --  whether to leave intermediate files -
                                           default cleanup
  -qc|noqc                           --  whether to run FastQC
  -cluster|nocluster                 --  -cluster doesn't run anything but gives
                                           commands you can use to submit to the
                                           cluster
  -maxreads                          --  max number of reads to use for assembly
                                           default $default_options->{MAXREADS}
                                           if set to 0, then use everything
  -seed                              --  seed used for seqtk to downsample if
                                           we need to use fewer reads
                                           default $default_options->{SEED}
  -species <species>                 --  If set, won't rerun Kraken to get a
                                           species identification
                                           default $default_options->{KRAKEN}->{classification} - meaning run Kraken

Advanced options:
-----------------------------------
  -soap option=value
  -velvetoptimizer option=value

    For each of these, specify the exact option for that program.  This script does no checking on these options/values.  If you need to specify more than one option, use multiple switches. For example, the defaults for SOAP and Velvetoptimizer are currently set as follows:

    -soap k=41 -soap avg_ins=300
    -velvetoptimizer s=51 -velvetoptimizer e=99

    N.B. For velvetoptimizer there are special options for start and end. The default start is 0, which indicates to do auto-selection and optimization of the k-mer based on the read length. This also will try to use less temp drive space by running velvetoptimizer multiple times if needed. If this fails, the fallback is to just take whatever comes out of velvetoptimizer based on the value of the end parameter (default 99; the start value will be 20 less than this). (i.e. default behavior is --velvetoptimizer s=0 -velvetoptimizer e=99)
__USAGE__
}
