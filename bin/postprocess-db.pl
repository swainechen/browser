#!/usr/bin/perl -w
#
# need a full set of files - look for these based on read
# .srst2.gz - goes into Resistance and Genes
# .gcov.gz - some info for Fastq
# .lofreq.gz - not really needed here
# .tgz - info for Fastq and Assembly
#
use warnings;
use strict;
use slchen;
use Archive::Tar;
use DBI;
use DateTime::Format::DBI;
use DateTime::Format::DateParse;
use Data::Dumper;
use File::Spec;
use File::Basename;
use Getopt::Long;
use JSON;
use GERMS;

# Needs configuration
my $S3_BUCKET = "";
my $S3_SPECIES = "";
my $S3_REFDIR = "";
my $S3_PROFILE = "";
my $secondary_S3 = "";
my $secondary_PROFILE = "";

my $S3_OP = "cp";
my $OUTBREAK_BASE = $ENV{"OUTBREAK_BASE"};
$OUTBREAK_BASE =~ s/\/$//;

my $USE_DB = 0;
my $verbose = 0;
my $MLSTDatabase = "";
my $SOURCE = "";
my $force = "";	# update no matter what - ignore date stamps
my @data_dirs = ();
my $push_s3 = 1;
my $delete = 0;
my $reference_fna = "";
&Getopt::Long::Configure("pass_through");
GetOptions (
  'mlst_database=s' => \$MLSTDatabase,
  'source=s' => \$SOURCE,
  'db!' => \$USE_DB,
  'verbose!' => \$verbose,
  'force=s' => \$force,
  'data_dir=s' => \@data_dirs,
  's3!' => \$push_s3,
  'delete!' => \$delete,
  'reference=s' => $reference_fna
);

push @data_dirs, ".";
my $data;
my $parsed_data;
my @f;
my @g;
my $class;
my $data_keys;
my @out;
my $id;
my $i;
my $j;
my $sth;
my $sql;
my @array;
my $param_number;
my $param_map;
my @return;
my $operation;
my $fastq_data;
my $mlst_data;
my $resistance_data;
my $genes_data;
my $status;
my $s3_object;
my $s3_date;
my $s3_MD5;
my $s3_temp;
my $s3_command;
my $secondary;
my $secondary_MD5;
my @order = qw(
  AGly
  Bla
  Colistin
  Fcyn
  Flq
  Gly
  MLS
  Phe
  Rif
  Sul
  Tet
  Tmt
);
my $classes = {
  "AGly" => 1,
  "Bla" => 1,
  "Colistin" => 1,
  "Fcyn" => 1,
  "Flq" => 1,
  "Gly" => 1,
  "MLS" => 1,
  "Phe" => 1,
  "Rif" => 1,
  "Sul" => 1,
  "Tet" => 1,
  "Tmt" => 1
};

sub print_usage {
  print "Usage: $0 <Run ID> -mlst_database <organism> -source <source> [ -data <dir> [ -data <dir> [ ... ] ] ] [ -db ] [ -verbose ] [ -force <table> ] [ -s3|nos3 ] [ -delete|nodelete ] [ -reference <reference fna> ]\n";
  print "  Multiple -data_dir options can be specified and will look through all of those\n";
  print "  Need to specify MLST database because it's not present in the srst2.gz file\n";
  print "  Also need to specify the source - this should be FASTQ or ASSEMBLY\n";
  print "  With -db do the database import / update\n";
  print "  With -force <table> then always update that table\n";
  print "    Choices are: All, Tips, Files, Fastq, MLST, Resistance, Genes, Assembly\n";
  print "  With -verbose print out all the information collected\n";
  print "  With -s3 push files to S3, configuration hardcoded for now in the script\n";
  print "  With -delete then delete files after processing (probably only use with -s3)\n";
  print "  -reference specifies the reference .fna file used if not using the default - this needs to be in the ReferenceGenome table\n";
  if ($MLSTDatabase eq "") {
    print "ERROR: No -mlst_database specified\n";
  }
  if ($SOURCE eq "") {
    print "ERROR: No -source specified\n";
  }
}

if (!defined $ARGV[0] || !length($ARGV[0]) || $MLSTDatabase eq "" || $SOURCE eq "") {
  &print_usage;
  exit;
}
# check for files
my $runID = $ARGV[0];
# don't store lofreq.gz.tbi because it's usually fast to calculate
# but gcov.png takes a while to calculate
my @extensions = qw(lofreq.gz gcov.gz gcov.png srst2.gz tgz);
my $inputs = ();
my $dates = ();
my $dir;
my $ext;
foreach $dir (@data_dirs) {
  foreach $ext (@extensions) {
    if (-f "$dir/$runID.$ext") {
      $inputs->{$ext} = File::Spec->rel2abs("$dir/$runID.$ext");
      $inputs->{$ext . "_MD5"} = `md5sum $inputs->{$ext} | cut -f1 -d' '`;
      chomp($inputs->{$ext . "_MD5"});
      $dates->{$ext} = DateTime->from_epoch(epoch => (stat($inputs->{$ext}))[9]);
    }
  }
}
my $missing_file = 0;
foreach $ext (@extensions) {
  if (!defined $inputs->{$ext} || !-f $inputs->{$ext}) {
    print STDERR "Can't find $runID.$ext\n";
    $missing_file++;
  }
}
if ($missing_file == scalar @extensions) {
  print "Data directories:\n", join ("\n", @data_dirs), "\n";
  &print_usage;
  exit;
}
if ($verbose) {
  print "Found these files:\n";
}
my $files_data;
my $DBH = GERMS::dbconnect();
my $DB_PARSER = DateTime::Format::DBI->new($DBH);
my $TIP = GERMS::get_browser_tip($runID, $DBH, $USE_DB);
$inputs->{TIP} = $TIP;
if ($verbose) {
  print "RunID: $runID\n";
  print "TIP: $TIP\n";
  if ($USE_DB) {
    print "Doing database operations (-db is on)\n";
  } else {
    print "No database operations (no -db flag)\n";
  }
}
foreach $ext (@extensions) {
  $files_data = ();
  if (defined $inputs->{$ext}) {
    $files_data->{TIP} = $TIP;
    $files_data->{Type} = $ext;
    $files_data->{Filename} = $inputs->{$ext};
    $files_data->{Filename} =~ s/^$OUTBREAK_BASE\///;
    $files_data->{MD5} = $inputs->{$ext . "_MD5"};
    $files_data->{DateStamp} = $dates->{$ext};
    $status = GERMS::do_db_withSourceFile($files_data, "Files", $force, $DBH, $USE_DB);
    if ($verbose) {
      print "  File: $inputs->{$ext} ($dates->{$ext}) $status\n";
    }
  }
}
$status = ();

# SRST2 output gives these columns
# 0	Sample
# 1	DB
# 2	gene
# 3	allele
# 4	coverage
# 5	depth
# 6	diffs
# 7	uncertainty
# 8	divergence
# 9	length
# 10	maxMAF
# 11	clusterid
# 12	seqid
# 13	annotation
# all these should match database fields
# resistances have one extra field which is class
#
# MLST has MLST at the beginning, and
# mismatches, uncertainty, depth, maxMAF at the end
# will join everything else with dashes for profile
#
if (defined $inputs->{'srst2.gz'}) {
  foreach $i ("Resistance", "Genes") {
    foreach $j ("INSERT", "UPDATE", "NOTHING", "ERROR") {
      $status->{$i}->{$j} = 0;
    }
  }
  $mlst_data->{MLST} = "";
  open (SRST2, "zcat $inputs->{'srst2.gz'} |");
  while (<SRST2>) {
    chomp;
    if (s/^# // && $mlst_data->{MLST} eq "") {
      @f = split /\t/, $_;
      die "Mismatch run ID: expecting $runID but got $f[0] in $inputs->{'srst2.gz'}\n" if $runID ne $f[0];
      $mlst_data->{MLST} = $f[1];
      if (scalar(@f) > 6) {
        $mlst_data->{Profile} = join ("-", @f[2..($#f-4)]);
      }
      $mlst_data->{maxMAF} = pop @f;
      $mlst_data->{depth} = pop @f;
      $mlst_data->{uncertainty} = pop @f;
      $mlst_data->{mismatches} = pop @f;
      $mlst_data->{TIP} = $TIP;
      $mlst_data->{SOURCE} = $SOURCE;
      $mlst_data->{MLSTDatabase} = $MLSTDatabase;
      $mlst_data->{SourceFileType} = 'srst2.gz';
      $mlst_data->{SourceFileMD5} = $inputs->{'srst2.gz_MD5'};
      $status->{MLST} = GERMS::do_db_withSourceFile($mlst_data, "MLST", $force, $DBH, $USE_DB);
      next;
    }
    next if /^$/;
    if (/^>/) {
      # fasta sequence
      <SRST2>;
      next;
    }
    next if !/\t/;
    @f = split /\t/, $_;
    next if $#f == 0;
    if (defined $f[13] && ($f[13] =~ /DB:ARGannot/ || $f[13] =~ /DB:\s+/) ||
       (defined $f[1] && $f[1] =~ /ARGannot/)) {
      $resistance_data = ();
      @g = split /_/, $f[2];
      $class = pop @g;
      if ($class eq "Agly") { $class = "AGly"; }
      if ($class eq "AGlyFlqn") { $class = "AGly"; }
      if ($class eq "PheCmlA5") { $class = "Phe"; }
      if (!defined $classes->{$class} && !/['"]/) {
        print STDERR "Found unknown class $class in $ARGV[0] - skipping this...\n";
      } else {
        $resistance_data->{DB} = $f[1];
        $resistance_data->{gene} = $f[2];
        $resistance_data->{allele} = $f[3];
        $resistance_data->{coverage} = $f[4];
        $resistance_data->{depth} = $f[5];
        $resistance_data->{diffs} = $f[6];
        $resistance_data->{uncertainty} = $f[7];
        $resistance_data->{divergence} = $f[8];
        $resistance_data->{length} = $f[9];
        $resistance_data->{maxMAF} = $f[10];
        $resistance_data->{clusterid} = $f[11];
        $resistance_data->{seqid} = $f[12];
        $resistance_data->{annotation} = $f[13];
        $resistance_data->{class} = $class;
        $resistance_data->{TIP} = $TIP;
        $resistance_data->{SOURCE} = $SOURCE;
        $resistance_data->{SourceFileType} = 'srst2.gz';
        $resistance_data->{SourceFileMD5} = $inputs->{'srst2.gz_MD5'};
        $status->{Resistance}->{GERMS::do_db_withSourceFile($resistance_data, "Resistance", $force, $DBH, $USE_DB)}++;
      }
    } else {
      # everything else except resistances
      $genes_data = ();
      $genes_data->{DB} = $f[1];
      $genes_data->{gene} = $f[2];
      $genes_data->{allele} = $f[3];
      $genes_data->{coverage} = $f[4];
      $genes_data->{depth} = $f[5];
      $genes_data->{diffs} = $f[6];
      $genes_data->{uncertainty} = $f[7];
      $genes_data->{divergence} = $f[8];
      $genes_data->{length} = $f[9];
      $genes_data->{maxMAF} = $f[10];
      $genes_data->{clusterid} = $f[11];
      $genes_data->{seqid} = $f[12];
      $genes_data->{annotation} = $f[13];
      $genes_data->{TIP} = $TIP;
      $genes_data->{SOURCE} = $SOURCE;
      $genes_data->{SourceFileType} = 'srst2.gz';
      $genes_data->{SourceFileMD5} = $inputs->{'srst2.gz_MD5'};
      $status->{Genes}->{GERMS::do_db_withSourceFile($genes_data, "Genes", $force, $DBH, $USE_DB)}++;
    }
  }
  close SRST2;
  if ($verbose) {
    print "SRST2 file ($inputs->{'srst2.gz'}) done\n";
    print "  MLST: $status->{MLST}\n";
    foreach $i ("Resistance", "Genes") {
      print "  $i: INSERT for $status->{$i}->{INSERT}, UPDATE for $status->{$i}->{UPDATE}, NOTHING for $status->{$i}->{NOTHING} line(s)";
      if ($USE_DB) {
        print " with $status->{$i}->{ERROR} Errors\n";
      } else {
        print "\n";
      }
    }
  }
} # SRST2

# Assembly file
my $assembly_data;
if (defined $inputs->{'tgz'}) {
  my $tar = Archive::Tar->new;
  $tar->read($inputs->{'tgz'});

  my @files = $tar->get_files;
  my $summary = -1;
  my $log = -1;
  my @s;
  foreach $i (0..$#files) {
    if ($files[$i]->name =~ /-log.txt/) {
      $log = $i;
    } elsif ($files[$i]->name =~ /-summary.txt/) {
      $summary = $i;
    }
  }
  if ($log >= 0) {
    @s = split /\n/, $files[$log]->get_content;
    if ($s[0] =~ /^# GERMS-wgs.pl ## (.*)/) {
      $assembly_data->{date_run} = $1;
    }
  }
  if ($summary >= 0) {
    @s = split /\n/, $files[$summary]->get_content;
    foreach $i (0..$#s) {
      chomp $s[$i];
      if ($s[$i] =~ /Your original command line was (.*)/) {
        $assembly_data->{command} = $1;
      } elsif ($s[$i] =~ /STARTED (.*) using file.* in (.*) by (.*) on (.*)/) {
        $assembly_data->{name} = $1;
        $assembly_data->{run} = $1;
        $assembly_data->{temp} = $2;
        $assembly_data->{user} = $3;
        $assembly_data->{host} = $4;
      } elsif ($s[$i] =~ /^Files:$/) {
        foreach (1,2,3) {
          $i++;
          if ($s[$i] =~ /^You can follow/) {
            $i--;
            last;
          }
          if ($s[$i] =~ /([0-9a-f]+)\s+(.*)$/i) {
            my $md5 = $1;
            my $file = $2;
            if ($file =~ /GERMS-wgs.pl$/) {
              $assembly_data->{GERMS_MD5} = $md5;
            } elsif (!defined $assembly_data->{R1}) {
              $assembly_data->{R1} = $file;
              $assembly_data->{R1_MD5} = $md5;
            } else {
              $assembly_data->{R2} = $file;
              $assembly_data->{R2_MD5} = $md5;
            }
          }
        }
        if (!defined $assembly_data->{R2}) {
          $assembly_data->{R2} = "None";
          $assembly_data->{R2_MD5} = "None";
        }
      } elsif ($s[$i] =~ /^velveth kmer: (\d+)/) {
        $assembly_data->{velveth_k} = $1;
      } elsif ($s[$i] =~ /^velvetg exp_cov: (\d+)/) {
        $assembly_data->{velvetg_exp_cov} = $1;
      } elsif ($s[$i] =~ /^velvetg cov_cutoff: (\d+)/) {
        $assembly_data->{velvetg_cov_cutoff} = $1;
      } elsif ($s[$i] =~ /== Final Summary ==/) {
        $i++;
        while ($i <= $#s) {
          if ($s[$i] =~ /Number of .*-end reads: (\d.*)$/) {
            # sometimes get large numbers in scientific notation - add 0 to fix
            $assembly_data->{num_reads} = $1 + 0;
          } elsif ($s[$i] =~ /Read length .*: (\d+)/) {
            $assembly_data->{read_length} = $1;
          } elsif ($s[$i] =~ /Assembly: (.*)/) {
            $assembly_data->{assembly} = $1;
          } elsif ($s[$i] =~ /Scaffolding: (.*)/) {
            $assembly_data->{scaffolding} = $1;
          } elsif ($s[$i] =~ /Finishing: (.*)/) {
            $assembly_data->{finishing} = $1;
          } elsif ($s[$i] =~ /number_contigs: (\d+)/) {
            $assembly_data->{num_contigs} = $1;
          } elsif ($s[$i] =~ /total_length: (\d+)/) {
            $assembly_data->{total_length} = $1;
          } elsif ($s[$i] =~ /avg_length: (\d+)/) {
            $assembly_data->{average_length} = $1;
          } elsif ($s[$i] =~ /max_length: (\d+)/) {
            $assembly_data->{max_length} = $1;
          } elsif ($s[$i] =~ /min_length: (\d+)/) {
            $assembly_data->{min_length} = $1;
          } elsif ($s[$i] =~ /N50_length: (\d+)/) {
            $assembly_data->{N50} = $1;
          } elsif ($s[$i] =~ /N50_number: (\d+)/) {
            $assembly_data->{N50_number} = $1;
          } elsif ($s[$i] =~ /N90_length: (\d+)/) {
            $assembly_data->{N90} = $1;
          } elsif ($s[$i] =~ /N90_number: (\d+)/) {
            $assembly_data->{N90_number} = $1;
          } elsif ($s[$i] =~ /Reference sequence: (.*)/) {
            $assembly_data->{closest_reference} = $1;
          } elsif ($s[$i] =~ /Reads used for assembly: (.*)/) {
            $assembly_data->{assembly_reads} = $1;
          } elsif ($s[$i] =~ /Kraken classification: (.*)/) {
            $assembly_data->{Kraken_classification} = $1;
          }
          $i++;
        }
      }
    }
  }
  if (!defined $assembly_data->{assembly_reads} && defined $assembly_data->{num_reads}) {
    $assembly_data->{assembly_reads} = $assembly_data->{num_reads};
  }
  $assembly_data->{TIP} = $TIP;
  $assembly_data->{SourceFileType} = 'tgz';
  $assembly_data->{SourceFileMD5} = $inputs->{'tgz_MD5'};
  $status = GERMS::do_db_withSourceFile($assembly_data, "Assembly", $force, $DBH, $USE_DB);
  print "  Assembly: $status\n" if $verbose;
} # Assembly

if ($push_s3) {
  if ($delete) {
    $S3_OP = "mv";
  } else {
    $S3_OP = "cp";
  }
  # get all the info
  if ($reference_fna ne "") {
    $reference_fna = File::Basename::basename($reference_fna);
  }
  if ($reference_fna ne "") {
    $sql = "SELECT ReferenceFile, DefaultReference FROM ReferenceGenomes WHERE Species = ?";
  } else {
    $sql = "SELECT ReferenceFile, DefaultReference FROM ReferenceGenomes WHERE DefaultReference = 1 AND Species = ?";
  }
  my $sth = $DBH->prepare($sql);
  $S3_SPECIES = $MLSTDatabase;
  $sth->execute($MLSTDatabase);
  while (@f = $sth->fetchrow_array) {
    if ($reference_fna eq File::Basename::basename($f[0]) ||
         ($reference_fna eq "" && $f[1] == 1)) {
      @g = split /\//, File::Basename::dirname($f[0]);
      $S3_REFDIR = pop @g;
      my $species_check = pop @g;
      if (defined $species_check && length($species_check) && $species_check ne $S3_SPECIES) {
        # whitelist this for now - Ccoli uses Cjejuni MLST and reference
        if ($species_check = 'Cjejuni' && $S3_SPECIES eq 'Ccoli') {
          $S3_SPECIES = $species_check;
        }
      }
    }
  }

  # secondary files go to web browser
  # These are .resistance.fasta, .gbk, .gcov.png
  # do these first because .gocv.png might get deleted (moved) below
  if (length($S3_SPECIES)) {
    $secondary_S3 .= uc($S3_SPECIES);
    foreach $dir (@data_dirs) {
      if (-f "$dir/$runID.gcov.png") {
        $secondary = File::Spec->rel2abs("$dir/$runID.gcov.png");
        $secondary_MD5 = `md5sum $secondary | cut -f1 -d' '`;
        chomp($secondary_MD5);
        $s3_command = "aws --profile $secondary_PROFILE s3 cp $secondary $secondary_S3/coverage/ --metadata MD5sum=$secondary_MD5";
        if ($verbose) {
          print "Running: $s3_command\n";
        }
        $j = `$s3_command`;
        if ($? != 0) {
          print STDERR "error on $secondary_S3 on $secondary\n";
        }
      }
      if (-f "$dir/$runID.resistance.fasta") {
        $secondary = File::Spec->rel2abs("$dir/$runID.resistance.fasta");
        $secondary_MD5 = `md5sum $secondary | cut -f1 -d' '`;
        chomp($secondary_MD5);
        $s3_command = "aws --profile $secondary_PROFILE s3 cp $secondary $secondary_S3/resistance/ --metadata MD5sum=$secondary_MD5";
        if ($verbose) {
          print "Running: $s3_command\n";
        }
        $j = `$s3_command`;
        if ($? != 0) {
          print STDERR "error on $secondary_S3 on $secondary\n";
        }
      }
      if (-f "$dir/$runID.gbk") {
        $secondary = File::Spec->rel2abs("$dir/$runID.gbk");
        $secondary_MD5 = `md5sum $secondary | cut -f1 -d' '`;
        chomp($secondary_MD5);
        $s3_command = "aws --profile $secondary_PROFILE s3 cp $secondary $secondary_S3/assembly/ --metadata MD5sum=$secondary_MD5";
        if ($verbose) {
          print "Running: $s3_command\n";
        }
        $j = `$s3_command`;
        if ($? != 0) {
          print STDERR "error on $secondary_S3 on $secondary\n";
        }
      }
    }
  }

  # main archive - we are capturing .lofreq.gz, .gcov.gz, .gcov.png, .srst2.gz, .tgz
  # note we are not doing .lofreq.gz.tbi, since that's quick to generate
  if (length($S3_REFDIR)) {
    # check for existence, timestamp, MD5 if present to decide on upload
    foreach $ext (keys %$inputs) {
      next if $ext =~ /_MD5$/;
      next if $ext eq "TIP";
      if (defined $inputs->{$ext . "_MD5"}) {
        $i = "--metadata MD5sum=" . $inputs->{$ext . "_MD5"};
      } else {
        $i = "";
      }
      $s3_object = "s3://$S3_BUCKET/$S3_SPECIES/$S3_REFDIR/" . File::Basename::basename($inputs->{$ext});
      $s3_command = "aws s3 ls $s3_object";
      $j = `$s3_command`;
      if ($? == 0) {
        $s3_temp = "$S3_SPECIES/$S3_REFDIR/" . File::Basename::basename($inputs->{$ext});
        $s3_command = "aws s3api head-object --bucket $S3_BUCKET --key $s3_temp";
        $j = `$s3_command`;
        $j = JSON::decode_json($j);
        if (defined $j->{LastModified}) {
          $s3_date = DateTime::Format::DateParse->parse_datetime($j->{LastModified});
        }
        if (defined $j->{Metadata}->{md5sum}) {
          $s3_MD5 = $j->{Metadata}->{md5sum};
        }
      }
      if ($force eq "S3" ||
           ( (!defined($s3_MD5) || $s3_MD5 ne $inputs->{$ext . "_MD5"}) &&
             (!defined($s3_date) || ($s3_date < $dates->{$ext}) ) )
         ) {
        $s3_command = "aws s3 $S3_OP $inputs->{$ext} s3://$S3_BUCKET/$S3_SPECIES/$S3_REFDIR/ $i";
        if ($USE_DB) {
          if ($verbose) {
            print "Running: $s3_command\n";
          }
          $j = `$s3_command`;
          if ($? != 0) {
            print STDERR "error on S3 $S3_OP on $inputs->{$ext}\n";
          }
        } else {
          if ($verbose) {
            print "Would run: $s3_command\n";
          }
        }
      }
    }
  }
}
