#!/usr/bin/perl
#
# some useful functions
#

use warnings;
use Cwd;
use File::Temp;
use File::Spec;
use File::Basename;
use File::Copy;
use File::Type;
use Compress::Zlib;
use JSON -support_by_pp;
use Data::Dumper;
use DBI;
use DateTime::Format::DBI;
use DateTime::Format::DateParse;
use strict;

package GERMS;
require Exporter;

use vars qw(@ISA @EXPORT @EXPORT_OK);
use vars qw(%orgcode_translate);

@ISA = ("Exporter");
@EXPORT = ();
@EXPORT_OK = ();

sub dbconnect {
  my ($database, $host, $user, $pass, $port) = @_;
  # preconfig should give us a hash reference here as $config
  my @paths = ("/etc/GERMS.conf",
               "/usr/etc/GERMS.conf",
               "/usr/local/etc/GERMS.conf",
               "~/.GERMS.conf");
  use vars qw($config);
  my $config_file;

  # default values first
  foreach my $config_file (@paths) {
    if (-f $config_file) {
      do $config_file;
    }
  }

  if (defined $config->{DEFAULT} && (!defined $database || $database eq "")) {
    $database = $config->{DEFAULT};
  }
  if (defined $config->{$database}) {
  $database = 'positive_selection' if !defined $database || $database eq "";
    $host = $config->{$database}->{HOST} if !defined $host || $host eq "";
    $user = $config->{$database}->{USER} if !defined $user || $user eq "";
    $pass = $config->{$database}->{PASS} if !defined $pass || $pass eq "";
    $port = $config->{$database}->{PORT} if !defined $port || $port eq "";
  } else {
    $host = "localhost" if !defined $host || $host eq "";
    $port = 3306 if !defined $port || $port eq "";
  }

  my $dbh;
  $dbh = DBI->connect('DBI:mysql:database='.$database.';host='.$host.';port='.$port, $user, $pass, { LongReadLen => 100000000 });
  return $dbh;
}

sub make_key {
  my @k = @_;
  return (join ("~~", @k));
}

sub break_key {
  my ($k) = @_;
  return (split /~~/, $k);
}

# sort unique, be smart about numbers and non numbers
sub sortu {
  my @a = @_;
  my @defined = ();
  my @undefined = ();
  my $number = 1;
  foreach my $i (0..$#a) {
    if (!defined $a[$i]) {
      push @undefined, $i;
    }
    push @defined, $i;
    if (!isfloat($a[$i])) {
      $number = 0;
    }
  }
  if ($number) {
    @a = sort { $a <=> $b } @a[@defined];
  } else {
    @a = sort { $a cmp $b } @a[@defined];
  }
  foreach my $i (reverse 1..$#a) {
    if ($number) { 
      if ($a[$i] == $a[$i-1]) { splice @a, $i, 1; }
    } else {
      if ($a[$i] eq $a[$i-1]) { splice @a, $i, 1; }
    }
  }
  foreach my $i (@undefined) {
    push @a, undef;
  }
  return @a;
}

# take 1,2,3,5-9,14..18,22 and give array with all numbers
sub parse_list {
  my ($list) = @_;
  my @numbers = ();
  my @temp = ();
  chomp $list;
  my @atoms = split /,/, $list;
  foreach my $i (0..$#atoms) {
    if (isfloat($atoms[$i])) {
      push @numbers, $atoms[$i];
    } elsif ($atoms[$i] =~ /\.\./) {
      @temp = split /\.\./, $atoms[$i];
      if (isint($temp[0]) && isint($temp[1])) {
        if ($temp[0] < $temp[1]) {
          push @numbers, ($temp[0]..$temp[1]);
        } else {
          push @numbers, reverse ($temp[1]..$temp[0]);
        }
      }
    } elsif ($atoms[$i] =~ /-/) {
      @temp = split /-/, $atoms[$i];
      if (isint($temp[0]) && isint($temp[1])) {
        if ($temp[0] < $temp[1]) {
          push @numbers, ($temp[0]..$temp[1]);
        } else {
          push @numbers, reverse ($temp[1]..$temp[0]);
        }
      }
    }
  }
  return @numbers;
}

# take list of integers, give a shortened list of ranges
sub make_list {
  my (@a) = @_;
  my $ranger = "..";
  my $separater = ",";
  my $i;
  my @return;
  foreach $i (0..$#a) {
    if (!isint($a[$i])) {
      print STDERR "Noninteger $a[$i] found in make_list call\n";
      return @a;
    }
  }
  @a = sortu (@a);
  my $start;
  my $end;
  undef $start;
  undef $end;
  foreach $i (0..$#a) {
    if (!defined $start) {
      $start = $a[$i];
      $end = $a[$i];
      next;
    }
    if ($a[$i] == $end + 1) {
      $end = $a[$i];
    } else {
      if ($end == $start) {
        push @return, $start;
      } elsif ($end == $start+1) {
        push @return, $start, $end;
      } else {
        push @return, join ($ranger, $start, $end);
      }
      $start = $a[$i];
      $end = $a[$i];
    }
  }
  if (defined $start && defined $end) {
    if ($end == $start) {
      push @return, $start;
    } elsif ($end == $start+1) {
      push @return, $start, $end;
    } else {
      push @return, join ($ranger, $start, $end);
    }
  }
  return join ($separater, @return);
}

sub isfloat {
  my ($string) = @_;
  if (!defined $string) {
    return 0;
  }
  chomp $string;
  if ($string =~ /^([+-])?(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
    return 1;
  } else {
    return 0;
  }
}

sub isint {
  my ($string) = @_;
  if (!defined $string) {
    return 0;
  }
  chomp $string;
  if ($string =~ /^[+-]?\d+$/) {
    return 1;
  } else {
    return 0;
  }
}

sub align {
  # take a fasta array of protein sequence, make sure fasta names are short
  # run it through clustalw
  # return the results
  my @f = @_;
  my @return = ();
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $tempfile = $tempdir."/clustal-input";
  my $clustalout = $tempfile.".phy";
#  my $clustalout = $tempfile.".aln";
  my $namehash = to_phylip_names(\@f);
  my $in;
  open T, ">$tempfile";
  foreach my $f (@f) {
    chomp $f;
    next if $f =~ /^#/;
    next if $f =~ /^$/;
    print T $f, "\n";
  }
  close T;
  if (scalar @f > 2) {
    system "clustalw -output=PHYLIP -type=PROTEIN -infile=$tempfile > /dev/null"
;
    open T, "phy2fasta -interleaved $clustalout |";
  } else {
    open T, $tempfile;
  }
  while ($in = <T>) {
    chomp $in;
    next if $in =~ /^#/;
    next if $in =~ /^$/;
    push @return, $in;
  }
  close T;
  &from_phylip_names (\@return, $namehash);
  unlink $tempfile;
  unlink "$tempfile.dnd";
  unlink $clustalout;
  system "rm -rf $tempdir";
  return @return;
}

sub dnaalign {
  # take a fasta array of dna sequence
  # run it through clustalw
  # return the results as an array
  my (@g) = @_;
  my @f = ();
  my @return = ();
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $tempfile = $tempdir."/clustal-input";
#  my $clustalout = $tempfile.".phy";
  my $clustalout = $tempfile.".aln";
  my $in;
  if (ref($g[0]) eq "ARRAY") {
    @f = @{$g[0]};
  } elsif (ref($g[0]) eq "HASH") {
    @f = hash2fasta($g[0]);
  } else {
    @f = @g;
  }
  open T, ">$tempfile";
  my $namehash = to_phylip_names(\@f);
  foreach my $f (@f) {
    chomp $f;
    next if $f =~ /^#/;
    next if $f =~ /^$/;
    print T $f, "\n";
  }
  close T;
  if (scalar @f > 2) {
#    system "clustalw -output=PHYLIP -type=DNA -infile=$tempfile > /dev/null"
#    system "clustalw -type=DNA -infile=$tempfile -gapopen=1000 > /dev/null"
    system "clustalw -type=DNA -infile=$tempfile > /dev/null"
;
    open T, "aln2phy $clustalout | phy2fasta -interleaved |";
  } else {
    open T, $tempfile;
  }
  while ($in = <T>) {
    chomp $in;
    next if $in =~ /^#/;
    next if $in =~ /^$/;
    push @return, $in;
  }
  close T;
  &from_phylip_names (\@return, $namehash);
  unlink $tempfile;
  unlink "$tempfile.dnd";
  unlink $clustalout;
  system "rm -rf $tempdir";
  return @return;
}

sub paln2daln {
  # take an aligned fasta array of aligned protein sequence
  # also need fasta array of dna sequence with the same headers
  # return a DNA version with codons kept together
  # look for the gid in the second field of the header line
  my ($pref, $dref, $errorref) = @_;
  $$errorref = 0;
  my @prot = @$pref;
  my @dna = @$dref;
  my (@return, $pos, %dseq, $name, $dseq);
  foreach my $i (0..$#dna) {
    chomp $dna[$i];
    if ($dna[$i] =~ /^>/) {
      $name = $dna[$i];
    } else {
      $dseq{$name} = $dna[$i];
    }
  }
  foreach my $i (0..$#prot) {
    chomp $prot[$i];
    if ($prot[$i] =~ /^>/) {
      $name = $prot[$i];
      $return[$i] = $name;
    } else {
      my $pseq = $prot[$i];
      if (!defined $dseq{$name}) {
        die "Can't find DNA sequence for $name\n";
      }
      $dseq = $dseq{$name};
      my @p = split //, $pseq;
      $return[$i] = "";
      $pos = 0;
      my $d;
      foreach my $p (@p) {
        if ($p eq '-') {
          $return[$i] .= "---";
        } else {
          $d = substr ($dseq, $pos, 3);
          if (aa($d) ne uc $p && $pos > 0 && $d =~ m/[^GATCgatc]/) {
            print STDERR "Translation error at sequence $name, pos $pos, dna $d, prot $p, trans ", aa($d), "\n";
            $$errorref += 1;
          }
          if (aa($d) eq 'X') {
            $return[$i] .= '---';
          } else {
            $return[$i] .= $d;
          }
          $pos += 3;
        }
      }
    }
  }
  return @return;
}

sub align_stripgaps {
  # take aligned fasta sequence
  # strip out gaps from alignment
  my (@s) = @_;
  my ($i, $j, $k);
  my $sref;
  my $name;
  my @a;
  my @return;
  $j = 0;
  foreach $i (0..$#s) {
    if ($s[$i] =~ /^>/) {
      $sref->[$j]->{NAME} = $s[$i];
    } else {
      $sref->[$j]->{SEQ} = $s[$i];
      $j++;
    }
  }
  foreach $j (0..$#$sref) {
    @a = split //, $sref->[$j]->{SEQ};
    foreach $i (reverse 0..$#a) {
      if ($a[$i] eq '-') {
        foreach $k (0..$#$sref) {
          substr ($sref->[$k]->{SEQ}, $i, 1) = "";
        }
      }
    }
  }
  foreach $j (0..$#$sref) {
    push @return, $sref->[$j]->{NAME};
    push @return, $sref->[$j]->{SEQ};
  }
  return @return;
}

# take orgcode, references to hashes, fill them
# hashes are:
#   gid2triv, gid2syst, gid2annot, gid2faa, gid2ffn, gid2loc, gid2org, loc2ffn
sub setup_hash {
  my ($org, $gid2triv_r, $gid2syst_r, $gid2annot_r, $gid2faa_r, $gid2ffn_r, $gid2loc_r, $gid2org_r, $loc2ffn_r) = @_;
  my (@ptt, @faa, @ffn) = ();
  my (@subseq) = ();
  my ($gi, $i, $j, @f, @g, $start, $stop, $strand, $switch, $loc, @sub);
  @ptt = `print-genome-file.pl $org -suf ptt`;
  @faa = `print-genome-file.pl $org -suf faa | chompnewline.pl`;
  @ffn = `print-genome-file.pl $org -suf ffn | chompnewline.pl`;
  foreach my $p (@ptt) {
    chomp $p;
    if ($p =~ /^\s*(\d+)\.\.(\d+)\s+/) {
      ($start, $stop) = ($1, $2);
      @f = split /\t/, $p;
      $strand = $f[1];
      $gi = $f[3];
      next if ($gi < 0 || $gi == 9999);
      if (defined $gid2triv_r) {
        $$gid2triv_r{$gi} = $f[4];	# this is usually ok
      }
      if (defined $gid2syst_r) {
        $$gid2syst_r{$gi} = $f[5];
      }
      if (defined $gid2annot_r) {
        $$gid2annot_r{$gi} = $f[8];
      }
      if (defined $gid2loc_r) {
        $$gid2loc_r{$gi} = "$org-$start..$stop";
      }
      if (defined $gid2org_r) {
        $$gid2org_r{$gi} = $org;
      }
      if ($strand =~ /-/) {
        $switch = $start;
        $start = $stop;
        $stop = $switch;
      }
      push @subseq, ">~~$gi~~\n$start..$stop\n";
    }
  }
  foreach my $s (@faa) {
    chomp $s;
    if ($s =~ /^>gi\|(\d+)\|/) {
      $gi = $1;
    } else {
      if (defined $gid2faa_r) {
        $$gid2faa_r{$gi} = $s;
      }
    }
  }
  if (defined $loc2ffn_r) {
    foreach my $s (@ffn) {
      chomp $s;
      if ($s =~ /^>/) {
        @f = split /:/, $s;
        $f[1] =~ s/\s.*//;
        $f[1] =~ s/\(//g;
        $f[1] =~ s/\)//g;
        $f[1] =~ s/c//g;
        $f[1] =~ s/\s//g;
        $f[1] =~ s/>//g;
        @g = split /[-,]/, $f[1];
        @g = sort { $a <=> $b } @g;
        $loc = "$org-$g[0]..$g[$#g]";
      } else {
        $$loc2ffn_r{$loc} = $s;
      }
    }
  }
  if (defined $gid2ffn_r) {
    my $temp = "/tmp/pull.temp.".rand();
    while (-f $temp) {
      $temp = "/tmp/pull.temp.".rand();
    }
    open SUBSEQ, "| subseq.pl $org > $temp";
    print SUBSEQ @subseq;
    close SUBSEQ;
    open SUBSEQ, $temp;
    @sub = <SUBSEQ>;
    close SUBSEQ;
    unlink ($temp);
    $gi = -1;
    foreach $j (0..$#sub) {
      chomp $sub[$j];
      if ($sub[$j] =~ /^>~~(\d+)~~$/) {
        $gi = $1;
      } elsif ($sub[$j] !~ /^>/) {	# if frameshift, use genbank ffn
        if ($gi > 0 && $$gid2loc_r{$gi} && $$loc2ffn_r{$$gid2loc_r{$gi}}
            && $$loc2ffn_r{$$gid2loc_r{$gi}} ne $sub[$j]) {
          $$gid2ffn_r{$gi} = $$loc2ffn_r{$$gid2loc_r{$gi}};
        } else {
          $$gid2ffn_r{$gi} = $sub[$j];
        }
        $gi = -1;
      }
    }
  }
}

sub translate {
  # translate sequence to aa
  # use frame if we have it
  my $sequence = shift;
  my $frame;
  if (scalar @_) {
    $frame = shift;
  }
  $sequence = uc($sequence);
  $frame = 1 if !defined $frame;
  if ($frame < 0) {
    $sequence = revcomp($sequence);
  }
  $frame = abs($frame);
  $sequence = substr($sequence, $frame-1);
  my $i;
  my $return = "";
  for ($i = 0; $i < length $sequence; $i += 3) {
    $return .= aa(substr($sequence, $i, 3));
  }
  return $return;
}

sub revcomp {
  my ($inline) = $_[0];
  my $outline = reverse ($inline);
  $outline =~ tr/ABCDGHKMNRSTVWXYabcdghkmnrstvwxy/TVGHCDMKNYSABWXRtvghcdmknysabwxr/;
  return $outline;
}

# translate codon to aa
sub aa {
  my ($d) = @_;
  $d = uc $d;
  if ($d eq "TTT") { return "F"; }
  if ($d eq "TTC") { return "F"; }
  if ($d eq "TTY") { return "F"; }
  if ($d eq "TTA") { return "L"; }
  if ($d eq "TTG") { return "L"; }
  if ($d eq "TTR") { return "L"; }
  if ($d eq "TCT") { return "S"; }
  if ($d eq "TCC") { return "S"; }
  if ($d eq "TCA") { return "S"; }
  if ($d eq "TCG") { return "S"; }
  if ($d =~ /^TC/) { return "S"; }
  if ($d eq "TAT") { return "Y"; }
  if ($d eq "TAC") { return "Y"; }
  if ($d eq "TAY") { return "Y"; }
  if ($d eq "TAA") { return "X"; }
  if ($d eq "TAG") { return "X"; }
  if ($d eq "TGT") { return "C"; }
  if ($d eq "TGC") { return "C"; }
  if ($d eq "TGY") { return "C"; }
  if ($d eq "TGA") { return "X"; }
  if ($d eq "TGG") { return "W"; }
  if ($d eq "CTT") { return "L"; }
  if ($d eq "CTC") { return "L"; }
  if ($d eq "CTA") { return "L"; }
  if ($d eq "CTG") { return "L"; }
  if ($d =~ /^CT/) { return "L"; }
  if ($d eq "YTA") { return "L"; }
  if ($d eq "YTG") { return "L"; }
  if ($d eq "YTR") { return "L"; }
  if ($d eq "CCT") { return "P"; }
  if ($d eq "CCC") { return "P"; }
  if ($d eq "CCA") { return "P"; }
  if ($d eq "CCG") { return "P"; }
  if ($d =~ /^CC/) { return "P"; }
  if ($d eq "CAT") { return "H"; }
  if ($d eq "CAC") { return "H"; }
  if ($d eq "CAY") { return "H"; }
  if ($d eq "CAA") { return "Q"; }
  if ($d eq "CAG") { return "Q"; }
  if ($d eq "CAR") { return "Q"; }
  if ($d eq "CGT") { return "R"; }
  if ($d eq "CGC") { return "R"; }
  if ($d eq "CGA") { return "R"; }
  if ($d eq "CGG") { return "R"; }
  if ($d =~ /^CG/) { return "R"; }
  if ($d eq "MGA") { return "R"; }
  if ($d eq "MGG") { return "R"; }
  if ($d eq "MGR") { return "R"; }
  if ($d eq "ATT") { return "I"; }
  if ($d eq "ATC") { return "I"; }
  if ($d eq "ATA") { return "I"; }
  if ($d eq "ATY") { return "I"; }
  if ($d eq "ATW") { return "I"; }
  if ($d eq "ATM") { return "I"; }
  if ($d eq "ATH") { return "I"; }
  if ($d eq "ATG") { return "M"; }
  if ($d eq "ACT") { return "T"; }
  if ($d eq "ACC") { return "T"; }
  if ($d eq "ACA") { return "T"; }
  if ($d eq "ACG") { return "T"; }
  if ($d =~ /^AC/) { return "T"; }
  if ($d eq "AAT") { return "N"; }
  if ($d eq "AAC") { return "N"; }
  if ($d eq "AAY") { return "N"; }
  if ($d eq "AAA") { return "K"; }
  if ($d eq "AAG") { return "K"; }
  if ($d eq "AAR") { return "K"; }
  if ($d eq "AGT") { return "S"; }
  if ($d eq "AGC") { return "S"; }
  if ($d eq "AGY") { return "S"; }
  if ($d eq "AGA") { return "R"; }
  if ($d eq "AGG") { return "R"; }
  if ($d eq "AGR") { return "R"; }
  if ($d eq "GTT") { return "V"; }
  if ($d eq "GTC") { return "V"; }
  if ($d eq "GTA") { return "V"; }
  if ($d eq "GTG") { return "V"; }
  if ($d =~ /^GT/) { return "V"; }
  if ($d eq "GCT") { return "A"; }
  if ($d eq "GCC") { return "A"; }
  if ($d eq "GCA") { return "A"; }
  if ($d eq "GCG") { return "A"; }
  if ($d =~ /^GC/) { return "A"; }
  if ($d eq "GAT") { return "D"; }
  if ($d eq "GAC") { return "D"; }
  if ($d eq "GAY") { return "D"; }
  if ($d eq "GAA") { return "E"; }
  if ($d eq "GAG") { return "E"; }
  if ($d eq "GAR") { return "E"; }
  if ($d eq "GGT") { return "G"; }
  if ($d eq "GGC") { return "G"; }
  if ($d eq "GGA") { return "G"; }
  if ($d eq "GGG") { return "G"; }
  if ($d =~ /^GG/) { return "G"; }
  return "X";
}

sub randomize_array {
  my (@a) = @_;
  my @index;
  foreach my $i (0..$#a) {
    $index[$i] = rand();
  }
  return (@a[sort { $index[$a] <=> $index[$b] } (0..$#a)]);
}

sub to_phylip_names {
  # return fasta dna sequence with names that are short enough (< 10 char) for 
  # phylip format
  # we're going to use numbers preceded by a hopefully unique string
  # this gives us a capacity of 10^6 or so
  # return a hash reference that can be used to convert back also
  # try to take both arrays and hashes
  # usage: $name_hash = to_phylip_names(\@dna_array);
  # or:    $name_hash = to_phylip_names(\%dna_hash);
  my ($dna_ref) = @_; 
  my $name;
  my $i = 0;
  my $order = 0;
  my $j;
  my $string;
  if (ref($dna_ref) eq "ARRAY") {
    foreach $j (0..$#{$dna_ref}) {
      $i++ if $dna_ref->[$j] =~ /^>/;
    }
    $order = length($i);
    $i = 0;
    foreach $j (0..$#{$dna_ref}) {
      chomp $dna_ref->[$j]; 
      if ($dna_ref->[$j] =~ /^>(.*)$/) {
        $i++;
        $string = name_key($i, $order);
        $name->{$string} = $1;
        $dna_ref->[$j] = ">$string";
      }
    }
  } elsif (ref($dna_ref) eq "HASH") {
    $order = length(scalar(keys(%$dna_ref)));
    foreach $j (keys %$dna_ref) {
      $i++;
      $string = name_key($i, $order);
      $name->{$string} = $j;
      $dna_ref->{$string} = $dna_ref->{$j};
      delete $dna_ref->{$j};
    }
  }
  return ($name); 

  sub name_key {
    my ($i, $order) = @_;
    return sprintf("__P__%." . $order . "d", $i);
  }
}

sub from_phylip_names {
  # take fasta sequence named by to_phylip_names and the name hash
  # change the sequences back
  # usage: &from_phylip_names(\@dna_array, $name_hash);
  # or:    &from_phylip_names(\%dna_hash, $name_hash);
  # actually hope the names are unique so we can take anything here, like trees
  # so:    $string1 = from_phylip_names($string2, $name_hash);
  my ($dna_ref, $name) = @_;
  my $j = 0;
  if (ref($dna_ref) eq "ARRAY") {
    foreach $j (0..$#{$dna_ref}) {
      chomp $dna_ref->[$j];
      if ($dna_ref->[$j] =~ /^>(.*)$/ && defined $name->{$1}) {
        $dna_ref->[$j] = ">" . $name->{$1};
      }
    }
  } elsif (ref($dna_ref) eq "HASH") {
    foreach $j (keys %$dna_ref) {
      if (defined $name->{$j}) {
        $dna_ref->{$name->{$j}} = $dna_ref->{$j};
        delete $dna_ref->{$j};
      }
    }
  } else {
    foreach $j (keys %$name) {
      $dna_ref =~ s/$j/$name->{$j}/;
    }
    return $dna_ref;
  }
}

sub dna_distances {
  my @dna = @_;
  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  my $current = `pwd`;
  chomp $current;
  chdir ($tempdir);
  my $control = "phylip.cmd";
  my $dnafile = "phylip.dna";
  open DNA, "| fasta2phy.pl > $dnafile";
  my $namehash = to_phylip_names (\@dna);
  print DNA join ("\n", @dna);
  close DNA;
  open CMD, ">$control";
  print CMD "$dnafile\n";
  print CMD "l\n";      # lower triangular output
  print CMD "y\n";      # start the program
  close CMD;

  if (-f "/usr/bin/phylip") {
    system "phylip dnadist < $control > /dev/null";
  } else {
    system "dnadist < $control > /dev/null";
  }

  open O, "outfile";
  my @out = <O>;
  close O;
  my $i;
  my $j;
  my @f;
  my $dist;
  my @name;
  foreach $i (1..$#out) {       # first line should be # of sequences
    @f = split /\s+/, $out[$i];
    $name[$i] = $f[0];
    foreach $j (1..$#f) {       # make a symmetric distance hash using original names
      $dist->{$namehash->{$f[0]}}->{$namehash->{$name[$j]}} = $f[$j];
      $dist->{$namehash->{$name[$j]}}->{$namehash->{$f[0]}} = $f[$j];
    }
  }

  chdir $current;
  system "rm -rf $tempdir";

  return $dist;
}

sub array_mean {
  # we're going to bascially throw out non number values
  my @a = @_;
  my ($t, $mean, $count);
  $mean = 0;
  $count = 0;
  foreach $t (@a) {
    if (isfloat($t)) {
      $mean += $t;
      $count++;
    }
  }
  if (!$count) { return 0; }
  return ($mean/$count);
}

sub array_median {
  # again, throwing out non number values
  my @a = @_;
  my ($t, $i, $j, @b);
  foreach $t (@a) {
    if (isfloat($t)) {
      push @b, $t;
    }
  }
  @b = sort { $a <=> $b } @b;
  return $b[0] if scalar @b == 1;
  return undef if !scalar @b;
  $i = int $#b/2;
  $j = int ($#b+1)/2;
  return ($b[$i] + $b[$j])/2;
}

sub array_mode {
  my @a = @_;
  my %count;
  my $t;
  foreach $t (@a) {
    if ($count{$t}) { $count{$t}++; }
    else { $count{$t} = 1; }
  }
  my $max = 0; 
  foreach $t (keys %count) {
    if ($count{$t} > $max) { $max = $count{$t}; }
  }
  my @max;
  foreach $t (keys %count) {
    if ($count{$t} == $max) {
      push @max, $t;
    }
  }
  return (array_mean(@max));
}

sub array_stdev {
  my @a = @_;
  my $total = 0;
  my $mean = array_mean(@a);
  my $i = 0;
  foreach my $a (@a) {
    next if !isfloat($a);
    $total += ($mean - $a) ** 2;
    $i++;
  }
  if ($i > 1) {
    return ($total/($i-1)) ** 0.5;
  } else {
    return 0;
  }
}

sub isintergenic {
  # return 1 if intergenic, 0 if not
  # require a setup first, then can reuse this hash
  # assume we have Orgmap and $orgcode and $pttfile and $sequence defined
  my ($ref, $pos) = @_;
  $pos = abs($pos);
  my $i;
  my $line;
  if (!length $Orgmap::sequence) {
    &Orgmap::read_sequence;
  }
  while ($pos > length $Orgmap::sequence) {
    $pos -= length $Orgmap::sequence;
  }
  if (!scalar keys %$ref) {
    foreach $i (1..length $Orgmap::sequence) {
      $ref->{$i} = 1;
    }
    open PTT, $Orgmap::pttfile;
    while ($line = <PTT>) {
      if ($line =~ /^\s*(\d+)\.\.(\d+)\s+/) {
        foreach $i ($1..$2) {
          $ref->{$i} = 0;
        }
      }
    }
    close PTT;
  }
  return $ref->{$pos};
}

sub orfupstream {
  # return genbankids we are upstream of if intergenic
  # do this in a hash with keys LEFT and RIGHT
  # return empty list if inside a gene
  # require a setup first, then can reuse this hash
  # assume we have Orgmap and $orgcode and $pttfile and $sequence defined
  my ($ref, $pos) = @_;
  $pos = abs($pos);
  my $i;
  my @f;
  my $line;
  my $r = {};
  if (!defined $Orgmap::sequence) {
    &Orgmap::read_sequence;
  }
  while ($pos > length $Orgmap::sequence) {
    $pos -= length $Orgmap::sequence;
  }
  if (!scalar keys %$ref) {
    print "  --Remake ptt hash...\n";
    foreach $i (1..length $Orgmap::sequence) {
      $ref->{$i} = 0;
    }
    open PTT, $Orgmap::pttfile;
    while ($line = <PTT>) {
      if ($line =~ /^\s*(\d+)\.\.(\d+)\s+/) {
        my ($s, $e) = ($1, $2);
        @f = split /\t/, $line;
        next if $f[3] <= 0;	# can't deal with RNAs right now
        $f[1] =~ s/\s//g;
        my $tag = $f[1].$f[3];
        foreach $i ($s..$e) {
          $ref->{$i} = $tag;
        }
      }
    }
    close PTT;
    print "done\n";
  }

  if (!$ref->{$pos}) {	# if we're intergenic
    # forward direction, looking for genes in + direction
    foreach $i ($pos..length($Orgmap::sequence)) {
      if ($ref->{$i}) {
        if (substr($ref->{$i}, 0, 1) eq '+') {
          $r->{RIGHT} = abs $ref->{$i};
        }
        last;
      }
    }
    # going backwards, look for genes in - direction
    foreach $i (reverse 1..$pos) {
      if ($ref->{$i}) {
        if (substr($ref->{$i}, 0, 1) eq '-') {
          $r->{LEFT} = abs $ref->{$i};
        }
        last;
      }
    }
  }
  return $r;
}

sub opengri {
  # take a filename
  # get a GRI handle
  # take care of some options also
  my ($filename, $options) = @_;
  my $o;
  if ($filename !~ /\.ps$/) {
    $filename .= ".ps";
  }
  open GRI, "| gri -batch -nowarn_offpage -output $filename";
  foreach $o (keys %$options) {
  }
  return *GRI;
}

sub ns_setup {
  # read ptt file, etc
  my ($pttfile, $fnafile) = @_;
  my ($s, $e, $d, $gi, @f);
  my $ptt;
  my $sequence = "";
  my $line;
  if (-f $pttfile) {
    open PTT, $pttfile;
    while ($line = <PTT>) {
      if ($line =~ /^\s*(\d+)\.\.(\d+)\s+/) {
        ($s, $e) = ($1, $2);
        chomp $line;
        @f = split /\t/, $line;
        $gi = $f[3];
        $f[1] =~ s/\s//g;
        $ptt->{$gi}->{STRAND} = $f[1];
        $ptt->{$gi}->{START} = $s;
        $ptt->{$gi}->{END} = $e;
        $ptt->{$gi}->{TRIVIAL} = $f[4];
        $ptt->{$gi}->{SYST} = $f[5];
        $ptt->{$gi}->{CODE} = $f[6];
        $ptt->{$gi}->{ANNOTATION} = $f[8];
      }
    }
    close PTT;
  }
  if (-f $fnafile) {
    open FNA, $fnafile;
    while ($line = <FNA>) {
      next if $line =~ /^>/;
      chomp $line;
      $sequence .= $line;
    }
  }
  return ($ptt, $sequence);
}

sub ns_setup_gbk {
  my ($gbkfile) = @_;
  my $gbk = read_gbk($gbkfile);
  my $sequence = $gbk->{__SEQUENCE__};
  my $ptt;
  my $locus;
  my $type;
  foreach $locus (keys %$gbk) {
    next if $locus =~ /__\w+__/;
    next if !defined $gbk->{$locus}->{gene};
    foreach $type (keys %{$gbk->{$locus}}) {
      if ($type eq "gene") {
        $ptt->{$locus}->{STRAND} = $gbk->{$locus}->{gene}->{STRAND};
        $ptt->{$locus}->{START} = $gbk->{$locus}->{gene}->{START};
        $ptt->{$locus}->{END} = $gbk->{$locus}->{gene}->{END};
        $ptt->{$locus}->{COMPLETE} = $gbk->{$locus}->{gene}->{COMPLETE};
        $ptt->{$locus}->{SYST} = $locus;
        if (defined $gbk->{$locus}->{gene}->{gene}) {
          $ptt->{$locus}->{TRIVIAL} = $gbk->{$locus}->{gene}->{gene};
        } else {
          $ptt->{$locus}->{TRIVIAL} = $locus;
        }
      } else {
        $ptt->{$locus}->{CODE} = "";
        $ptt->{$locus}->{TYPE} = $type;
        if (defined $gbk->{$locus}->{$type}->{product}) {
          $ptt->{$locus}->{ANNOTATION} = $gbk->{$locus}->{$type}->{product};
        }
      }
    }
    if (!defined $ptt->{$locus}->{TYPE}) {
      $ptt->{$locus}->{TYPE} = "Unknown";
    }
    if (!defined $ptt->{$locus}->{ANNOTATION}) {
      $ptt->{$locus}->{ANNOTATION} = "";
    }
  }
  return ($ptt, $sequence, $gbk);
}

sub read_gbk {
  # read it all from the genbank file
  my ($gbkfile) = @_;
  my $gbk;
  my $header;
  my $gbkline;
  my $margin;
  my $margin_length = 12;
  my $rest;
  my $i;
  my $mode = "";
  my $submode = "";
  my @f;
  open GBK, $gbkfile;
  # read in header lines
  while ($gbkline = <GBK>) {
    last if $gbkline =~ /^FEATURES\s/;
    chomp $gbkline;

    $margin = substr($gbkline, 0, $margin_length);
    $rest = substr($gbkline, $margin_length);

    $margin =~ s/\s*$//;
    if (length($margin) > 0) {
      if ($margin =~ s/^\s+// && length($mode)) {
        $submode = $margin;
        if ($mode eq "SOURCE" && $submode eq "ORGANISM") {
          $mode = "ORGANISM";
          $submode = "";
        }
      } else  {
        $mode = $margin;
        $submode = "";
      }
    }
    if ($mode eq "LOCUS") {
      # LOCUS   SCU49845   5028 bp   DNA    PLN   21-JUN-1999
      if ($rest =~ /\s*(\S+.*)\s+(\d+\s\w\w)\s+(\S+.*)\s+([A-Z]{3})\s+(\d\d-\w\w\w-\d\d\d\d)/) {
        $header->{$mode}->{Name} = $1;
        $header->{$mode}->{Length} = $2;
        $header->{$mode}->{Type} = $3;
        $header->{$mode}->{Division} = $4;
        $header->{$mode}->{Date} = $5;
      }
    } elsif ($mode eq "DBLINK") {
      @f = split /\: /, $rest;
      if ($#f == 1) {
        $header->{$mode}->{$f[0]} = $f[1];
      } else {
        $header->{$mode}->{text} = $rest;
      }
    } elsif ($mode eq "REFERENCE") {
      if ($rest =~ /^\s*(\d+)\s+\(/) {
        $i = $1;
        next;
      } else {
        if (defined $submode && length $submode) {
          $header->{$mode}->[$i]->{$submode} .= " $rest";
        } else {
          $header->{$mode}->[$i]->{$submode} = "$rest";
        }
      }
    } else {
      if (length $submode) {
        if (defined $header->{$mode}->{$submode}) {
          $header->{$mode}->{$submode} .= " $rest";
        } else {
          $header->{$mode}->{$submode} = "$rest";
        }
      } else {
        if (defined $header->{$mode}) {
          $header->{$mode} .= " $rest";
        } else {
          $header->{$mode} = "$rest";
        }
      }
    }
  }
  $gbk->{__HEADER__} = $header;

  # read in the rest
  $mode = "";
  $submode = "";
  my $text = "";
  my $temp = ();
  my $locus_tag = "";
  my $region;
  my $strand;
  my $complete;
  my $s;
  my $e;
  while ($gbkline = <GBK>) {
    last if $gbkline =~ /^ORIGIN/;
    chomp $gbkline;
    if ($gbkline =~ s/^\s+\/(\S+)="//) {
      $submode = $1;
      $text = $gbkline;
      while ($text !~ s/"$//) {
        $gbkline = <GBK>;
        chomp $gbkline;
        $gbkline =~ s/^\s+//;
        $text .= $gbkline;
      }
      if (length($mode)) {
        if ($mode eq "source") {
          $gbk->{__SOURCE__}->{$submode} = $text;
        } else {
          $temp->{$submode} = $text;
        }
      }
      if ($submode eq "locus_tag") {
        $locus_tag = $text;
      }
      $text = "";
      $submode = "";
    } else {
      if ($gbkline =~ /\d+\.\.\d+/) {
        if (defined $locus_tag && length($mode) && $mode ne "misc_feature" && $mode ne "source") {
          $gbk->{$locus_tag}->{$mode} = $temp;
          $temp = ();
          $mode = "";
          $locus_tag = "";
        }
        $gbkline =~ s/^\s+//;
        @f = split /\s+/, $gbkline;
        $region = pop @f;
        $mode = join (" ", @f);
        if ($region =~ /<\d/ || $region =~ /\d>/) {
          $complete = 0;
        } else {
          $complete = 1;
        }
        if ($region =~ /complement/) {
          $strand = "-";
        } else {
          $strand = "+";
        }
        $region =~ /(\d+)\.\.(\d+)/;
        ($s, $e) = ($1, $2);
        if ($mode eq "source") {
          $gbk->{__SOURCE__}->{START} = $s;
          $gbk->{__SOURCE__}->{END} = $e;
        } else {
          $temp->{START} = $s;
          $temp->{END} = $e;
          $temp->{COMPLETE} = $complete;
          $temp->{STRAND} = $strand;
        }
      }
    }
  }

  # get sequence
  while ($gbkline = <GBK>) {
    chomp $gbkline;
    $gbkline =~ s/^\s*\d+\s*//;
    $gbkline =~ s/\s//g;
    $gbk->{__SEQUENCE__} .= uc($gbkline);
  }

  return $gbk;
}

sub ns {
  # we'll return hash reference with fields:
  # INTERGENIC
  # SYNONYMOUS
  # AAPOSITION
  # ORIGINALAA
  # ORIGINALCODON
  # NEWAA
  # NEWCODON
  # CODONPOSITION
  # SYSTEMATIC
  # GID

  my ($p, $nt, $ptt, $sequence) = @_;
  my $tempnt = $nt;
  my $return = {};
  # first check if it's intergenic
  my @gi = ();
  my $gi;
  my ($seq, $coord, $pos, $oldcod, $newcod);
  my ($oldaa, $newaa);
  foreach $gi (keys %$ptt) {
    if ($p >= $ptt->{$gi}->{START} && $p <= $ptt->{$gi}->{END}) {
      push @gi, $gi;
    }
  }
  if (!scalar @gi) {
    $return->{INTERGENIC} = 1;
  } else {
    foreach $gi (@gi) {	# right now this only does the first of these
        $seq = substr($sequence, $ptt->{$gi}->{START} - 1, $ptt->{$gi}->{END} - $ptt->{$gi}->{START} + 1);
      if ($ptt->{$gi}->{STRAND} eq '+') {
        $coord = int (($p - $ptt->{$gi}->{START})/3);
        $pos = $p - $ptt->{$gi}->{START} - 3 * $coord;
        $oldcod = uc substr($seq, 3*$coord, 3);
        $newcod = $oldcod;
        substr ($newcod, $pos, 1) = uc $tempnt;
      } else {
        # the subseq function will take care of doing the reverse complement
        $seq = GERMS::revcomp($seq);
        $coord = int (($ptt->{$gi}->{END} - $p)/3);
        $pos = $ptt->{$gi}->{END} - $p - 3 * $coord;
        $oldcod = uc substr($seq, 3*$coord, 3);
        $newcod = $oldcod;
        $tempnt =~ tr/gatcGATC/ctagCTAG/;
        substr ($newcod, $pos, 1) = uc $tempnt;
      }
      $oldaa = aa($oldcod);
      $newaa = aa($newcod);
      if ($oldaa eq $newaa) {
        $return->{INTERGENIC} = 0;
        $return->{SYNONYMOUS} = 1;
        $return->{ORIGINALAA} = $oldaa;
        $return->{NEWAA} = $newaa;
        $return->{ORIGINALCODON} = $oldcod;
        $return->{NEWCODON} = $newcod;
        $return->{CODONPOSITION} = $pos + 1;
        $return->{SYSTEMATIC} = $ptt->{$gi}->{SYST};
        $return->{GID} = $gi;
        $return->{AAPOSITION} = $coord + 1;
      } else {
        $return->{INTERGENIC} = 0;
        $return->{SYNONYMOUS} = 0;
        $return->{ORIGINALAA} = $oldaa;
        $return->{NEWAA} = $newaa;
        $return->{ORIGINALCODON} = $oldcod;
        $return->{NEWCODON} = $newcod;
        $return->{CODONPOSITION} = $pos + 1;
        $return->{SYSTEMATIC} = $ptt->{$gi}->{SYST};
        $return->{GID} = $gi;
        $return->{AAPOSITION} = $coord + 1;
      }
    }
  }
  return $return;
}

sub null_bind_param {
  my ($sth, $field, $param) = @_;
  if (defined $param) {
    $sth->bind_param($field, $param);
  } else {
    $sth->bind_param($field, undef);
  }
}

sub sql_statement_replace {
  # take sql statement from dbi, i.e. with ? for values
  # substitute in values to get a valid SQL statement
  my ($sql) = shift;
  my @args = @_;
  my $num_slots = $sql =~ tr/?/?/;
  die "Value/slot mismatch, ", scalar @args, " values, $num_slots slots, statement $sql\n" if scalar @args != $num_slots;
  my $a;
  foreach $a (@args) {
    if (!defined $a) {
      $sql =~ s/\?/NULL/;
    } elsif (isfloat $a) {
      $sql =~ s/\?/$a/;
    } else {
      $sql =~ s/\?/'$a'/;
    }
  }
  return $sql;
}

sub get_root_name {
  # this is a support procedure for the do_*_tree procedures
  # if we're rooting, we need to know the name of the outgroup
  my $root = shift;
  my (@dna) = @_;
  my ($i, $j);
  $j = 0;
  foreach $i (0..$#dna) {
    if ($dna[$i] =~ s/^>//) {
      $j++;
      if ($j == $root) {
        return ($dna[$i]);
      }
    }
  }
}

sub tree_oneline {
  # this is a support procedure for the do_*_tree procedures
  # fix the trees so that they are one per line
  my (@t) = @_;
  my ($i, $j);
  my @r;
  $j = 0;
  foreach $i (0..$#t) {
    chomp $t[$i];
    if (!defined $r[$j]) {
      $r[$j] = $t[$i];
    } else {
      $r[$j] .= $t[$i];
    }
    if ($t[$i] =~ /;(?:\[\d\.\d+\])?$/) {
      $j++;
    }
  }
  return @r;
}

sub prune_tree {
  # this is a support procedure for the do_*_tree procedures
  # if we're rooting, we need to get the outgroup out of the tree to get a
  # fully rooted tree
  my $root_name = shift;
  my (@t) = @_;
  my $i;
  # the root_name should be at the top level, by itself
  # so it should precede or be followed by a comma
  foreach $i (0..$#t) {
    if (!($t[$i] =~ s/,$root_name(?::\d\.\d+)//) &&
        !($t[$i] =~ s/$root_name(?::\d\.\d+),//) ) {
      print STDERR "slchen.pm: sub prune_tree didn't find $root_name in $t[$i]\n";
    }
  }
  return @t;
}

sub do_pars_tree {
  # first parameter is whether we root or not (0 = no root, >0 = root)
  # we will root on the sequence specified by root (1-based, like phylip)
  # we take fasta sequence, it should already be aligned
  # we have to assume that the names are phylip-compatible already
  my $root = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n1000\n";	# randomize input order with $seed, 1000 reps
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnapars < $tempctl > /dev/null";
  } else {
    system "dnapars < $tempctl > /dev/null";
  }
  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub do_nj_tree {
  my $root = shift;
  my $kappa = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  # first need to calculate distance matrix
  open T, ">$tempctl";
  if ($kappa) {
    print T "t\n$kappa\n";	# set transition/transverion ratio
  }
  print T "2\n";		# no progress output
  print T "y\n";		# run the program
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnadist < $tempctl > /dev/null";
  } else {
    system "dnadist < $tempctl > /dev/null";
  }
  unlink $tempdna;
  system "mv $tempdir/outfile $tempdna";
  unlink $tempctl;

  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n";		# randomize input order with $seed
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip neighbor < $tempctl > /dev/null";
  } else {
    system "neighbor < $tempctl > /dev/null";
  }

  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub do_ml_tree {
  my $root = shift;
  my $kappa = shift;
  my @dna = @_;

  my $tempdir = File::Temp::tempdir( CLEANUP => 1 );

  my $tempdna = $tempdir."/infile";
  my $tempctl = $tempdir."/phylip.ctl";
  my $current = `pwd`;
  chomp $current;
  my $i;
  my $root_name;
  my @trees;
  my $temptree;
  my $seed;

  # make sure we control newlines properly
  foreach $i (0..$#dna) {
    chomp $dna[$i];
  }

  open T, "| fasta2phy.pl > $tempdna";
  print T join ("\n", @dna);
  close T;
  open T, ">$tempctl";
  $seed = int (rand(8192)) * 4 + 1;
  print T "j\n$seed\n10\n";	# randomize input order with $seed, 10 reps
  if ($kappa) {
    print T "t\n$kappa\n";	# set transition/transversion ratio
  }
  print T "s\ng\n";		# use better analysis, global rearrangements
  print T "2\n3\n";		# no progress output, no graphical tree
  if ($root) {
    print T "o\n$root\n";	# root the tree
  }
  print T "y\n";		# run the program
  close T;
  chdir $tempdir;
  if (-f "/usr/bin/phylip") {
    system "phylip dnaml < $tempctl > /dev/null";
  } else {
    system "dnaml < $tempctl > /dev/null";
  }

  $temptree = $tempdir."/outtree";
  open T, $temptree;
  @trees = <T>;
  close T;
  chdir $current;
  system "rm -rf $tempdir";
  @trees = tree_oneline(@trees);

  # if we're rooting, get the name of the outgroup
  # then clean out the tree to get a fully rooted tree
  if ($root) {
    $root_name = get_root_name($root, @dna);
    @trees = prune_tree($root_name, @trees);
  }

  return @trees;
}

sub combinations {
  # give all the combinations of $n elements of array @a
  # return an array ref where each element is a ref to another array with the
  # particular combination of elements of @a
  my $n = shift;
  my @a = @_;
  my $return = [];
  my $i;
  my $j;
  my $others;
  return if !scalar @a;
  if ($n == 1) {
    foreach my $i (0..$#a) {
      @{$return->[$i]} = ($a[$i]);
    }
    return $return;
  } else {
    foreach $i (0..$#a-$n+1) {
      $others = slchen::combinations($n-1, @a[$i+1..$#a]);
      foreach $j (0..$#$others) {
        push @$return, [ $a[$i], @{$others->[$j]} ];
      }
    }
    return $return;
  }
}

sub gri_boxwhisker {
  # take some parameters, return a string that has GRI code
  # we're going to use 25, 50, 75 percentile for box
  # mark the mean with a plus
  # if type = 0, whiskers at most 1.5*interquartile range
  # if type = 1, then whiskers at 5th and 95th percentile
  # mild outliers are x's - less than 3*interquartile range from box
  # extreme outliers are circles - more than 3*interquartile range from box
  # data should be an array reference
  my ($ref, $x, $width, $type, $draw_mean) = @_;
  $draw_mean = 1 if !defined $draw_mean;	# default to drawing the mean
  my @a;
  my $i;
  # check and clean data first
  foreach $i (@$ref) {
    if (isfloat($i)) {
      push @a, $i;
    }
  }
  @a = sort { $a <=> $b } @a;
  my $mean = array_mean(@a);
  my $median = array_median(@a);
  my @temp = @a[0..int($#a/2)];
  my $lower = array_median(@temp);
  @temp = @a[int(($#a+1)/2)..$#a];
  my $upper = array_median(@temp);
  my $iqr = $upper - $lower;
  my $lwhisker = $lower - 1.5*$iqr;
  my $hwhisker = $upper + 1.5*$iqr;
  if ($type) {
    $lwhisker = $a[int(0.05*scalar @a)];
    $hwhisker = $a[int(0.95*scalar @a)];
  }
  my @outlier = ();
  my @extreme = ();
  # lower and upper whisker must be at a data point
  foreach $i (@a) {
    if ($i >= $lwhisker) {
      $lwhisker = $i;
      last;
    }
    if ($i < $lwhisker && $i < $lower - 3*$iqr) {
      push @extreme, $i;
    } elsif ($i < $lwhisker) {
      push @outlier, $i;
    }
  }
  foreach $i (reverse @a) {
    if ($i <= $hwhisker) {
      $hwhisker = $i;
      last;
    }
    if ($i > $hwhisker && $i > $upper - 3*$iqr) {
      push @extreme, $i;
    } elsif ($i > $hwhisker) {
      push @outlier, $i;
    }
  }

  # generate the GRI code
  my $return = "";
  my $xleft = $x - $width/2;
  my $xright = $x + $width/2;
  $return .= "draw box $xleft $lower $xright $upper\n";
  $return .= "draw line from $xleft $median to $xright $median\n";
  if ($draw_mean) {
    $return .= "draw symbol plus at $x $mean\n";
  }
  $return .= "draw line from $xleft $lwhisker to $xright $lwhisker\n";
  $return .= "draw line from $xleft $hwhisker to $xright $hwhisker\n";
  $return .= "draw line from $x $lower to $x $lwhisker\n";
  $return .= "draw line from $x $upper to $x $hwhisker\n";
  foreach $i (@outlier) {
    $return .= "draw symbol times at $x $i\n";
  }
  foreach $i (@extreme) {
    $return .= "draw symbol circ at $x $i\n";
  }
  return $return;
}

sub poisson {
  my ($k, $lambda) = @_;
  my $fact = 1;
  foreach my $i (2..$k) {
    $fact *= $i;
  }
  return ($lambda ** $k * exp(-$lambda) / $fact);
}

sub shannon {
  my (@a) = @_;
  # take an array of things
  # return shannon entropy
  # undefined values or zero length values count as no measurements
  my %items;
  my $a;
  my $h = 0;
  my $p;
  foreach $a (@a) {
    $items{$a}++ if defined $a && length $a;
  }
  foreach $a (keys %items) {
    if (defined $items{$a} && $items{$a}) {
      $p = $items{$a} / scalar(@a);
      $h += -$p * log($p)/log(2);
    }
  }
  return $h;
}

sub moving_average {
  my $window = shift;
  my (@a) = @_;
  # return moving average
  # try to handle wrap around also
  return @a if $#a < $window;
  my @w = @a[$#a-$window+1..$#a];
  my @return = ();
  my $position = int($#a - $window/2);
  my $i = 0;
  while ($i <= $#a) {
    $return[$position] = array_mean(@w);
    $position++;
    $position -= scalar @a if $position >= scalar @a;
    push @w, $a[$i];
    shift @w;
    $i++;
  }
  return @return;
}

sub segregating_sites {
  my ($seq) = @_;
  my $nt;
  my $length;
  my $i;
  my $sites = 0;
  my $name;
  my $consensus;
  foreach $name (keys %$seq) {
    @{$nt->{$name}} = split //, $seq->{$name};
    $length = length $seq->{$name} if !defined $length;
    die "sequence length ", length $seq->{name}, " but expected $length\n" if $length != length $seq->{$name};
  }
  foreach $i (0..$length-1) {
    $consensus = "";
    foreach $name (keys %$nt) {
      $consensus = $nt->{$name}->[$i] if $consensus eq "";
      if ($nt->{$name}->[$i] ne $consensus) {
        $sites++;
        last;
      }
    }
  }
  return $sites;
}

sub average_difference {
  my ($seq) = @_;
  my @name = keys %$seq;
  my $i;
  my $j;
  my $n = 0;
  my $total = 0;
  foreach $i (0..$#name-1) {
    foreach $j ($i+1..$#name) {
      $n++;
      $total += differences($seq->{$name[$i]}, $seq->{$name[$j]});
    }
  }
  return $total/$n if $n;
}

sub differences {
  my ($a, $b) = @_;
  # take 2 sequences
  # give # of differences
  # they need to be the same length and aligned already
  my $diff = 0;
  my $i;
  my @a = split //, $a;
  my @b = split //, $b;
  die "sub differences, lengths different:\n$a\n$b\n" if scalar @a != scalar @b;
  foreach $i (0..$#a) {
    $diff++ if $a[$i] ne $b[$i];
  }
  return $diff;
}

sub Si {
  my ($seq, $i) = @_;
  # this is number of sites that occur a given number of times
  # supposed to be derived alleles but we just assume the first sequence is
  # ancestral
  #
  my $s;
  my $n = scalar keys %$seq;
  my $j;
  my $ancestral;
  my $k;
  my $name;
  my $length;
  my $Si = 0;
  foreach $name (keys %$seq) {
    $length = length $seq->{$name} if !defined $length;
    @{$s->{$name}} = split //, $seq->{$name};
  }
  foreach $k (0..$length-1) {
    $j = 0;
    $ancestral = "";
    foreach $name (keys %$seq) {
      $ancestral = $s->{$name}->[$k] if !length $ancestral;
      $j++ if $s->{$name}->[$k] ne $ancestral;
    }
    $Si++ if $j == $i;
  }
  return $Si;
}

sub fasta2hash {
  # try to chomp newlines so we can take any format
  # try to be smart about strings vs. arrays
  my @a = @_;
  my @b = ();
  my $seq;
  my $name;
  my $a;
  my $i;
  my $test;
  foreach $a (@a) {
    push @b, split (/\n/, $a);
  }
  foreach $a (@b) {
    next if $a =~ /^#/;
    chomp $a;
    if ($a =~ s/^>//) {
      $i = 1;
      $name = $a;
      $test = $name;
      while (defined $seq->{$test}) {
        $test = $name . $i;
        $i++;
      }
      $name = $test;
      $seq->{$name} = "";
    } elsif (defined $name) {
      $seq->{$name} .= $a;
#      undef $name;
    }
  }
  return $seq;
}

sub hash2fasta {
  # returns an array
  # doesn't have newlines
  my ($seq) = @_;
  my $name;
  my @return;
  foreach $name (keys %$seq) {
    push @return, ">$name";
    push @return, $seq->{$name};
  }
  return @return;
}

sub tajimaD {
  # no p-values
  # expect sequence as aligned fasta set up in a hash
  my ($seq) = @_;

  my $debug = 0;
  my $n = scalar(keys %$seq);
  my $S = segregating_sites($seq);
  my $k = average_difference($seq);

  return 0 if !$S;
  my ($a1, $a2, $b1, $b2, $c1, $c2, $e1, $e2, $Dmin, $Dmax, $D);
  my $i;
  
  # a1 is sum of 1/i for i from 1 to n-1
  # a2 is sum of 1/i^2 for i from 1 to n-1
  $a1 = 0;
  $a2 = 0;
  foreach $i (1..$n-1) {
    $a1 += 1/$i;
    $a2 += 1/($i*$i);
  }
  # $debug && print "a1 = $a1\na2 = $a2\n";

  # b1 is (n+1)/3(n-1)
  $b1 = ($n + 1) / (3 * ($n-1));
  # b2 is 2(n^2 + n + 3)/9n(n-1)
  $b2 = 2 * ($n*$n + $n + 3) / (9 * $n * ($n-1));
  # $debug && print "b1 = $b1\nb2 = $b2\n";

  # c1 is b1 - 1/a1
  $c1 = $b1 - 1/$a1;
  # c2 is b2 - (n+2)/a1*n + a2/a1^2
  $c2 = $b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1);
  # $debug && print "c1 = $c1\nc2 = $c2\n";

  # e1 is c1/a1
  $e1 = $c1/$a1;
  # e2 is c2/(a1^2 + a2)
  $e2 = $c2/($a1*$a1 + $a2);
  # $debug && print "e1 = $e1\ne2 = $e2\n";

  # min is (2/n - 1/a1) / sqrt(e2)
  $Dmin = (2/$n - 1/$a1) / sqrt($e2);
  # max is ( n+1/2n - 1/a1 ) / sqrt(e2)
  $Dmax = (($n+1)/(2*$n) - 1/$a1) / sqrt($e2);
  # $debug && print "Dmin = $Dmin\nDmax = $Dmax\n";

  # now we can calculate D
  $D = ($k - $S/$a1) / sqrt($e1*$S + $e2*$S*($S-1));

  return ($D);
}

sub fuliDF {
  my ($seq) = @_;
  my $n;
  my $eta;
  my $eta_s;
  my $Dstar;
  my $pi_n;
  my $Fstar;

  $n = scalar (keys %$seq);
  ($eta, $eta_s) = eta($seq);
  return (0,0) if !$eta;
  $pi_n = average_difference($seq);

  $Dstar = $n / ($n-1) * $eta - a($n) * $eta_s;
  $Dstar /= sqrt(u_D($n) * $eta + v_D($n) * $eta * $eta);

  $Fstar = $pi_n - ($n - 1) / $n * $eta_s;
  $Fstar /= sqrt(u_F($n) * $eta + v_F($n) * $eta * $eta);

  return ($Dstar, $Fstar);

  sub a {
    my ($n) = @_;
    my $i;
    my $a = 0;
    foreach $i (1..$n-1) {
      $a += 1/$i;
    }
    return $a;
  }
  sub b {
    my ($n) = @_;
    my $i;
    my $b = 0;
    foreach $i (1..$n-1) {
      $b += 1/($i*$i);
    }
    return $b;
  }
  sub c {
    my ($n) = @_;
    my $c;
    if ($n == 2) {
      $c = 1;
    } elsif ($n > 2) {
      $c = 2 * ($n * a($n) - 2*($n-1)) / ($n-1) / ($n-2);
    }
    return $c;
  }
  sub d {
    my ($n) = @_;
    my $d = c($n);
    $d += ($n - 2) / ($n - 1) / ($n - 1);
    $d += 2 / ($n - 1) * (3/2 - (2*a($n+1) - 3)/($n-2) - 1/$n);
    return $d;
  }
  sub v_D {
    my ($n) = @_;
    my $v;
    $v = $n*$n/($n-1)/($n-1) * b($n);
    $v += a($n) * a($n) * d($n);
    $v -= 2 * $n * a($n) * (a($n)+1) / ($n - 1) / ($n - 1);
    $v /= a($n) * a($n) + b($n);
    return $v;
  }
  sub u_D {
    my ($n) = @_;
    my $u;
    $u = $n / ($n-1);
    $u *= a($n) - $n/($n-1);
    $u -= v_D($n);
    return $u;
  }
  sub v_F {
    my ($n) = @_;
    my $v;
    $v = 2*$n*$n*$n + 110*$n*$n - 255*$n + 153;
    $v /= 9 * $n * $n * ($n-1);
    $v += 2 * ($n-1) * a($n) / $n / $n;
    $v -= 8 * b($n) / $n;
    $v /= a($n) * a($n) + b($n);
    return $v;
  }
  sub u_F {
    my ($n) = @_;
    my $u;
    $u = 4*$n*$n + 19*$n + 3 - 12*($n+1)*a($n+1);
    $u /= 3 * $n * ($n - 1);
    $u /= a($n);
    $u -= v_F($n);
    return $u;
  }
  sub eta {
    my ($seq) = @_;
    my $name;
    my $i;
    my $data;
    my $eta = 0;
    my $eta_s = 0;
    foreach $name (keys %$seq) {
      foreach $i (0..length($seq->{$name})-1) {
        $data->{$i}->{uc substr($seq->{$name}, $i, 1)}++;
      }
    }
    foreach $i (keys %$data) {
      foreach $name (keys %{$data->{$i}}) {
        $eta_s++ if $data->{$i}->{$name} == 1;
        $eta++;
      }
      $eta-- if keys %{$data->{$i}};	# eta is number of nt at each position minus 1
    }
    return ($eta, $eta_s);
  }
}

sub faywuH {
  my ($seq) = @_;
  my $S = segregating_sites($seq);
  my $k = average_difference($seq);
  my $n = scalar keys %$seq;
  my $i;
  my $H;
  my $thetaH = 0;
  foreach $i (1..$n-1) {
    $thetaH += Si($seq, $i) * $i * $i;
  }
  $thetaH *= 2;
  $thetaH /= $n * ($n-1);
  $H = $k - $thetaH;
  return $H;
}

sub numhaplotypes {
  my ($seq) = @_;
  my $K;
  my $name;
  foreach $name (keys %$seq) {
    $K->{$seq->{$name}} = 1;
  }
  return scalar keys %$K;
}

sub homozygosity {
  # 1 - homozygosity is heterozygosity
  my ($seq) = @_;
  my $H = 0;
  my $hap;
  my $name;
  my $n = scalar keys %$seq;
  foreach $name (keys %$seq) {
    $hap->{$seq->{$name}} = 0 if !defined $hap->{$seq->{$name}};
    $hap->{$seq->{$name}}++;
  }
  foreach $name (keys %$hap) {
    $H += $hap->{$name} * $hap->{$name};
  }
  $H /= $n * $n;
  return $H;
}

sub hash_key_match {
  my ($a, $b) = @_;
  # expect two hash references
  # return the number of keys in common
  my $k;
  my $count = 0;
  foreach $k (keys %$a) {
    $count++ if defined $b->{$k};
  }
  return $count;
}

sub filter_redundant {
  # take a sequence hash
  # filter out redundant sequences
  # look for sequences that are identical or proper subsequences of another
  my ($seq) = @_;
  my @name = keys %$seq;
  my $filter;
  my $return;
  my $i;
  my $j;
  foreach $i (0..$#name-1) {
    next if $filter->{$name[$i]};
    foreach $j ($i+1..$#name) {
      next if $filter->{$name[$j]};
      if ($seq->{$name[$i]} eq $seq->{$name[$j]}) {
        $filter->{$name[$j]} = 1;
      } elsif ($seq->{$name[$i]} =~ /$seq->{$name[$j]}/) {
        $filter->{$name[$j]} = 1;
      } elsif ($seq->{$name[$j]} =~ /$seq->{$name[$i]}/) {
        $filter->{$name[$i]} = 1;
      }
    }
  }
  foreach $i (0..$#name) {
    next if $filter->{$name[$i]};
    $return->{$name[$i]} = $seq->{$name[$i]};
  }
  return $return;
}

my %orgcode_translate = ();
sub orgcode_translate {
  my @return = @_;
  if (!scalar keys(%orgcode_translate)) {
    eval "use Orgmap";
    my $orgmap = "$Orgmap::LIBPATH/org-map";
    $orgmap = "$Orgmap::LIBPATH/org-map";
    my @f;
    my @g;
    my @h;
    my $parsed_string;
    open ORGMAP, $orgmap;
    while (<ORGMAP>) {
      chomp;
      next if /^$/;
      next if /^$/;
      @f = split /\t/, $_;
      @g = split /\//, $f[1];
      if ($g[0] =~ /GENOMESPATH/) {
        @h = split /_/, $g[1];
      } else {
        @h = split /_/, $g[0];
      }
      $parsed_string = join ("_", @h[2..$#h-1]);
      $orgcode_translate{$f[0]} = $parsed_string;
      $orgcode_translate{$g[$#g]} = $parsed_string;
    }
    close ORGMAP;
  }
  my ($i, $code);
  foreach $i (0..$#return) {
    foreach $code (keys %orgcode_translate) {
      $return[$i] =~ s/$code/$orgcode_translate{$code}/g;
    }
  }
  return @return;
}

sub assembly_stats {
  # we will delete N's in the assembly and not count them - scaffold style
  my $fasta = shift;
  my $contig_cutoff;
  if (scalar @_) {
    $contig_cutoff = shift;
  } else {
    $contig_cutoff = 0;
  }

  my $key;
  my $seq;
  my @a;
  my $subtotal;
  my $total;
  my $return = ();
  my @length = ();

  $return->{number_contigs} = 0;
  $return->{total_length} = 0;
  $return->{avg_length} = 0;
  $return->{max_length} = 0;
  $return->{min_length} = 0;
  $return->{N50_length} = 0;
  $return->{N50_number} = 0;
  $return->{N90_length} = 0;
  $return->{N90_number} = 0;
  $return->{text} = "";

  if (-f $fasta) {
    open FASTA, $fasta;

    # length filter first
    @a = <FASTA>;
    close FASTA;
    $seq = fasta2hash(@a);
    $total = 0;
    foreach $key (keys %$seq) {
      $seq->{$key} =~ tr/nN//d;
      if (length($seq->{$key}) >= $contig_cutoff) {
        push @length, length($seq->{$key});
        $total += length($seq->{$key});
      }
    }
    @length = sort { $a <=> $b } @length;
    $subtotal = 0;
    foreach my $i (reverse 0..$#length) {
      $subtotal += $length[$i];
      if ($subtotal >= $total * 0.9 && $return->{N90_number} == 0) {
        $return->{N90_number} = $#length - $i + 1;
        $return->{N90_length} = $length[$i];
      }
      if ($subtotal >= $total/2 && $return->{N50_number} == 0) {
        $return->{N50_number} = $#length - $i + 1;
        $return->{N50_length} = $length[$i];
      }
    }
    $return->{number_contigs} = scalar(@length);
    $return->{total_length} = $total;
    $return->{avg_length} = int($total/scalar(@length)) if scalar(@length) > 0;
    $return->{max_length} = $length[$#length];
    $return->{min_length} = $length[0];
    foreach $key (qw(number_contigs total_length avg_length max_length min_length N50_length N50_number N90_length N90_number)) {
      $return->{text} .= "$key: $return->{$key}\n";
    }
  }

  return $return;
}

sub check_files {
  my ($q1, $q2) = @_;
  my ($abs1, $abs2);
  if (-f $q1) {
    $abs1 = File::Spec->rel2abs($q1);
  }
  if (-f $q2) {
    $abs2 = File::Spec->rel2abs($q2);
  }
  return ($abs1, $abs2);
}

sub info_from_mux {
  my ($mux) = @_;
  my $return;
  my $data;
  my $data2;
  my $content;
  my $content2;
  my $URL;
  my $URL2;
  my $json;
  my $basedir;
  my @dir;
  my $dir;
  my $machine;
  my $library;
  my $csvfile;
  my $csvdata;
  my $hash;
  my @f;
  my $start;
  my $end;
  my $plex;
  my $mux_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__MUX__/solexaRun/json";
  my $library_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__LIB__/expand/json";
  # need to fix this
  my $elm_client = "/mnt/software/lib/elmClient.cml.jar";

  if ($mux =~ /^\d+$/) {
    $mux = "MUX$mux";
  }
  if ($mux !~ /^MUX\d+$/) {
    die "Invalid MUX identifier $mux\n";
  }
  # get template first with insert length
  $json = new JSON;
  $URL = $library_URL;
  $URL =~ s/__LIB__/$mux/;
  $data = `java -jar $elm_client $URL`;
  $content = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data);
  $start = $content->{templates}->{BP_Range_Start};
  $end = $content->{templates}->{BP_Range_End};
  $return->{INSERTLENGTH} = ($start + $end) / 2;

  # now plexes and other run info
  $URL = $mux_URL;
  $URL =~ s/__MUX__/$mux/;
  $data = `java -jar $elm_client $URL`;
  $content = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data);

  $return->{RUNID} = $content->{runs}->{runId};
  $return->{LANEID} = $content->{runs}->{lanes}->{laneId};
  @f = split /-/, $content->{runs}->{runId};
  $machine = $f[0];
  @dir = glob("/mnt/SolexaPool/$machine/".$return->{RUNID}."*");
  foreach $dir (@dir) {
    $basedir = "$dir/Data/Intensities/BaseCalls/Demultiplexed";
    $csvfile = "$basedir/SamplesDirectories.csv";
    last if -f $csvfile;
  }
  if (-f $csvfile) {
    open F, $csvfile;
    while (<F>) {
      chomp;
      @f = split /,/, $_;
      # key on lane then barcode
      $csvdata->{$f[1]}->{$f[4]}->{SampleID} = $f[2];
      $csvdata->{$f[1]}->{$f[4]}->{Directory} = $f[9];
    }
    close F;
  } else {
    print "can't find $csvfile\n";
    return (0);
  }

  foreach $plex (0..$#{$content->{plexes}}) {
    $hash = $content->{plexes}->[$plex];
    $library = $hash->{libraryId};
    $URL2 = $library_URL;
    $URL2 =~ s/__LIB__/$library/;
    $data2 = `java -jar $elm_client $URL2`;
    $content2 = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data2);
    $return->{$plex}->{LIBRARY} = $hash->{libraryId};
    $return->{$plex}->{BARCODE} = $hash->{barCode};
    $return->{$plex}->{DESCRIPTION} = $content2->{description};
    $return->{$plex}->{INDEXDIR} = $csvdata->{$return->{LANEID}}->{$hash->{barCode}}->{Directory};
    $return->{$plex}->{FULLDIR} = "$basedir/$return->{$plex}->{INDEXDIR}";
  }

  return $return;
}

sub get_data_bydir {
  my ($dir, $lane, $final_sample_name) = @_;
  my @q1files = ();
  my @q2files = ();
  my $q1;
  my $q2;
  my $f;
  my $expected_file;
  my $command;
  my $output;
  my $tempdir = set_option("tempdir");
  my $tempfile = "$tempdir/temp_qseq.txt";
  # need to fix this
  my $programs;
  $programs->{qseq2fastq} = "/mnt/software/bin/illumina2fastq";

  foreach $f (glob("$dir/s_".$lane."_1_*_qseq.txt")) {
    push @q1files, $f;
  }
  foreach $f (glob("$dir/s_".$lane."_2_*_qseq.txt")) {
    push @q2files, $f;
  }
  if ($#q1files >= 0 && $#q1files == $#q2files) {
    # first $q1
    $expected_file = "$tempdir/$final_sample_name.1_fastq.txt";
    $command = "cat " . join (" ", @q1files) . " > $tempfile && $programs->{qseq2fastq} -i $tempfile -o $expected_file 2>&1 && rm -f $tempfile";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (!-f $expected_file) {
      &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
      return (0,0);
    }
    $q1 = $expected_file;

    # now $q2
    $expected_file = "$tempdir/$final_sample_name.2_fastq.txt";
    $command = "cat " . join (" ", @q2files) . " > $tempfile && $programs->{qseq2fastq} -i $tempfile -o $expected_file 2>&1 && rm -f $tempfile";
    &shortlog($command);
    $output = `$command 2>&1`;
    &log($output);
    if (!-f $expected_file) {
      &log("Couldn't find $expected_file after ".(caller(0))[3]." run");
      return (0,0);
    }
    $q2 = $expected_file;
  }
  if (defined $q1 && -f $q1 && -f $q2) {
    return ($q1, $q2);
  } else {
    return (0,0);
  }
}


sub info_from_mux_182 {
  my ($mux) = @_;
  my $return;
  my $data;
  my $data2;
  my $content;
  my $content2;
  my $URL;
  my $URL2;
  my $json;
  my $basedir;
  my $sdir;
  my $sampledir;
  my @dir;
  my $dir;
  my $machine;
  my $library;
  my $csvfile;
  my @csvfile;
  my $csvdata;
  my $hash;
  my @f;
  my $start;
  my $end;
  my $plex;
  my $i;
  my %template_nums;
  my $max_template_num;
  my $mux_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__MUX__/solexaRun/json";
  my $library_URL = "http://elm.gis.a-star.edu.sg/rest/libinfo/__LIB__/expand/json";
  # need to fix this
  my $programs;
  my $elm_client = "/mnt/software/lib/elmClient.cml.jar";

  if ($mux =~ /^\d+$/) {
    $mux = "MUX$mux";
  }
  if ($mux !~ /^MUX\d+$/) {
    die "Invalid MUX identifier $mux\n";
  }
  # get template first with insert length
  $json = new JSON;
  $URL = $library_URL;
  $URL =~ s/__LIB__/$mux/;
  $data = `java -jar $elm_client $URL`;
  $content = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data);
  if (ref($content->{templates}) eq "HASH") {
    $start = $content->{templates}->{BP_Range_Start};
    $end = $content->{templates}->{BP_Range_End};
  } else {
    # should be an array ref and should have multiple templates
    undef %template_nums;
    $max_template_num = 0;
    foreach $i (0..$#{$content->{templates}}) {
      $content->{templates}->[$i]->{template_name} =~ /T(\d+)$/;
      if ($1 > $max_template_num) {
        $max_template_num = $1;
      }
      $template_nums{$1} = $i;
    }
    if (defined $template_nums{$max_template_num} && defined $content->{templates}->[$template_nums{$max_template_num}]->{BP_Range_Start}) {
      $start = $content->{templates}->[$template_nums{$max_template_num}]->{BP_Range_Start};
      $end = $content->{templates}->[$template_nums{$max_template_num}]->{BP_Range_End};
    }
  }
  $return->{INSERTLENGTH} = ($start + $end) / 2;

  # now plexes and other run info
  $URL = $mux_URL;
  $URL =~ s/__MUX__/$mux/;
  $data = `java -jar $elm_client $URL`;
  $content = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data);

  if (ref($content->{runs}) eq "HASH") {
    $return->{RUNID} = $content->{runs}->{runId};
    $return->{LANEID} = $content->{runs}->{lanes}->{laneId};
    @f = split /-/, $content->{runs}->{runId};
  } else {
    # should be an array ref and then have multiple runs
    # for now let's just take the last one
    my $last_index = $#{$content->{runs}};
    $return->{RUNID} = $content->{runs}->[$last_index]->{runId};
    $return->{LANEID} = $content->{runs}->[$last_index]->{lanes}->{laneId};
    @f = split /-/, $content->{runs}->[$last_index]->{runId};
  }
  $machine = $f[0];

  my $all_dir_path = "/mnt/SolexaPool/$machine/".$return->{RUNID}."*/output_".$return->{RUNID}."*/CASAVA*";
   @dir = glob($all_dir_path);

  if (scalar(@dir) == 0){
      print "Missing files in the directory $all_dir_path\n";
      print "Exiting....\n";
      exit(1);
  }
  foreach $dir (@dir) {
    next if $dir =~ /toDelete/;
    $basedir = "$dir/Project_$mux";
    next if !-d $basedir;
    # as of 28 Oct 2013 - don't seem to have these .csv files any more
    # we just need this for the directory - essentially the csv file tells us
    # there is data in the directory - so have to change this
    @csvfile = glob("$basedir/*/SampleSheet.csv");
    foreach $csvfile (@csvfile) {
      if (-f $csvfile) {
        open F, $csvfile;
        while (<F>) {
          chomp;
          @f = split /,/, $_;
          # key on lane then barcode
          $csvdata->{$f[1]}->{$f[4]}->{SampleID} = $f[2];
          $csvdata->{$f[1]}->{$f[4]}->{Directory} = $basedir;
          # make another with lane then library id
          # seem to not always be able to trust info in ELM - this prioritizes
          # the info we find in the csv files
          $csvdata->{$f[1]}->{$f[2]}->{barCode} = $f[4];
          $csvdata->{$f[1]}->{$f[2]}->{Directory} = $basedir;
        }
        close F;
      }
    }
    opendir(PROJECTDIR, $basedir);
    while ($sdir = readdir PROJECTDIR) {
      if ($sdir =~ /^Sample_(.*)$/) {
        $sampledir->{$1} = "$basedir/$sdir";
      }
    }
    closedir(PROJECTDIR);
  }

  foreach $plex (0..$#{$content->{plexes}}) {
    $hash = $content->{plexes}->[$plex];
    $library = $hash->{libraryId};
    $URL2 = $library_URL;
    $URL2 =~ s/__LIB__/$library/;
    $data2 = `java -jar $elm_client $URL2`;
    $content2 = $json->allow_nonref->utf8->relaxed->escape_slash->loose->allow_singlequote->allow_barekey->decode($data2);
    $return->{$plex}->{LIBRARY} = $hash->{libraryId};
    $return->{$plex}->{DESCRIPTION} = $content2->{description};
    if (defined $csvdata->{$return->{LANEID}}->{$hash->{libraryId}}) {
      # default now to use info in SampleSheet.csv files
      $return->{$plex}->{BARCODE} = $csvdata->{$return->{LANEID}}->{$hash->{libraryId}}->{barCode};
      $return->{$plex}->{INDEXDIR} = $csvdata->{$return->{LANEID}}->{$hash->{libraryId}}->{Directory}  ."/Sample_$hash->{libraryId}";
    } elsif (defined $csvdata->{$return->{LANEID}}->{$hash->{barCode}}) {
      # fallback is to try to use info in ELM
      $return->{$plex}->{BARCODE} = $hash->{barCode};
      $return->{$plex}->{INDEXDIR} = $csvdata->{$return->{LANEID}}->{$hash->{barCode}}->{Directory}  ."/Sample_$hash->{libraryId}";
    } elsif (defined $sampledir->{$hash->{libraryId}}) {
      # we get here if we have no CSV file, which seems to be in newer runs
      # as of 28 Oct 2013
      $return->{$plex}->{BARCODE} = $hash->{barCode};
      $return->{$plex}->{INDEXDIR} = $sampledir->{$hash->{libraryId}};
    } else {
      # we're in trouble we can't figure out the directories
      print STDERR "Can't find data directory for $hash->{libraryId}, barcode $hash->{barCode} in $all_dir_path\n";
      return(0);
    }
    $return->{$plex}->{FULLDIR} = $return->{$plex}->{INDEXDIR};
  }

  return $return;
}

sub get_data_bydir_182 {
  my ($index_dir, $tempdir, $final_sample_name, $gz) = @_;
  my ($q1, $q2);
  my $f;
  my $command;
  my $output;
  my $expected_file;
  my $read;
  my @allfiles;
  my @sortedfiles;
  my @expected_files;
  my $rtype;
  my $line;
  my $file_cnt;
  my $tempfile;
  my $filename;
  my $read_cnt_before_filtering;
  my $read_cnt_after_filtering;

  if (-d $index_dir) {
    # first check in $index_dir, we're looking for gzipped fastq files now
    opendir(CASAVA182, $index_dir) or die "Couldn't open the index directory\n";

    @allfiles = readdir(CASAVA182);
    foreach $f (@allfiles){
      next if ($f =~ m/^\./);
      if ($f =~ m/R([12])_\d+.fastq/) {
        push(@sortedfiles, $f);
      }
    }
    $file_cnt = scalar(@sortedfiles);
    if ($file_cnt == 0){
      print "No files found in the directory: $index_dir\n";
      print "Exiting....\n";
      exit(1);
    }

    if (defined $gz && $gz eq "nogz") {
      $expected_files[0] = "$tempdir/$final_sample_name" . "_R1" . ".fastq";
      $expected_files[1] = "$tempdir/$final_sample_name" . "_R2" . ".fastq";
    } else {
      $expected_files[0] = "$tempdir/$final_sample_name" . "_R1" . ".fastq.gz";
      $expected_files[1] = "$tempdir/$final_sample_name" . "_R2" . ".fastq.gz";
    }

    foreach $rtype (1, 2) {
      # Keep only the reads with good quality values (i.e. with N flag),
      # bad ones have Y flag.
      # grep -A3 is used to force separate the 3 lines associated with each
      # read; this will add a line starting '--'. Need to remove this line!

      $filename = "*_R" .$rtype. "_*.gz";
      $tempfile = "$tempdir/$final_sample_name" . "_temp" . ".txt";

      # report filtered out reads (i.e. the ones with Y flag) =  number of reads before filtering out 'Y' flags - after filtering 'Y' flag 
      if (defined $gz && $gz eq "nogz") {
        system "zcat $index_dir/$filename | tee '$tempfile' | grep -A3 '^@.* [12]:N:' | grep -v -e '^--\$' > '$expected_files[$rtype-1]'";
      } else {
        system "zcat $index_dir/$filename | tee '$tempfile' | grep -A3 '^@.* [12]:N:' | grep -v -e '^--\$' | gzip > '$expected_files[$rtype-1]'";
      }
      print "End of creating $expected_files[$rtype-1]\n";
      $read_cnt_before_filtering = get_read_cnt($tempfile);
      $read_cnt_after_filtering = get_read_cnt($expected_files[$rtype-1]);
      print "Number of reads filtered out = " . ($read_cnt_before_filtering - $read_cnt_after_filtering) . "\n";
    }
    unlink $tempfile;
    if ((!-f $expected_files[0])  || (!-f $expected_files[1])){
      return (0,0);
    }
    $q1 = $expected_files[0];
    $q2 = $expected_files[1];
  }
  if (defined $q1 && -f $q1 && -f $q2) {
    return ($q1, $q2);
  } else {
    return (0,0);
  }
}

sub get_read_cnt {
  my ($filename) = @_;
  if ($filename =~ /\.gz$/) {
    open(FILE, "gunzip -c '$filename' |")  || die ("Couldn't open pipe to file $filename\n") ;
  } else {
    open(FILE, $filename)  || die ("Couldn't open the file $filename\n") ;
  }

  my @lines = <FILE>;
  my $line_cnt = scalar(@lines);	
  my $read_cnt = $line_cnt/4;
  close(FILE);
  return($read_cnt);
}

sub annotate_by_ref {
  my ($faa, $gbk, $ref) = @_;
  my $tempdir = File::Temp::tempdir(CLEANUP => 1);
  my $currentdir = Cwd::getcwd();
  my @annotated;
  my $extra;
  my @f;
  my @g;
  my $key;
  my $hit;
  my $qlen;
  my $slen;
  my $align;
  my $ident;
  my $string;
  my $gi;
  my $basename;
  my $abs_faa;
  my $ref_ptt = $ref;
  my $ref_faa = $ref;
  $ref_ptt =~ s/\.fna$/.ptt/;
  $ref_faa =~ s/\.fna$/.faa/;
  my $ptt;
  open P, $ref_ptt;
  while (<P>) {
    chomp;
    next if !/^\s*\d+\.\.\d+/;
    @f = split /\t/, $_;
    $ptt->{$f[3]}->{TRIVIAL} = $f[4];
    $ptt->{$f[3]}->{SYSTEMATIC} = $f[5];
    $ptt->{$f[3]}->{ANNOTATION} = $f[8];
  }
  close P;

  $abs_faa = File::Spec->rel2abs($faa);
  chdir($tempdir);

  # blast
  $basename = File::Basename::basename($ref_faa);
  File::Copy::copy($ref_faa, $tempdir);
  system "makeblastdb -in $basename -dbtype prot 2>&1 > /dev/null";

  # prodigal gives us this kind of fasta header:
  # >GERMS-WGS.pl_5 # 2651 # 4111 # -1 # ID=1_5;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.327 
  # we want the ID=1_5 part to show up in the blast results because the
  # the .gbk file will have that in there also - so mangle the headers first
  system "sed -e '/^>/s/^>.*ID=/>/' -e '/^>/s/;/ /' $abs_faa | blastp -query - -db $basename -outfmt '6 qseqid sseqid qlen slen length pident' -max_target_seqs 1 > $tempdir/blast_temp_file";

  # parse blast
  $extra = ();
  open B, "$tempdir/blast_temp_file";
  while (<B>) {
    chomp;
    @f = split /\t/, $_;
    $key = $f[0];
    $hit = $f[1];
    @g = split /\|/, $hit;
    $gi = $g[1];
    $qlen = $f[2];
    $slen = $f[3];
    $align = $f[4];
    $ident = $f[5];
    $string = "";
    if (defined $ptt->{$gi}) {
      if ($ptt->{$gi}->{TRIVIAL} eq "-") {
        $string .= "                     /gene=\"$ptt->{$gi}->{SYSTEMATIC}\"\n";
      } else {
        $string .= "                     /gene=\"$ptt->{$gi}->{TRIVIAL}\"\n";
      }
      $string .= "                     /locus_tag=\"$ptt->{$gi}->{SYSTEMATIC}\"\n";
      $string .= "                     /product=\"$ptt->{$gi}->{ANNOTATION}; ";
      $string .= "$ident\% ID over $align aa, ";
      $string .= sprintf("%.2f of Query, %.2f of Reference\"\n", $align/$qlen, $align/$slen);
    } else {
    }
    $extra->{$key} = $string;
  }
  close B;
  unlink("$tempdir/blast_temp_file");
  unlink("$tempdir/$basename");
  unlink("$tempdir/$basename.pin");
  unlink("$tempdir/$basename.phr");
  unlink("$tempdir/$basename.psq");

  chdir($currentdir);
  
  # fix gbk
  @annotated = ();
  open G, $gbk;
  while (<G>) {
    if (/\/note="ID=(.*?);/) {
      $key = $1;
      if (defined($extra->{$key})) {
        push @annotated, $extra->{$key};
      }
    }
    push @annotated, $_;
  }
  close G;
  unlink($gbk);

  open G, ">$gbk";
  print G @annotated;
  close G;
  return ($gbk);
}

sub get_browser_tip {
  my ($Run, $DBH, $USE_DB) = @_;
  return 0 if !length($Run);
  my @data;
  my $check_sql = "SELECT TIP FROM Tips WHERE Run = ?";
  my $insert_sql = "INSERT INTO Tips (Run) VALUES (?)";
  my $sth;
  # first check if this exists
  $sth = $DBH->prepare($check_sql);
  $sth->execute($Run);
  while (@data = $sth->fetchrow_array()) {
    if ($data[0] > 0) {
      return($data[0]);
    }
  }
  # otherwise we're making a new one
  if ($USE_DB) {
    $sth = $DBH->prepare($insert_sql);
    $sth->execute($Run);
    $sth = $DBH->prepare($check_sql);
    $sth->execute($Run);
    while (@data = $sth->fetchrow_array()) {
      if ($data[0] > 0) {
        return($data[0]);
      }
    }
  } else {
    return(0);
  }
  # if we still don't have a tip number, then we have a problem
  die "Can't seem to get a TIP number for $Run\n";
}

sub get_DateStamp {
  my ($id, $type, $md5, $DBH, $USE_DB) = @_;
  my $tip = $id;
  my @data;
  my ($sql, $sth);
  my $DB_PARSER = DateTime::Format::DBI->new($DBH);
  if ($tip !~ /^\d+$/) {
    $tip = get_browser_tip($tip, $DBH, $USE_DB);
  } 
  if ($tip == 0) {
    die ("Can't find good Tip number for $id, $type, $md5\n");
  }
  $sql = "SELECT DateStamp FROM Files WHERE TIP = ? AND Type = ? AND MD5 = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($tip, $type, $md5);
  while (@data = $sth->fetchrow_array()) {
    if (length($data[0])) {
      return(DateTime::Format::DateParse->parse_datetime($data[0]));
    }
  } 
  return "NULL";
}

sub get_DateStampStudies {
  # this is different from the postprocess one
  my ($run, $DBH) = @_;
  my @data;
  my ($sql, $sth);
  my $DB_PARSER = DateTime::Format::DBI->new($DBH);
  $sql = "SELECT DateStamp FROM Studies WHERE Run = ?";
  $sth = $DBH->prepare($sql);
  $sth->execute($run);
  while (@data = $sth->fetchrow_array()) {
    if (length($data[0])) {
      return(DateTime::Format::DateParse->parse_datetime($data[0]));
    }
  }
  return "NULL";
}

sub get_columns {
  my ($table, $DBH) = @_;
  my $sql = "DESCRIBE $table";
  my @data;
  my $return;
  my $sth = $DBH->prepare($sql);
  $sth->execute();
  while (@data = $sth->fetchrow_array()) {
    $return->{$data[0]}->{TYPE} = $data[1];
    $return->{$data[0]}->{NULL} = $data[2];
  }
  return $return; 
}

sub get_indexes {
  my ($table, $DBH) = @_;
  my $sql = "SHOW INDEXES FROM $table";
  my @data;
  my $return;
  my $sth = $DBH->prepare($sql);
  $sth->execute();
  while (@data = $sth->fetchrow_array()) {
    $return->{$data[4]} = 1;
  }
  return $return;
}

sub do_db_withSourceFile {
  my ($data, $table, $force, $DBH, $USE_DB) = @_;
  my (@data, $i, $j, $col);
  my $r;
  my $mode = "INSERT";
  my $expected = get_columns($table, $DBH);
  my @exp = sort keys %$expected;
  my $keylist = get_indexes($table, $DBH);
  my @keyl = sort keys %$keylist;
  my @update_exp = @exp;
  my $db_date;
  my $nodot;
  my $newdate;
  my $existingdate;
  my $sth;
  my $key_condition;
  my $check_sql;
  my $insert_sql;
  my $update_sql;

  # for filenames substitute "." to "_"
  foreach $i (keys %$data) {
    $nodot = $i;
    $nodot =~ s/\./_/g;
    if (!defined $data->{$nodot}) {
      $data->{$nodot} = $data->{$i};
    }
  }
  foreach $i (reverse 0..$#update_exp) {
    foreach $j (@keyl) {
      if ($update_exp[$i] eq $j) {
        splice @update_exp, $i, 1;
        last;
      }
    }
  }
  $key_condition = join (" = ? AND ", @keyl) . " = ?";
  $check_sql = "SELECT * FROM $table WHERE $key_condition";
  $insert_sql = "INSERT INTO $table (" . join (", ", @exp) . ") VALUES ";
  $insert_sql .= "(" . join (", ", ("?") x scalar(@exp)) . ")";
  $update_sql = "UPDATE $table SET ";
  $update_sql .= join (" = ?, ", @update_exp) . " = ?";
  $update_sql .= " WHERE $key_condition";

  # check we have all the data we need
  foreach $col (@exp) {
    if (!defined $data->{$col}) {
      if ($expected->{$col}->{NULL} eq "NO") {
        die "Have null for non-null column $col in $table\n" . Dumper(\$data);
      }
    }
  }

  # check if we have something there already
  $sth = $DBH->prepare($check_sql);
  foreach $i (0..$#keyl) {
    $sth->bind_param($i+1, $data->{$keyl[$i]});
  }
  $sth->execute;
  $r = $sth->fetchrow_hashref();
  if (defined $r) {
    if (defined $expected->{SourceFileType}) {
      $newdate = get_DateStamp($data->{TIP}, $data->{SourceFileType}, $data->{SourceFileMD5}, $DBH, $USE_DB);
      $existingdate = get_DateStamp($r->{TIP}, $r->{SourceFileType}, $r->{SourceFileMD5}, $DBH, $USE_DB);
      if ($existingdate eq "NULL" || ($newdate > $existingdate) ||
          uc($force) eq "ALL" ||
          uc($force) eq uc($table)) {
        $mode = "UPDATE";
        } else {
        return "NOTHING";
      }
    } else {
      if ($force) {
        $mode = "UPDATE";
      } else {
        return "NOTHING";
      }
    }
  } else {
    $mode = "INSERT";
  }

  # do the actual imports
  if ($USE_DB) {
    if ($mode eq "INSERT") {
      $sth = $DBH->prepare($insert_sql);
      foreach $i (0..$#exp) {
        $sth->bind_param($i+1, $data->{$exp[$i]});
      }
      if (defined $sth->execute()) {
        return "INSERT";
      } else {
        return "ERROR";
      }
    } elsif ($mode eq "UPDATE") {
      $sth = $DBH->prepare($update_sql);
      foreach $i (0..$#update_exp) {
        $sth->bind_param($i+1, $data->{$update_exp[$i]});
      }
      foreach $i (0..$#keyl) {
        $sth->bind_param($i+1+scalar(@update_exp), $data->{$keyl[$i]});
      }
      if (defined $sth->execute()) {
        return "UPDATE";
      } else {
        return "ERROR";
      }
    }
  } else {
    return $mode;
  }
}

sub do_db {
  my ($data, $table, $force, $DBH, $USE_DB) = @_;
  my (@data, $i, $j, $col);
  my $r;
  my $mode = "INSERT";
  my $expected = get_columns($table, $DBH);
  my @exp = sort keys %$expected;
  my $keylist = get_indexes($table, $DBH);
  my @keyl = sort keys %$keylist;
  my @update_exp = @exp;
  my $db_date;
  my $nodot;
  my $newdate;
  my $existingdate;
  my $sth;
  my $key_condition;
  my $check_sql;
  my $insert_sql;
  my $update_sql;
  my $change = 0;

  # for filenames substitute "." to "_"
  foreach $i (keys %$data) {
    $nodot = $i;
    $nodot =~ s/\./_/g;
    if (!defined $data->{$nodot}) {
      $data->{$nodot} = $data->{$i};
    }
  }
  foreach $i (reverse 0..$#update_exp) {
    foreach $j (@keyl) {
      if ($update_exp[$i] eq $j) {
        splice @update_exp, $i, 1;
        last;
      }
    }
  }
  $key_condition = join (" = ? AND ", @keyl) . " = ?";
  $check_sql = "SELECT * FROM $table WHERE $key_condition";
  $insert_sql = "INSERT INTO $table (" . join (", ", @exp) . ") VALUES ";
  $insert_sql .= "(" . join (", ", ("?") x scalar(@exp)) . ")";
  $update_sql = "UPDATE $table SET ";
  $update_sql .= join (" = ?, ", @update_exp) . " = ?";
  $update_sql .= " WHERE $key_condition";

  # check we have all the data we need
  foreach $col (@exp) {
    if (!defined $data->{$col}) {
      if ($expected->{$col}->{NULL} eq "NO") {
        die "Have null for non-null column $col in $table\n" . Dumper(\$data);
      }
    }
  }

  # check if we have something there already
  $sth = $DBH->prepare($check_sql);
  foreach $i (0..$#keyl) {
    $sth->bind_param($i+1, $data->{$keyl[$i]});
  }
  $sth->execute;
  $r = $sth->fetchrow_hashref();
  if (defined $r) {
    foreach $i (keys %$r) {
      next if $i eq "DateStamp";
      if (defined $r->{$i} || defined $data->{$i}) {
        if (!defined $r->{$i} || !defined $data->{$i}) {
          $change++;
        } else {
          $change++ if $r->{$i} ne $data->{$i};
        }
      }
    }
    if (defined $data->{DateStamp}) {
      $newdate = $data->{DateStamp};
      $existingdate = get_DateStampStudies($data->{Run}, $DBH);
      if ( ( $existingdate eq "NULL" || ($newdate > $existingdate && $change) )
           || $force) {
        $mode = "UPDATE";
      } else {
        return "NOTHING";
      }
    } elsif ($change || $force) {
      $mode = "UPDATE";
    } else {
      return "NOTHING";
    }
  } else {
    $mode = "INSERT";
  }

  # do the actual imports
  if ($USE_DB) {
    if ($mode eq "INSERT") {
      $sth = $DBH->prepare($insert_sql);
      foreach $i (0..$#exp) {
        $sth->bind_param($i+1, $data->{$exp[$i]});
      }
      if (defined $sth->execute()) {
        return "INSERT";
      } else {
        return "ERROR";
      }
    } elsif ($mode eq "UPDATE") {
      $sth = $DBH->prepare($update_sql);
      foreach $i (0..$#update_exp) {
        $sth->bind_param($i+1, $data->{$update_exp[$i]});
      }
      foreach $i (0..$#keyl) {
        $sth->bind_param($i+1+scalar(@update_exp), $data->{$keyl[$i]});
      }
      if (defined $sth->execute()) {
        return "UPDATE";
      } else {
        return "ERROR";
      }
    }
  } else {
    return $mode;
  }
}

sub downsample {
  # $file should be absolute path
  my ($file, $SEQTK, $num_reads, $seed, $tempdir) = @_;
  my $reads;
  my $command;
  my $output;
  my $outfile;
  # /mnt/software libs are messing up local workstation awk
  my $ENV_save = $ENV{LD_LIBRARY_PATH};
  my @ft = GERMS::file_type($file);
  if (!defined $tempdir) {
    $tempdir = ".";
  }
  $outfile = "$tempdir/$num_reads-" . File::Basename::basename($file);
  $ENV{LD_LIBRARY_PATH} = "";
  if (lc($ft[1]) eq "gzip") {
    $reads = `zcat $file | wc -l | awk '{printf "\%.0f", \$1/4}'`;
  } else {
    $reads = `wc -l $file | awk '{printf "\%.0f", \$1/4}'`;
  }
  $ENV{LD_LIBRARY_PATH} = $ENV_save;
  if ($num_reads < $reads) {
    if (lc($ft[1]) eq "gzip") {
      $command = "$SEQTK sample -s $seed $file $num_reads | gzip > $outfile";
    } else {
      $command = "$SEQTK sample -s $seed $file $num_reads > $outfile";
    }
    $output = `$command 2>&1`;
  } else {
    return($file);
  }
  if (-f $outfile && -s $outfile) {
    return($outfile);
  } else {
    return("");
  }
}

sub file_type {
  my ($file) = @_;
  # figure out fasta/fastq, figure out if gzipped
  # return will be 2 element array: fasta/fastq and gzip/text
  # may also get none for either return value
  my $ft = File::Type->new();
  my $fh;
  my @head;
  my $string;
  my @return = ("none", "none");
  if (-f $file && -s $file) {
    $file = File::Spec->rel2abs($file);
  } else {
    return(@return);
  }
  $string = $ft->checktype_filename($file);
  if ($string eq "application/x-gzip" || $string =~ /gzip$/) {
    $return[1] = "gzip";
  } elsif ($string eq "application/octet-stream") {
    $return[1] = "text";
  } elsif ($file =~ /\.gz$/) {
    $return[1] = "gzip";
  } else {
    $return[1] = "text";
  }

  if ($return[1] eq "gzip") {
    $fh = Compress::Zlib::gzopen($file, "rb");
    foreach $string (0..3) {
      $fh->gzreadline($head[$string]);
    }
    $fh->gzclose;
  } else {
    @head = `head $file -n 4`;
  }
  if ($head[0] =~ /^>/) {
    $return[0] = "fasta";
  } elsif ($head[0] =~ /^@/) {
    if ($head[2] =~ /^\+/) {
      $return[0] = "fastq";
    }
  }
  return(@return);
}

# return 1 to make Perl happy
1;
