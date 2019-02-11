#!/bin/bash
INPUT=$1	# this basically needs to be like GCA_000157635.1
		# or GCA_000157635.1_ASM15763v1_genomic.fna
		# without a path
if [ -z $INPUT ]; then
  echo "Usage: $0 <assembly ID>"
  echo "Will download if fna file can't be found in the species directory"
  echo "assembly ID is of the form GCA_000157635.1 or GCA_000157635.1_ASM15763v1_genomic.fna"
  exit
fi

FIELD=19	# GCF
FIELD=20	# GCA

# get real script directory - from https://stackoverflow.com/questions/59895/get-the-source-directory-of-a-bash-script-from-within-the-script-itself
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

# set some common variables
if [ -f $DIR/../etc/browser-conf.sh ]; then
  . $DIR/../etc/browser-conf.sh
fi

function process_fna {
  nucmer --maxmatch --coords $REF $FNA -p $RUN
  show-snps -CT $RUN.delta > $RUN.snps
  coords2gcov.pl $RUN.coords | gzip > $RUN.gcov.gz
  snp2vcf.pl $RUN.snps | bgzip -c > $RUN.lofreq.gz
  tabix -p vcf $RUN.lofreq.gz
  rm -f $RUN.delta $RUN.snps $RUN.coords
  blast-mlst.pl -species $SPECIES -name $RUN $FNA | grep -v '^# Sample' | sed -e "s/^/# /" > $RUN.srst2
  blast-vf.pl -species $SPECIES -name $RUN $FNA -report_all_alleles >> $RUN.srst2
  gzip $RUN.srst2
  gcov-graph.R $RUN.gcov.gz
}

# figure some stuff out
# could be GCA_000157635.1, GCA_000157635.1_ASM15763v1_genomic.fna
INBASE=$( basename $INPUT )
RUN=$( echo $INBASE | cut -f1,2 -d'_' )
# genomes_proks.txt has MS-DOS line feeds at the end of each line
URL=$( grep $RUN $GPROKS | cut -f$FIELD | perl -ne 'chomp; s/\r$//; print;' )
FNA=$( basename ${URL#ftp:} )
FNA=${FNA}_genomic.fna
URL=$URL/$FNA.gz
SPECIES=$( grep $RUN $GPROKS | cut -f1,2 -d' ' )
RUN=$( echo $FNA | cut -f1,2 -d'_' )

# first get some info
# $MYSQL_REF provides data as environment variables
# WORK is the working directory
# REF is the reference genome
# SPECIES is the short species name
# FULLSPECIES is the long species name
eval $( $MYSQL_REF $RUN -species "$SPECIES" )

if [ ! -d $WORK ]; then
  sleep 5
  echo "Current WORK ($WORK) doesn't exist."
  echo "Sleeping 5 to wait for database."
  eval $( $MYSQL_REF $RUN -species "$SPECIES" )
fi

if [ ! -d $WORK ]; then
  exit
fi

CURRENT=$PWD
cd $WORK

if [ "$REF" = "NONE" ]; then
  echo "Can't find reference for $1"
  exit
fi

if [ ! -f $FNA ]; then
  wget $URL
  gunzip $FNA.gz
  process_fna
  rm -f $FNA
else
  # the INPUT should really be the fna file
  process_fna
fi

postprocess-db.pl -db -v $RUN -source ASSEMBLY -mlst $SPECIES -s3 -delete
cd $CURRENT
