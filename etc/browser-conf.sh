# This is file browser-conf.sh - should be sourced to set some variables

# targets depend on whether we have a reference or not
export NOREF_TARGETS="noref_final"
export TARGETS="final wgs"

# binaries and files
export SNAKE=<path_to_snakemake>/snakemake
# qsub defaults to a csh for submission which has different redirection syntax
# set this to use bash
QSUB="/opt/uge-8.1.7p3/bin/lx-amd64/qsub -S /bin/bash"
QSTAT=/opt/uge-8.1.7p3/bin/lx-amd64/qstat
QALTER=/opt/uge-8.1.7p3/bin/lx-amd64/qalter
export GERMS_DATA=<path_to_sequence_files>
export STAGING=<path_to_staging_dir>
export SRADB=<path_to_SRAdb>

# usually BROWSER_BASE probably will be something like /usr/local
export BROWSER_BASE=<path_to_install_location>
export FQSNAKE=$BROWSER_BASE/lib/browser-fastq.snakefile
export SNPSNAKE=$BROWSER_BASE/lib/browser-snpcall.snakefile
export CONF=$BROWSER_BASE/etc/browser-snakemake.json
export MYSQL_REF=$BROWSER_BASE/bin/mysql-get-reference.pl

# scale threads to the machine
export THREADS=$(grep -c '^processor' /proc/cpuinfo)
# but don't go over some maximum (4 here)
if [ $THREADS -gt 4 ]; then
  export THREADS=4
fi
