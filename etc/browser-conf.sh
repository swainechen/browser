# This is file browser-conf.sh - should be sourced to set some variables

# targets depend on whether we have a reference or not
NOREF_TARGETS="noref_final"
TARGETS="final wgs"

# binaries and files
SNAKE=<path_to_snakemake>/snakemake
# qsub defaults to a csh for submission which has different redirection syntax
# set this to use bash
QSUB="/opt/uge-8.1.7p3/bin/lx-amd64/qsub -S /bin/bash"
QSTAT=/opt/uge-8.1.7p3/bin/lx-amd64/qstat
QALTER=/opt/uge-8.1.7p3/bin/lx-amd64/qalter
GERMS_DATA=<path_to_sequence_files>
STAGING=<path_to_staging_dir>

# usually BROWSER_BASE probably will be something like /usr/local
BROWSER_BASE=<path_to_install_location>
FQSNAKE=$BROWSER_BASE/lib/browser-fastq.snakefile
SNPSNAKE=$BROWSER_BASE/lib/browser-snpcall.snakefile
CONF=$BROWSER_BASE/etc/browser-snakemake.json
MYSQL_REF=$BROWSER_BASE/bin/mysql-get-reference.pl

# scale threads to the machine
THREADS=$(grep -c '^processor' /proc/cpuinfo)
# but don't go over some maximum (4 here)
if [ $THREADS -gt 4 ]; then
  THREADS=4
fi
