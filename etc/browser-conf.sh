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
STAGING=

BROWSER_BASE=
FQSNAKE=$BROWSER_BASE/lib/browser-fastq.snakefile
SNPSNAKE=$BROWSER_BASE/lib/browser-snpcall.snakefile
CONF=$BROWSER_BASE/etc/browser-snakemake.json
MYSQL_REF=$BROWSER_BASE/bin/mysql-get-reference.pl

THREADS=4
# GIS cluster config - for main jobs
MEM="32G"
TIME="48:00:00"
# for Fastq and downloading
FMEM="4G"
FTIME="2:00:00"
