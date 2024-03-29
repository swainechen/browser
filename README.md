# browser
Bacterial genome browser scripts

Most of these heavily depend on specifics of the cluster / servers I'm using at GIS.

# Old docs - a little out of date
GERMS browser
Version 2
SLC 171228

Main purpose is to eliminate any manual intervention with running assembly and SNP calling on sequencing data.

General workflow:
1. Prepping FASTQ files
  -- Download from GIS or Genbank SRA
  -- Collect some info and run Kraken
2. Running WGS and SNP calling
  -- Based on Kraken result, select a reference genome
  -- Run WGS, SNP caling, SRST2 and put the files in the right place
  -- Do a database import of everything
3. Cleaning up
  -- Check for completion
  -- Delete the FASTQ files

All meant to be stowed in /mnt/software/stow/GERMS_BROWSER
Subdirectories and files:
  bin
    browser-staging.sh
    browser-species.sh
  lib
    browser-fastq.snakefile
    browser-snpcall.snakefile
    browser-database.schema
  etc
    browser-conf.sh
    browser-snakemake.json

Generally should be run as:
for i in <list of IDs>; do
  while [ `df $GERMS_DATA | grep -v '^Filesystem' | awk '{print $3}' ` -lt 50000000 ]; do
    # this checks we have more than 50G available
    sleep 600
  done
  /mnt/software/stow/GERMS_BROWSER/bin/browser-staging.sh $i && date
  sleep 60
done
