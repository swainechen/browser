#!/bin/bash
RUN=$1

# set some common variables
if [ -f /mnt/software/stow/GERMS_BROWSER/etc/browser-conf.sh ]; then
  . /mnt/software/stow/GERMS_BROWSER/etc/browser-conf.sh
fi

# first get some info
# $MYSQL_REF provides data as environment variables
# WORK is the working directory
# REF is the reference genome
# SPECIES is the short species name
# FULLSPECIES is the long species name
eval $( $MYSQL_REF $RUN )

if [ ! -d $WORK ]; then
  sleep 5
  echo "Current WORK ($WORK) doesn't exist."
  echo "Sleeping 5 to wait for database."
  eval $( $MYSQL_REF $RUN )
fi

RUN_TARGETS=$TARGETS
if [ $REF = "NONE" ]; then
  RUN_TARGETS=$NOREF_TARGETS
fi

CURRENTHOST=$( hostname -s )

case $CURRENTHOST in
aquilaln1|aquilaln2)
  if $QSTAT -j $RUN >& /dev/null; then
    true
  else
    cd $STAGING
    # we generally should have a species set
    # if we don't have a reference that will come out as NONE
    # this is when we want to do the no reference version
    $QSUB -pe OpenMP $THREADS -l mem_free=$MEM,h_rt=$TIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN "$SNAKE -d $WORK -s $SNPSNAKE $RUN_TARGETS -p -j $THREADS --ri --configfile $CONF --config ACC=$RUN REF=$REF SPECIES=$SPECIES FULLSPECIES=\"$FULLSPECIES\" 2>> $RUN.log && qrls $RUN-stage3"
    # we need to hold the cleanup for this new job
  fi
;;
*)
  if ssh aquila "$QSTAT -j $RUN >& /dev/null"; then
    true
  else
    ssh aquila "export GERMS_DATA=$GERMS_DATA && cd $STAGING && $QSUB -pe OpenMP $THREADS -l mem_free=$MEM,h_rt=$TIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN \"$SNAKE -d $WORK -s $SNPSNAKE $RUN_TARGETS -p -j $THREADS --ri --configfile $CONF --config ACC=$RUN REF=$REF SPECIES=$SPECIES FULLSPECIES=\"\\\"$FULLSPECIES\\\"\" 2>> $RUN.log && qrls $RUN-stage3\""
  fi
;;
esac
sleep 1
