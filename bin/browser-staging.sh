#!/bin/bash
RUN=$1

# get directory where the script sits - needed to find the config file
# because the script might be linked, but the config file should be edited
# so we want the final link location of the script
# from https://unix.stackexchange.com/questions/17499/get-path-of-current-script-when-executed-through-a-symlink
DIR="$(cd "$(dirname "$0")" && pwd)"

# set some common variables
if [ -f $DIR/../etc/browser-conf.sh ]; then
  . $DIR/../etc/browser-conf.sh
fi

CURRENTHOST=$( hostname -s )

# we use $FMEM and $FTIME here for cluster jobs - should get us the short queue

# check if on AWS - https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/identify_ec2_instances.html
if GET http://169.254.169.254/latest/dynamic/instance-identity/ > /dev/null; then
  cd $STAGING
  browser-track.pl $RUN start
  $SNAKE -s $FQSNAKE fastq -p -j 1 --ri --configfile $CONF --config ACC=$RUN 2>> $RUN.log
  $BROWSER_BASE/bin/browser-species.sh $RUN 2>> $RUN.log
  tail -n 1 $RUN.log | grep -q 100...done && $SNAKE -s $FQSNAKE cleanup -p -j 1 --ri --configfile $CONF --config ACC=$RUN 2>> $RUN.log
  browser-track.pl $RUN finish
else
  case $CURRENTHOST in
  aquilaln1|aquilaln2)
    if $QSTAT -j $RUN >& /dev/null || $QSTAT -j $RUN-* >& /dev/null; then
      true
    else
      cd $STAGING
      # we should also really check for completion first here before doing this
      $QSUB -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage1 "$SNAKE -s $FQSNAKE fastq -p -j 1 --ri --configfile $CONF --config ACC=$RUN 2>> $RUN.log"
      # this next command starts a job with name $RUN
      $QSUB -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage2 -hold_jid $RUN-stage1 "$BROWSER_BASE/bin/browser-species.sh $RUN 2>> $RUN.log"
      # this needs better checking for completion - we wait first on $RUN
      # but $RUN doesn't start til the other job does as well
      # the other hold will have to get set by qalter in browser-species.sh
      $QSUB -h -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage3 -hold_jid $RUN-stage2 "tail -n 1 $RUN.log | grep -q 100...done && $SNAKE -s $FQSNAKE cleanup -p -j 1 --ri --configfile $CONF --config ACC=$RUN 2>> $RUN.log"
    fi
  ;;
  *)
    if ssh aquila "$QSTAT -j $RUN >& /dev/null" || ssh aquila "$QSTAT -j $RUN-* >& /dev/null"; then
      true
    else
      ssh aquila "export GERMS_DATA=$GERMS_DATA && cd $STAGING && $QSUB -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage1 \"$SNAKE -s $FQSNAKE fastq -p -j 1 --ri --configfile $CONF --config ACC=$RUN 2>> $RUN.log\""
      ssh aquila "export GERMS_DATA=$GERMS_DATA && cd $STAGING && $QSUB -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage2 -hold_jid $RUN-stage1 \"$BROWSER_BASE/bin/browser-species.sh $RUN\""
      ssh aquila "export GERMS_DATA=$GERMS_DATA && cd $STAGING && $QSUB -h -pe OpenMP 1 -l mem_free=$FMEM,h_rt=$FTIME -cwd -V -b y -o /dev/null -e /dev/null -N $RUN-stage3 -hold_jid $RUN \"tail -n 1 $RUN.log | grep -q 100...done && $SNAKE -s $FQSNAKE cleanup -p -j 1 --ri --configfile $CONF --config ACC=$RUN -p 2>> $RUN.log\""
    fi
  ;;
  esac
fi
sleep 1
