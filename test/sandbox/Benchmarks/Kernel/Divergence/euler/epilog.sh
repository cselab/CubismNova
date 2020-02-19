#!/bin/bash -l
source './common.sh'

IFS=$'\r\n' GLOBIGNORE='*' command eval 'ARY=($(cat $1))'
JC_JOBDATE="${ARY[0]}"
JC_JOBTIME="${ARY[1]}"
JC_JOBID="${ARY[2]}"

# removeFromList "${JC_LOGHOME}/running" "$JC_JOBID" "$LSB_JOBID" "$LSB_JOBNAME" "$JC_JOBWD"
rm -f "$1"

SIGNALDUMP="${JC_LOGHOME}/.<${LSB_JOBNAME}_${LSB_JOBID}>.tmp"
JOBSIG="$(cat $SIGNALDUMP)"
rm -f "$SIGNALDUMP"

RUNTIME=$(($LSB_JOB_END_TIME-$LSB_JOB_START_TIME))

JOBLOG="${JC_JOBID}_${LSB_JOBNAME}_<${LSB_JOBID}>_<${JC_JOBTIME}>"
JOBLOGDIR="${JC_LOGHOME}/${JC_JOBDATE}/${JOBLOG}"
JOBINFO="${JOBLOGDIR}/${LSB_JOBNAME}.info"
cat << EOF >> "$JOBINFO"
-------------------------------------------------------------------------------
EPILOG_EXEC : `date`
JOB_SUB     : `date -d@$LSB_JOB_SUBMIT_TIME`
JOB_PEND    : `date -u -d@$LSB_JOB_PEND_TIME '+%T'`
JOB_START   : `date -d@$LSB_JOB_START_TIME`
JOB_END     : `date -d@$LSB_JOB_END_TIME`
JOB_RUNTIME : `date -u -d@$RUNTIME '+%T'`
JOB_STATUS  : $JOBSIG
-------------------------------------------------------------------------------
EOF

rm -f "$JC_JOBWD/${JC_LOGHOME_LINKNAME}"
cp "$JC_JOBWD/$LSB_OUTPUTFILE" "$JOBLOGDIR"
cp "$JC_JOBWD/$LSB_ERRORFILE" "$JOBLOGDIR"
ln -s "$JC_JOBWD" "$JOBLOGDIR/rundir"
# backupTextFiles "$JOBLOGDIR"

EXIT_STATE="SIGNAL=$JOBSIG"
JOBLOGDIR_EXIT="${JOBLOGDIR}.[${EXIT_STATE}]"
mv "$JOBLOGDIR" "$JOBLOGDIR_EXIT"
rm -f "$JC_LOGHOME/$JC_JOBDATE/current"
ln -s "$JOBLOGDIR_EXIT" "$JC_LOGHOME/$JC_JOBDATE/current"
ln -s "$JOBLOGDIR_EXIT" "$JC_JOBWD/${JC_LOGHOME_LINKNAME}"

exit 0
