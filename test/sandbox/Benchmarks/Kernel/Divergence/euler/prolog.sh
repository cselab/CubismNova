#!/bin/bash -l
source './common.sh'

IFS=$'\r\n' GLOBIGNORE='*' command eval 'ARY=($(cat $1))'
JC_JOBDATE="${ARY[0]}"
JC_JOBTIME="${ARY[1]}"
JC_JOBID="${ARY[2]}"

# removeFromList "${JC_LOGHOME}/pending" "$JC_JOBID" 'pend' "$LSB_JOBNAME" "$JC_JOBWD"

JOBLOG="${JC_JOBID}_${LSB_JOBNAME}_<${LSB_JOBID}>_<${JC_JOBTIME}>"
JOBLOGDIR="${JC_LOGHOME}/${JC_JOBDATE}/${JOBLOG}"
JOBINFO="${JOBLOGDIR}/${LSB_JOBNAME}.info"
mkdir -p "$JOBLOGDIR"
rm -f "$JC_LOGHOME/current"
rm -f "$JC_LOGHOME/$JC_JOBDATE/current"
ln -s "$JOBLOGDIR" "$JC_LOGHOME/$JC_JOBDATE/current"
ln -s "$JC_LOGHOME/$JC_JOBDATE/current" "$JC_LOGHOME/current"
cat << EOF > "$JOBINFO"
-------------------------------------------------------------------------------
PROLOG_EXEC : `date`
JOB_NAME    : $LSB_JOBNAME
JOB_ID      : $LSB_JOBID
JOB_QUEUE   : $LSB_QUEUE
SUB_USER    : $LSB_SUB_USER
SUB_CWD     : $LS_SUBCWD
SUB_HOST    : $LSB_SUB_HOST
SUB_RES_REQ : $LSB_SUB_RES_REQ
EXEC_HOSTS  : $LSB_MCPU_HOSTS
EOF

rm -f "${JC_JOBWD}/${JC_LOGHOME_LINKNAME}"
ln -s "$JOBLOGDIR" "${JC_JOBWD}/${JC_LOGHOME_LINKNAME}"

# addToList "${JC_LOGHOME}/running" "$JC_JOBID" "$LSB_JOBID" "$LSB_JOBNAME" "$JC_JOBWD"
exit 0