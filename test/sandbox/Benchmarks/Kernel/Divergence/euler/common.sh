#!/bin/bash -l
# VARS (static after submission)
ROOT="$(pwd -P)"
export JC_LOGHOME="${ROOT}/01_jobs"
export JC_LOGHOME_LINKNAME='0000_joblog'
export JC_JOBNAME='euler-div-kernel-512'
export JC_COMMON="$JC_JOBWD/common.sh"
export JC_SUBMIT="$JC_JOBWD/submit.py"
export JC_PROLOG="$JC_JOBWD/prolog.sh"
export JC_EPILOG="$JC_JOBWD/epilog.sh"
export JC_JOB="$JC_JOBWD/job.sh"

# FUNCTIONS (exported at the end)
addToList()
{
    LIST="$1"; shift
    LID="$1"; shift # local ID
    JID="$1"; shift # job ID
    JNA="$1"; shift # job name
    JWD="$1"; shift # job working directory
    NOW=$(date '+%s') # time now

    LOCKFILE="${LIST}.lock"
    LINE="${LID},${JID},${JNA},${JWD},${NOW}"
    (flock -x 200
    echo "$LINE" >> "$LIST"
    ) 200> "$LOCKFILE"
}

removeFromList()
{
    LIST="$1"; shift
    LID="$1"; shift # local ID
    JID="$1"; shift # job ID
    JNA="$1"; shift # job name
    JWD="$1"; shift # job working directory

    LOCKFILE="${LIST}.lock"
    TMPFILE="${LIST}.tmp"
    PATTERN="${LID},${JID},${JNA},${JWD}"
    (flock -x 200
    cat "$LIST" | awk -v lid="$LID" -v jid="$JID" -v jna="$JNA" -v jwd="$JWD" 'BEGIN { FS="," }; ($1==lid) && ($2==jid) && ($3==jna) && ($4==jwd) { next }; { print $0 }' > "$TMPFILE"
    mv "$TMPFILE" "$LIST"
    ) 200> "$LOCKFILE"
}

backupTextFiles()
{
    list=`find "$JC_JOBWD" -maxdepth 1 -type f ! -iname "*~" ! -iname "*.xmf" -exec grep -Iq . {} \; -and -print`
    if [ ! -z "$list" ]
    then
        mkdir -p "$1/top_level_backup"
        cp -t "$1/top_level_backup" $list
    fi
}

trap_this()
{
    trap_exec="$1"; shift;
    for signal
    do
        trap "${trap_exec} ${signal}; exit" "$signal"
    done
}

export -f addToList
export -f removeFromList
export -f backupTextFiles
export -f trap_this
