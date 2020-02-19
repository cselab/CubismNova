#!/bin/bash -l
#BSUB -W 00:40
#BSUB -J euler-div-kernel-512
#BSUB -o euler-div-kernel-512_<%J>.out
#BSUB -e euler-div-kernel-512_<%J>.err
#BSUB -n 36
#BSUB -N fabianw@mavt.ethz.ch
#disabled#BSUB -B fabianw@mavt.ethz.ch

sigdump()
{
    echo -e "$1" > "$SIGNALDUMP"
}

if [[ -e './common.sh' ]]
then
    source './common.sh'
    SIGNALDUMP="${JC_LOGHOME}/.<${LSB_JOBNAME}_${LSB_JOBID}>.tmp"
    touch "$SIGNALDUMP"
    trap_this sigdump SIGHUP SIGINT SIGTERM SIGKILL SIGUSR1 SIGUSR2 SIGPIPE
fi

PATH="$(pwd -P)":$PATH
export OMP_NUM_THREADS=1

ROOT="$(pwd -P)"
export RUN="${ROOT}/euler-$(date +%Y-%m-%d-%s)"
mkdir -p ${RUN}

source env2lmod.sh
module purge
module load gcc/8.2.0

cd ../legacy
make clean
make CC=g++
for (( i = 0; i < 10; i++ )); do
    ./legacy >> ${RUN}/legacy.gcc
done

cd ../nova
make clean
make CC=g++
for (( i = 0; i < 10; i++ )); do
    ./nova >> ${RUN}/nova.gcc
done

module load llvm/6.0.0

cd ../legacy
make clean
make CC=clang++
for (( i = 0; i < 10; i++ )); do
    ./legacy >> ${RUN}/legacy.clang
done

cd ../nova
make clean
make CC=clang++
for (( i = 0; i < 10; i++ )); do
    ./nova >> ${RUN}/nova.clang
done

if [[ -e './common.sh' ]]
then
    sigdump 0
fi

exit 0
