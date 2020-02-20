#!/bin/bash -l
#BSUB -W 01:00
#BSUB -J euler-div-kernel
#BSUB -o euler-div-kernel_<%J>.out
#BSUB -e euler-div-kernel_<%J>.err
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
POST="${ROOT}/.."
export RUN="${ROOT}/euler-$(date +%Y-%m-%d-%s)"
mkdir -p ${RUN}

source env2lmod.sh
module purge
module load gcc/8.2.0
module load python/3.6.4

blocks=(1 2 4 8 16)

cd ../legacy
make clean
make CC=g++
for b in ${blocks[@]}; do
    for (( i = 0; i < 10; i++ )); do
        ./legacy $b >> ${RUN}/legacy.$b.gcc
    done
done
python ${POST}/post.py --files ${RUN}/legacy.*.gcc --output ${RUN}/legacy.gcc

cd ../nova
make clean
make CC=g++
for b in ${blocks[@]}; do
    for (( i = 0; i < 10; i++ )); do
        ./nova $b >> ${RUN}/nova.$b.gcc
    done
done
python ${POST}/post.py --files ${RUN}/nova.*.gcc --output ${RUN}/nova.gcc

module load llvm/6.0.0

cd ../legacy
make clean
make CC=clang++
for b in ${blocks[@]}; do
    for (( i = 0; i < 10; i++ )); do
        ./legacy $b >> ${RUN}/legacy.$b.clang
    done
done
python ${POST}/post.py --files ${RUN}/legacy.*.clang --output ${RUN}/legacy.clang

cd ../nova
make clean
make CC=clang++
for b in ${blocks[@]}; do
    for (( i = 0; i < 10; i++ )); do
        ./nova $b >> ${RUN}/nova.$b.clang
    done
done
python ${POST}/post.py --files ${RUN}/nova.*.clang --output ${RUN}/nova.clang

if [[ -e './common.sh' ]]
then
    sigdump 0
fi

exit 0
