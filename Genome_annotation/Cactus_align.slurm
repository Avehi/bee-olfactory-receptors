#!/bin/bash
#PBS -l nodes=1:ppn=24:rhel7  #how many compute nodes; how many cores/node; what type of node
#PBS -l walltime=48:00:00     #how long to run; HH:MM:SS
#PBS -l pmem=10gb             #memory per core
#PBS -A open                  #what allocation to use
#PBS -N spades_denovo         #job name 
#PBS -j oe                    #print stdin/stderr to same file

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=24
export PATH=/gpfs/group/mml64/default/share/singularity/cactus-v1.2.3.img:/gpfs/group/mml64/default/share/singularity/spades-v3.15.2.img:$PATH

# The tree was broken up into clades and the script below was run for each clade

CLADE="clade1"
if [ ! -f ${CLADE}.hal ]; then
echo "
mkdir -p cactus_${CLADE}_work
cactus cactus_${CLADE} jobfile.cactus.txt ${CLADE}.hal --workDir cactus_${CLADE}_work --root ${CLADE} --maxCores ${OMP_NUM_THREADS}
halStats ${CLADE}.hal > ${CLADE}.stats
" > run_cactus_${CLADE}.sh
singularity exec /gpfs/group/mml64/default/share/singularity/cactus-v1.2.3.img /bin/sh run_cactus_${CLADE}.sh &>run_cactus_${CLADE}.log
cp ${CLADE}.hal /gpfs/group/mml64/default/share/
chgrp mml64_collab /gpfs/group/mml64/default/share/${CLADE}.hal
fi

echo "Finished"