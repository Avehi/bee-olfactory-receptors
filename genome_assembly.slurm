{bash}
#!/bin/bash
#PBS -l nodes=1:ppn=24:rhel7  #how many compute nodes; how many cores/node; what type of node
#PBS -l walltime=48:00:00     #how long to run; HH:MM:SS
#PBS -l pmem=10gb             #memory per core
#PBS -A open                  #what allocation to use
#PBS -N spades_denovo         #job name 
#PBS -j oe                    #print stdin/stderr to same file

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=24
set -e

INPUT_DIR=/storage/home/aps6366/scratch/genome_assembly_avehi/Mdes
WORK_DIR=/storage/home/aps6366/scratch/genome_assembly_avehi/Mdes/work
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

free >memory.info

# 1. PREPROCESSING
module purge
module load bbtools
ADAPTERS=$GROUP/../.local/bbtools/resources/adapters.fa
BBDUK="hdist=1 ktrim=r k=23 mink=11 qtrim=rl trimq=10 tbo=t tbe=t minlen=70"
BBDUK=${BBDUK}" ref=${ADAPTERS}"
FQ1=${WORK_DIR}/INPUT_1.fastq
FQ2=${WORK_DIR}/INPUT_2.fastq
FQDEDUP1=${WORK_DIR}/INPUT_1.dedup.fastq
FQDEDUP2=${WORK_DIR}/INPUT_2.dedup.fastq
FQTRIM1=${WORK_DIR}/INPUT_1.trim.fastq
FQTRIM2=${WORK_DIR}/INPUT_2.trim.fastq
TEMPDIR=${WORK_DIR}/bbtools-clumpify-temp
if [ ! -f preprocess.finished ]; then
  zcat ${INPUT_DIR}/*_1.fq.gz >${FQ1}
  zcat ${INPUT_DIR}/*_2.fq.gz >${FQ2}
  clumpify.sh t=${OMP_NUM_THREADS} Xmx=100g usetmpdir=true tmpdir=$TEMPDIR in1=${FQ1} in2=${FQ2} out1=${FQDEDUP1} out2=${FQDEDUP2} dedupe reorder &>${W
ORK_DIR}/clumpify.log
  bbduk.sh in1=${FQDEDUP1} in2=${FQDEDUP2} out1=${FQTRIM1} out2=${FQTRIM2} ${BBDUK} &>${WORK_DIR}/bbduk.log
  rm -r ${FQDEDUP1} ${FQDEDUP2} ${FQ1} ${FQ2} ${TEMPDIR}
  touch preprocess.finished
fi

# 2. ASSEMBLY
module purge
module load spades
SPADES_WORKDIR=$WORK_DIR/spades
SPADES_CONTIGS=$SPADES_WORKDIR/contigs.fasta
SPADES_SCAFFOLDS=$SPADES_WORKDIR/scaffolds.fasta
if [ ! -f spades.finished ]; then
  if [ ! -d ${SPADES_WORKDIR} ]; then
    spades.py -t ${OMP_NUM_THREADS} -1 ${FQTRIM1} -2 ${FQTRIM2} --checkpoints all -o ${SPADES_WORKDIR} &>$WORK_DIR/spades.log
  else
    spades.py --continue -o ${SPADES_WORKDIR} &>>$WORK_DIR/spades.log
  fi
  touch spades.finished
fi

# 3. CORRECTION
module purge
module load ragtag
REF=$WORK_DIR/eucera_pruinosa.fasta
RAGTAG_SCAFFOLDS_IN=$WORK_DIR/spades_scaffolds.fa
RAGTAG_CORRECT_SCAFF=$WORK_DIR/ragtag_correct
if [ ! -f correction.finished ]; then
  zcat /storage/home/nsp5229/group/nsp5229/pepo/genome/15JUL2020/eucera*.fasta.gz >${REF}
  cat ${SPADES_SCAFFOLDS} | bioawk -cfastx 'length($seq) >= 1000 {print ">"$name"\n"$seq}' >${RAGTAG_SCAFFOLDS_IN}
  samtools faidx ${RAGTAG_SCAFFOLDS_IN}
  echo ${FQTRIM1} >reads.fofn
  echo ${FQTRIM2} >>reads.fofn
  ragtag.py correct -t ${OMP_NUM_THREADS} -o ${RAGTAG_CORRECT_SCAFF} --debug --remove-small -f 1000 -F reads.fofn -T sr ${REF} ${RAGTAG_SCAFFOLDS_IN} &
>$WORK_DIR/ragtag_correct.log
  touch correction.finished
fi

# 4. SCAFFOLDING
module purge
module load ragtag
RAGTAG_CORRECTED=$RAGTAG_CORRECT_SCAFF/ragtag.correct.fasta
RAGTAG_SCAFFOLD=$WORK_DIR/ragtag_scaffold
if [ ! -f $WORK_DIR/scaffolding.finished ]; then
  ragtag.py scaffold -t ${OMP_NUM_THREADS} -o ${RAGTAG_SCAFFOLD} ${REF} ${RAGTAG_CORRECTED} &>$WORK_DIR/ragtag_scaffold.log
  touch scaffolding.finished
fi

# 6. CONTAMINATION
if [ ! -f $RAGTAG_SCAFFOLD/ragtag.scaffold.coverage.png ]; then
  bioawk -cfastx '{print $name}' $RAGTAG_SCAFFOLD/ragtag.scaffold.fasta | perl -pe 's/(.+)_length_(.+)_cov_(.+)\..+/$1 $2 $3/g' | perl -pe 's/:.+//g' | awk 'NF==3' >$RAGTAG_SCAFFOLD/ragtag.scaffold.coverage
  echo "
  library(ggplot2)
  cov<-read.table('$RAGTAG_SCAFFOLD/ragtag.scaffold.coverage')
  plt<-ggplot(cov, aes(x=log10(V2),y=log10(V3), col=V3<4|V2<1000)) + geom_point()
  ggsave('$RAGTAG_SCAFFOLD/ragtag.scaffold.coverage.png', plt, height=8, width=8)" | R --slave
  filtering coverage >=4 and length >= 1000 seems to get rid of most of the rubbish
fi

FILTERED_SCAFFOLD=$WORK_DIR/ragtag_scaffolds.filtered.fa
if [ ! -f $FILTERED_SCAFFOLD ]; then
  # filtering based on given criteria
  bioawk -cfastx '{print ">"$name,$seq}' $RAGTAG_SCAFFOLD/ragtag.scaffold.fasta | \
    gawk '$1 ~ /epru/ {print $1"\n"$2} !($1 ~ /epru/) {match($1, /length_([0-9]+)/, ln); match($1, /cov_([0-9]+)/, cv); if(cv[1]>=4 && ln[1]>=1000) print $1"\n"$2}' | \
    fold >$FILTERED_SCAFFOLD
fi

if [ ! -f $FILTERED_SCAFFOLD.sketch ]; then
  # assess contamination via minhash (excluding insects)
  module load bbtools
  stats.sh in=$FILTERED_SCAFFOLD >$FILTERED_SCAFFOLD.stats
  sendsketch.sh address=nt size=100000 records=50 colors=f in=$FILTERED_SCAFFOLD exclude=50557 >$FILTERED_SCAFFOLD.sketch
  sendsketch.sh address=nt size=100000 records=50 colors=f in=$FILTERED_SCAFFOLD exclude=50557 include=2 >$FILTERED_SCAFFOLD.sketch.bacterial
  sendsketch.sh address=nt size=100000 records=50 colors=f in=$FILTERED_SCAFFOLD exclude=50557 include=4751 >$FILTERED_SCAFFOLD.sketch.fungal
  sendsketch.sh address=nt size=100000 records=50 colors=f in=$RAGTAG_SCAFFOLD/ragtag.scaffold.fasta exclude=50557 >$FILTERED_SCAFFOLD.sketch.prefilter
fi

# 5. POLISHING
POLISH=$WORK_DIR/nextpolish
SORT_THREADS=3
MAP_THREADS=$((OMP_NUM_THREADS - SORT_THREADS - 1))
INPUT=$POLISH/genome.polish.fa
FILTERED_POLISHED=$WORK_DIR/ragtag_scaffolds.filtered.polished.fa
if [ ! -f $FILTERED_POLISHED ]; then
  mkdir -p $POLISH
  #NB: change header prefix \/ to correct 4-letter abbreviation, here "eang" is for Eucera angelica
  bioawk -cfastx '{print ">mdes_"NR"\n"$seq}' $FILTERED_SCAFFOLD | fold >$INPUT
  for ROUND in 0 1 2 3 4; do
      if [ ! -f ${POLISH}/round${ROUND}.done ]; then for STEP in 1 2; do
        module purge
        module load minimap2
        module load samtools
        # map reads
        minimap2 -t ${MAP_THREADS} -ax sr ${INPUT} ${FQTRIM1} ${FQTRIM2} 2>${POLISH}/minimap.log | \
          samtools view -@ ${SORT_THREADS} -F 0x4 -b - | \
          samtools sort -m 2G -@ ${SORT_THREADS} -o ${POLISH}/sgs.sort.bam -
        samtools index -@ ${MAP_THREADS} ${POLISH}/sgs.sort.bam
        samtools stats ${POLISH}/sgs.sort.bam >${POLISH}/stats.round${ROUND}.step${STEP}.txt
        # polish
        samtools faidx ${INPUT}
        module purge
        module load nextpolish
        nextpolish1.py -g ${INPUT} -t ${STEP} -p ${OMP_NUM_THREADS} -s ${POLISH}/sgs.sort.bam 1>${POLISH}/genome.polishtemp.fa 2>${POLISH}/nextpolish.log
        mv ${INPUT} ${POLISH}/genome.backup.round${ROUND}.step${STEP}.fa
        mv ${POLISH}/genome.polishtemp.fa ${INPUT}
      done; touch ${POLISH}/round${ROUND}.done; fi
  done
  cp $INPUT $FILTERED_POLISHED
  grep "error" ${POLISH}/stats* >nextpolish.stats
  rm -r ${POLISH}
  module purge
fi

# 6. GAP CLOSING
module purge
module load abyss
SEAL_DIR=$WORK_DIR/seal
FILTERED_SEALED=$WORK_DIR/ragtag_scaffolds.filtered.polished.sealed.fa
if [ ! -f $FILTERED_SEALED ]; then
  mkdir -p $SEAL_DIR
  bioawk -cfastx '{print ">"$name"\n"toupper($seq)}' ${FILTERED_POLISHED} | perl -pe 's/^N+//g' | perl -pe 's/N+$//g' | fold >$SEAL_DIR/genome.fa
  abyss-sealer -b20G -k48 -k64 -k80 -k96 -k112 -k128 \
    -o ${SEAL_DIR}/genome.sealed \
    -S ${SEAL_DIR}/genome.fa \
    -j ${OMP_NUM_THREADS} \
    -L 150 \
    --verbose \
    ${FQTRIM1} ${FQTRIM2} &>$SEAL_DIR/sealer.log
  mv ${SEAL_DIR}/genome.sealed_scaffold.fa $FILTERED_SEALED
  mv ${SEAL_DIR}/genome.sealed_log.txt ${FILTERED_SEALED//fa/log}
  rm -r ${SEAL_DIR}
fi

## 7. BUSCO
module purge
module load anaconda3
source activate ~/work/busco
if [ ! -f $WORK_DIR/busco.finished ]; then
  busco -m genome -c ${OMP_NUM_THREADS} -i ${FILTERED_SEALED} -o busco -l hymenoptera_odb10 --augustus_species honeybee1 &>busco.log
  cp $WORK_DIR/busco/short_summary* $WORK_DIR/final.fasta.busco
  cp ${FILTERED_SEALED} $WORK_DIR/final.fasta
  gzip $WORK_DIR/final.fasta
  touch $WORK_DIR/busco.finished
fi
conda deactivate 
