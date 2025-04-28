
#!/bin/bash

set -e

module load metaeuk
mkdir -p mte

SPECIES=$1
BUFFER=$2
PROTEINS=$3

metaeuk easy-predict fna/$SPECIES.fna ${PROTEINS} mte/${SPECIES} mte/${SPECIES}-tmp
rm -r mte/${SPECIES}-tmp

#coordinates are 0-based on both ends of interval, change to GFF-style
grep ">" mte/$SPECIES.fas | \
  awk 'BEGIN{FS="|"}{sub(">","",$1); print $2"\tmetaeuk\tgene\t"$7+1"\t"$8+1"\t.\t"$3"\t.\tQUERY="$1";BITSCORE="$4";EVAL="$5";EXONS="$6";"}' | \
  sort -k1,1d -k4,4n >mte/$SPECIES.gff

#add buffer and merge overlaps
#switching to bed-style coordinates
module load gcc/8.3.1
module load bedtools
module load samtools
if [ ! -f fna/${SPECIES}.fna.fai ]; then
  samtools faidx fna/${SPECIES}.fna
fi
while read CONTIG SOURCE TYPE START END STUFF; do
  CONTIG_LENGTH=`grep -w "${CONTIG}" fna/${SPECIES}.fna.fai | cut -f2`
  LOWER=$((START - BUFFER - 1))
  UPPER=$((END + BUFFER))
  if [ $LOWER -lt 0 ]; then LOWER=0; fi
  if [ $UPPER -gt $CONTIG_LENGTH ]; then UPPER=$CONTIG_LENGTH; fi
  if [ $((UPPER - LOWER)) -ge 2000 ]; then echo -e "$CONTIG\t$LOWER\t$UPPER"; fi
#done<mte/$SPECIES.gff | sort | uniq >mte/$SPECIES.bed #keep overlapping intervals separate
done<mte/$SPECIES.gff | bedtools merge -i - >mte/$SPECIES.bed #merge overlappling intervals

echo "---FINISHED---"

