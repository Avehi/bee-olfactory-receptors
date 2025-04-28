
#!/bin/bash

set -e

PARAM_MAXHIT=5000

#-------generate hints for augustus cgp for a given chunk
#-------using output of find_candidate_intervals.sh
#-------run using catk singularity container

#using BED-style coordinates
SPECIES=${1}
CONTIG=${2}
LOWER=${3}
UPPER=${4}
PROTEIN_QUERY=${5}
HAL=${6}

#pull out chunk as maf and as reference fasta
#not using ancestral sequences right now in an effort to speed things up
CHUNK=${CONTIG}__${LOWER}__${UPPER}
mkdir -p maf/${SPECIES}/${CHUNK}
hal2maf --refGenome ${SPECIES} --refSequence ${CONTIG} --noAncestors --start $LOWER --length $((UPPER - LOWER)) \
  ${HAL} maf/${SPECIES}/${CHUNK}/maf
hal2fasta --sequence ${CONTIG} --start $LOWER --length $((UPPER - LOWER)) \
  ${HAL} ${SPECIES} >maf/${SPECIES}/${CHUNK}/fna

#prescreen with mmseqs2
#use top hits to reduce work for exonerate
util/mmseqs/bin/mmseqs easy-search --threads 1 --sort-results 1 --max-seqs $PARAM_MAXHIT \
   maf/${SPECIES}/${CHUNK}/fna ${PROTEIN_QUERY} maf/${SPECIES}/${CHUNK}/m8 maf/${SPECIES}/${CHUNK}/tmp
rm -rf maf/${SPECIES}/${CHUNK}/tmp
awk '{if(!($2 in _)) _[$2]=$12}END{for(i in _) print i"	"_[i]}' maf/${SPECIES}/${CHUNK}/m8 | \
  sort -k2,2nr | head -${PARAM_MAXHIT} | cut -f1 >maf/${SPECIES}/${CHUNK}/hits

if [ -s maf/${SPECIES}/${CHUNK}/hits ]; then

grep --no-group-separator -A1 -f maf/${SPECIES}/${CHUNK}/hits ${PROTEIN_QUERY} >maf/${SPECIES}/${CHUNK}/faa

#get gff via exonerate
exonerate --model protein2genome --showvulgar no --showalignment no --showtargetgff yes \
   maf/${SPECIES}/${CHUNK}/faa maf/${SPECIES}/${CHUNK}/fna >maf/${SPECIES}/${CHUNK}/gtf

#fix coordinates
awk -voffset=$LOWER 'BEGIN{FS="	"; OFS="	"} NR>2 && !($0 ~ /^#/) {$4=$4+offset; $5=$5+offset; print $0}' maf/${SPECIES}/${CHUNK}/gtf >maf/${SPECIES}/${CHUNK}/gtf.tmp

#convert to hints
exonerate2hints.pl --in=maf/${SPECIES}/${CHUNK}/gtf.tmp --source=XNT --CDSpart_cutoff=0 --out=maf/${SPECIES}/${CHUNK}/hints
sort -n -k 4,4 maf/${SPECIES}/${CHUNK}/hints | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl >maf/${SPECIES}/${CHUNK}/gtf.tmp
mv maf/${SPECIES}/${CHUNK}/gtf.tmp maf/${SPECIES}/${CHUNK}/hints

else
echo "No hits"
fi

echo "---FINISHED---"

