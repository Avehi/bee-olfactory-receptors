
#!/bin/bash

set -e

echo "`date`"

WORKDIR=$1
DB=$2
CFG=$3

#"exactlyone" requires a complete gene which might not be possible :-|
#could also use --alternatives-from-sampling to try to get a bunch of possible transcripts for debugging purposes

if [ -s ${WORKDIR}/maf ]; then
grep "# hal" ${WORKDIR}/maf | sed 's/# hal //' >${WORKDIR}/tree
augustus --species=honeybee1 --softmasking=1 --treefile=${WORKDIR}/tree --alnfile=${WORKDIR}/maf --dbaccess=$DB \
  --allow_hinted_splicesites=atac --extrinsicCfgFile=$CFG --dbhints=1 \
  --uniqueGeneId=true \
  --speciesfilenames=db/genomes.tbl --alternatives-from-evidence=0 --/CompPred/outdir=${WORKDIR}/cgp
fi

echo "---FINISHED---"

