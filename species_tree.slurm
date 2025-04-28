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

# Here, we assume that the BUSCO output files for each genome are in the current working directory.
# get list of gene names of complete genes
for file in $(find . -iname "full_table.tsv"); do
  grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> complete_busco_ids.txt;
done
# filter out BUSCOs that are present in less than 3 genomes
sort complete_busco_ids.txt |uniq -c > complete_busco_ids_with_counts.txt
awk '$1 > 2 {print $2}' complete_busco_ids_with_counts.txt > final_busco_ids.txt

# put all amino acid and nucleotide sequences into a single folder and add species names into the fasta headers
mkdir -p busco_nt
mkdir -p busco_aa
for dir in $(find . -type d -name "single_copy_busco_sequences"); do 
  sppname=$( dirname $dir | cut -f 1-2 -d "_" | sed 's@./@@g');
  for file in ${dir}/*.faa; do
    filename=`basename ${file}`
    cp $file busco_aa/${sppname}_${filename}
    sed -i 's/^>/>'${sppname}'|/g' busco_aa/${sppname}_${filename}
      cut -f 1 -d ":" busco_aa/${sppname}_${filename} | tr '[:lower:]' '[:upper:]' > busco_aa/${sppname}_${filename}.1
  	mv busco_aa/${sppname}_${filename}.1 busco_aa/${sppname}_${filename}
  done
  for file in ${dir}/*.fna; do
    filename=`basename ${file}`
    cp $file busco_nt/${sppname}_${filename}
    sed -i 's/^>/>'${sppname}'|/g' busco_nt/${sppname}_${filename}
      cut -f 1 -d ":" busco_nt/${sppname}_${filename} | tr '[:lower:]' '[:upper:]' > busco_nt/${sppname}_${filename}.1
	    mv busco_nt/${sppname}_${filename}.1 busco_nt/${sppname}_${filename}
	done
done

# sort all files by gene
while read line; do
cat busco_aa/*_${line}.faa >> ${line}_aa.fasta;
cat busco_nt/*_${line}.fna >> ${line}_nt.fasta;
done<final_busco_ids.txt

# make alignments
# mafft script maybe run in different folder
for file in *_aa.fasta; do
  mafft.bat ${I} > ${I%.fasta}_aln.fasta
done

for file in *_nt.fasta; do
  mafft.bat ${I} > ${I%.fasta}_aln.fasta
done

#change to new directory too
cd aa_alignment
for aln in *.aln; do
  iqtree -s $aln
done

cd nt_alignment
for aln in *.aln; do
  iqtree -s $aln
done

# make a single file with all ML trees -- check this code first
cat *.treefile >> busco_genes_ML_nt.tree 

# remove busco gene IDs
# pick qany one RAXML best tree from teh raxml_aa folder
sed -i 's/GENOMIC.FNA[^:]*/GENOMIC.FNA/g' busco_genes_ML_nt.tree
sed -i 's/VANCOUVERENSIS[^:]*/VANCOUVERENSIS_NEARCTICUS.GENOMIC.FNA/g' busco_genes_ML_nt.tree
sed -i 's/VENTRALIS_[^:]*/VENTRALIS_HOOSANA.GENOMIC.FNA/g' busco_genes_ML_nt.tree
q_nt.tree

# make species tree
~/Biotools/ASTER/bin/astral -C -i busco_genes_ML_aa.trees -c busco_genes_ASTRAL_aa_spp.tree -o busco_genes_ASTRAL_aa_annotated.tree > annotated.tre
python3 ~/Biotools/CASTLES-1.3.0/castles.py -t annotated.tre -g busco_genes_ML_aa.trees -o CASTLES_busco_genes_aa_spp.tree

~/Biotools/ASTER/bin/astral -C -i busco_genes_ML_nt.trees -c busco_genes_ASTRAL_nt_spp.tree -o busco_genes_ASTRAL_nt_annotated.tree > annotated.tre
python3 ~/Biotools/CASTLES-1.3.0/castles.py -t annotated.tre -g busco_genes_ML_nt.trees -o CASTLES_busco_genes_nt_spp.tree

