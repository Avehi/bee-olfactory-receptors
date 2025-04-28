#/usr/bin/bash

# Step 0: Use R to create files with names for protein sequences for each orthogroup
# I generated these files and put them into a folder called "subfamily HMMER"

# Step 1: concatenate all the protein sequences into one big file (from initial set of alignments Nate sent)
cat /home/avehi/Work/sensory_orthogroups_22DEC21/OR/*.fa > OR_all_protein_alignments.fa

# Step 2: Pull out the protein sequence for each name and make a fasta file of each orthogroup
for FILE in *.txt; do
	seqtk subseq OR_all_protein_alignments.fa $FILE > $(basename $FILE .txt).fa;
done

# Step 3: Align fasta files with MAFFT
for I in `seq 1 1 46`; do
	/usr/bin/mafft  --genafpair  --maxiterate 16 --inputorder OR_subfamily_labels_${I}.fa > OR_subfamily_${I}_aligned.fa;
done

# Step 4: convert alignments to stockholm format
for I in *aligned.fa; do
	/home/avehi/Biotools/easel/miniapps/esl-reformat stockholm $I > $(basename $I .fa).stk;
done

# Step 5: run HMMER on the alignment to create a protein profile per subfamily
for I in *aligned.stk; do
	hmmbuild --wblosum $(basename $I _aligned.stk).hmm $I;
done

# Step 6: align sequences back to HMMs to produce final alignments
for I in `seq 1 1 46`; do 
	hmmalign --outformat afa OR_subfamily_${I}.hmm OR_subfamily_${I}_aligned.stk > OR_subfamily_${I}_final.fa; 
done

# Step 7: put all final alignments into a new folder
mkdir -p tree 
cd tree 
cp ../*final.fa .

# Step 8: filter out insertions (lowercase letters and .)
for I in `seq 1 1 46`; do
	 bioawk -cfastx '{gsub(/[a-z.]+/,"", $seq); print ">" $name "\n" $seq}' OR_subfamily_${I}_final.fa > OR_subfamily_${I}_noindel.fa;
done

# Step 9: remove gappy sequences
for I in `seq 1 1 46`; do
	trimal -in OR_subfamily_${I}_noindel.fa -out OR_subfamily_${I}_trim.fa -gappyout;
done
# trim header labels from fasta
bioawk -cfastx '{print ">"$name"\n"$seq}' Dmel.OR.fasta > out.fasta

# Step 10: align drosophila sequences to the OR subfamily alignments
for I in `seq 1 1 46`; do
	mafft  --keeplength  --add out.fasta OR_subfamily_${I}_trim.fa > OR_subfamily_${I}_dmel.fa;
done

#Step 11: make trees
for I in `seq 1 1 46`; do
    iqtree -s OR_subfamily_${I}_dmel.fa -m MFP -bb 1000 -nt AUTO;
done