#!/bin/bash
set -e

SPP=$1

##filtering criterion (for first pass):
#-no single exon (-U)
#-must be complete (-J)

grep "/$SPP/" pred/$SPP.list >pred/$SPP.list.primary

joingenes -f pred/$SPP.list.primary -j -o pred/$SPP.primary.gtf
fix_joingenes_gtf.pl < pred/$SPP.primary.gtf >pred/$SPP.primary.fixed.gtf
gtf2gff.pl --gff3 --includeStopInCDS <pred/$SPP.primary.fixed.gtf --out=pred/$SPP.primary.gff3
gffread -g fna/$SPP.fna -x pred/$SPP.primary.cds -y pred/$SPP.primary.aa -S -J -U pred/$SPP.primary.gff3

##this make be useful in a second pass to detect pseudogenes: combine secondary predictions
#joingenes -f pred/$SPP.list -j -o pred/$SPP.gtf
#fix_joingenes_gtf.pl < pred/$SPP.gtf >pred/$SPP.fixed.gtf
#gtf2gff.pl --gff3 --includeStopInCDS <pred/$SPP.fixed.gtf --out=pred/$SPP.gff3
#gffread -g fna/$SPP.fna -x pred/$SPP.cds -y pred/$SPP.aa -S pred/$SPP.gff3

rm pred/$SPP*gtf

