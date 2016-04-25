#!/usr/bin/env bash

datadir=$PWD/INPUT/athamap
logodir=$PWD/INTERMEDIATE/logos
outdir=$PWD/OUTPUT/motif-reports
rcode=$PWD/plot-motif.Rnw

[[ -d $outdir ]] || mkdir -p $outdir

cp motif-Makefile $outdir/Makefile

cd $outdir 

for f in $logodir/*png
do
    m=$(basename $f | sed 's/\.png//')
    logopath=$f
    datapath=$datadir/${m}.txt
    if [[ ! -r $m.pdf ]]
    then
        sed "s;LOGOPATH;$f;" $rcode |
            sed "s;DATAPATH;$datapath;" > motif.Rnw
        make
        mv motif.pdf $m.pdf
        rm -rf motif* figure
    fi
done

pdfunite *pdf final.pdf
