#!/usr/bin/env bash

for f in ../data/logos/*png
do
    m=$(basename $f | sed 's/\.png//')
    if [[ ! -r $m.pdf ]]
    then
        sed "s/MOTIF_NAME/$m/" ../plot-promoters.Rnw > motif.Rnw
        make
        mv motif.pdf $m.pdf
        rm -rf motif* figure
    fi
done

pdfunite *pdf final.pdf
