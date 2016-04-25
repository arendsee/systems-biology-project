#!/usr/bin/env bash

usage (){
    echo "Retrieve all logo images from athamap"
    exit 0
}

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

outdir=INTERMEDIATE/logos

mkdir -p $outdir
ls INPUT/athamap/*txt |
    sed -r 's/.*\/([^.]+).*/\1/' |
    while read j
    do
        wget -O $outdir/$j.png http://www.athamap.de/Logos/pic${j}.png
    done
