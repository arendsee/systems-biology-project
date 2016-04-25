#!/usr/bin/env bash

usage (){
cat << EOF
This is a wrapper for process-sra.sh.
Requires 1 argument: a file of SRA runids
EOF
    exit 0
}

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

[[ -z $1 ]] && usage

outdir=INTERMEDIATE/kallisto-runs
index=INTERMEDIATE/at.index

while read id
do
    echo $id
    ./process-sra.sh -t $index -r $id -o $outdir -d 8 -c
done < $1
