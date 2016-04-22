#!/usr/bin/bash
set -u

sradir=$HOME/ncbi/public/sra

usage (){
cat << EOF
Process one run with kallisto starting with a run accession
REQUIRED ARGUMENTS
  -t IND  Transcriptome index (*.index) or fasta file (*.fa)
  -r ACC  A run accession
  -o OUT  Output folder
  -n SRA  NCBI public folder (default = $sradir)
  -d THR  Number of threads to use in Kallisto
  -c      Clean up files after using them
EOF
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

nthreads=1 clean=0
while getopts "ht:r:o:n:d:c" opt; do
    case $opt in
        h)
            usage ;;
        t) 
            trans=$OPTARG ;;
        r) 
            runid=$OPTARG ;;
        o) 
            outdir=$OPTARG ;;
        n) 
            sradir=$OPTARG ;;
        d)
            nthreads=$OPTARG ;;
        c)
            clean=1 ;;
    esac 
done

mkdir -p $outdir

# build an index of the transcriptome for use by kallisto
if [[ $trans =~ fa$ ]]
then
    echo "Building " $trans " index"
    trans_idx=$(sed 's/\.fa$/.index/' <<< $trans)
    kallisto index --index="$trans_idx" "$trans"
    trans="$trans_idx"
fi

time prefetch        \
    --max-size 100G  \
    --transport ascp \
    --ascp-path "/opt/aspera/bin/ascp|/opt/aspera/etc/asperaweb_id_dsa.openssh" \
    $runid


# Options descriptions
# split-files    - split paired-end data into files suffixed with _1 and _2
# readids        - append read id (.1, .2) after spot id
# dumpbase       - output as ACGT bases rather than color-base (e.g. from SOLiD)
# clip           - remove left and right tags
# skip-technical - skip technical reads (not useable by Kallisto, also is
#                  specific to Illumina multiplexing library construction
#                  protocol)
time fastq-dump        \
    --readids          \
    --split-files      \
    --dumpbase         \
    --skip-technical   \
    --clip             \
    --qual-filter-1    \
    --outdir "$outdir" \
    $sradir/${runid}.sra

# remove SRA archives from download folder
[[ $clean -eq 1 ]] && rm $sradir/${runid}*

# align reads to transcriptome index
kallisto_outdir=${outdir}/${runid}

if [[ -f ${runid}_2.sra ]]
then
    # paired end reads
    time kallisto quant                 \
        --index="$trans"                \
        --bootstrap-samples=100         \
        --output-dir="$kallisto_outdir" \
        --threads=$nthreads             \
        $outdir/${runid}_*
else
    # single reads
    time kallisto quant                 \
        --index="$trans"                \
        --bootstrap-samples=100         \
        --output-dir="$kallisto_outdir" \
        --single                        \
        --fragment-length=200           \
        --sd=20                         \
        --threads=$nthreads             \
        $outdir/${runid}_*
fi

# Remove fastq files
[[ $clean -eq 1 ]] && rm -rf $outdir/*fastq
