while read id
do
    echo $id
    ./process-sra.sh -t data/at.index -r $id -o output -d 8 -c
done < $1
