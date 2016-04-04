LANC=C
cut -f1,3 data/athamap-promoters.tab |
    sed 's/\.[0-9]\+//' |
    sort -u

# cut -f1,3,5 data/athamap-promoters.tab |
#     sed 's/\.[0-9]\+//' |
#     sort -t $'\t' -k1,1 -k3,3rn |
#     cut -f1,2 |
#     uniq
