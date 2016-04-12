mkdir -p data/logos
ls data/athamap/*txt |
    sed -r 's/.*\/([^.]+).*/\1/' |
    while read j
    do
        wget -O data/logos/$j.png http://www.athamap.de/Logos/pic${j}.png
    done
