#!/bin/bash

convert figure/degree_chart-1.pdf -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_degree.png
convert figure/graph_stats-1.pdf  -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_graph_stats.png
convert figure/non_postcut-1.pdf  -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_non_net-2.png
convert figure/non_precut-1.pdf   -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_non_net-1.png
convert figure/orph_postcut-1.pdf -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_net-2.png
convert figure/orph_precut-1.pdf  -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_net-1.png
convert figure/score-1.pdf        -density 300 -quality 100 -trim -sharpen 0x1.0 -resize 200% /home/rnz/Dropbox/sysbio-img_score.png
