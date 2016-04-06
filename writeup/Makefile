TARGET=writeup

${TARGET}.pdf: ${TARGET}.Rnw
	Rscript -e "library(knitr); knit('${TARGET}.Rnw')"
	latexmk -pdf ${TARGET}.tex

.PHONY: clean
clean:
	latexmk -c ${TARGET}.tex
	rm -rf cache ${TARGET}.tex
