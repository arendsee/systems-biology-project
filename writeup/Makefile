TARGET=report

${TARGET}.pdf: ${TARGET}.Rnw ${TARGET}.bib
	Rscript -e "library(knitr); knit('${TARGET}.Rnw')"
	latexmk -bibtex -pdf ${TARGET}.tex
	cp ${TARGET}.pdf ../OUTPUT

.PHONY: clean
clean:
	latexmk -f -c ${TARGET}.tex
	rm -rf cache ${TARGET}.tex ${TARGET}.bbl

.PHONY: deepclean
deepclean:
	latexmk -f -c ${TARGET}.tex
	rm -rf cache figure ${TARGET}.tex ${TARGET}.bbl ${TARGET}.pdf
