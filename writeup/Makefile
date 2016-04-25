TARGET=report

${TARGET}.pdf: ${TARGET}.Rnw ${TARGET}.bib
	Rscript -e "library(knitr); knit('${TARGET}.Rnw')"
	latexmk -bibtex -pdf ${TARGET}.tex
	cp ${TARGET}.pdf ../OUTPUT

.PHONY: clean
clean:
	latexmk -f -c ${TARGET}.tex
	rm -rf cache ${TARGET}.tex ${TARGET}.bbl

.PHONY: deep-clean
deep-clean:
	latexmk -f -c ${TARGET}.tex
	rm -rf cache ${TARGET}.tex ${TARGET}.bbl ${TARGET}.pdf
