all: doc.pdf slides.pdf

doc.pdf: doc.aux
	pdflatex doc

doc.aux: doc.tex doc.bib
	pdflatex doc
	bibtex doc
	pdflatex doc

slides.pdf: slides.tex
	pdflatex slides
	pdflatex slides

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.nav *.snm *.toc *.vrb *.synctex.gz *.pyc *~

