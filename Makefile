all: paper.pdf

paper.aux: paper.tex
	pdflatex paper.tex

paper.bbl: paper.aux paper.bib
	bibtex paper
	pdflatex paper.tex

paper.pdf: paper.bbl
	pdflatex paper.tex

paper.ps: paper.dvi
	dvips paper

paper.dvi: paper.tex paper.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	latex paper.tex

.PHONY: spellcheck
spellcheck: aspell.conf
	aspell --conf ./aspell.conf --check paper.tex

clean:
	rm -f *.pdf
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f *.eps *.[1-9]
	rm -f src/*.mpx *.mpx

mrproper: clean
	rm -f *.ps *.pdf


plot_data/data-scaling.csv:
	python3 src/collect_data.py process-files 'scaling/data/chr21_n10*.ts' $@

# TODO make some substitution rules for this later
figures/data-scaling.pdf: plot_data/data-scaling.csv src/plot.py
	python3 src/plot.py data-scaling $< $@
