all: paper.pdf

paper.aux: paper.tex
	pdflatex -shell-escape paper.tex

paper.bbl: paper.aux paper.bib
	bibtex paper
	pdflatex -shell-escape paper.tex

paper.pdf: paper.bbl
	pdflatex -shell-escape paper.tex

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
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -fR _minted*

mrproper: clean
	rm -f *.ps *.pdf



# TODO these rules for creating figures and figure data are out of date.
# The plan is to systematise things once the basic structure has settled
# down a bit more

plot_data/data-scaling.csv:
	python3 src/collect_data.py file-size 'scaling/data/chr21_n10*.ts' $@

# TODO make rule for time-scaling

# TODO make some substitution rules for this later
figures/data-scaling.pdf: plot_data/data-scaling.csv src/plot.py
	python3 src/plot.py data-scaling plot_data/data-scaling.csv  \
		plot_data/time-scaling.csv figures/data-scaling.pdf
