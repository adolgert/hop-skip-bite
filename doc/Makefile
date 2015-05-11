PYTHON=python3

hopskip.pdf: hopskip.tex excess.pdf
	pdflatex hopskip.tex
	pdflatex hopskip.tex

excess.pdf: excess_plot.py
	$(PYTHON) excess_plot.py

