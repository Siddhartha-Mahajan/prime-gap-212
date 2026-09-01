PYTHON ?= python3
PDFLATEX ?= pdflatex
THREADS ?= 4

.PHONY: all paper verify verify-float verify-integer test clean
all: paper

paper: paper.pdf

paper.pdf: paper.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error paper.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error paper.tex
	$(PDFLATEX) -interaction=nonstopmode -halt-on-error paper.tex

verify:
	$(PYTHON) ancillary/run_certificate.py --threads $(THREADS)

verify-float: verify

verify-integer:
	$(PYTHON) ancillary/run_integer_certificate.py --threads $(THREADS)

test: verify verify-integer
	$(PYTHON) ancillary/test_formulas.py

clean:
	rm -f paper.aux paper.log paper.out paper.toc paper.fls paper.fdb_latexmk paper.synctex.gz
	rm -f ancillary/certify ancillary/certify_integer
	rm -f ancillary/verified_output.json ancillary/verified_integer_output.json
	rm -rf ancillary/__pycache__
