## paper

The paper was written in a Word document (tb-suscept.docx), converted
to LaTeX (tb-suscept.tex) using the Python script build-paper.py, and
rendered to PDF (tb-suscept.pdf) using pdflatex. For submission it was
further split into a main PDF file (blischak-et-al.pdf), supplement
PDF file (blischak-et-al-supplement.pdf), and main tex file
(blischak-et-al.tex) using split-pdf.py and subset-main-tex.py.

## Build

Run `make` from within `paper/` to build the paper.

## Requirements

To build the paper, you'll need to install GNU Make, LaTeX, and Python3
(including the Python packages Python-Docx and PyPDF2).
