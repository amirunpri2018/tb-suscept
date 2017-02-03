
NAME := tb-suscept
BIB := references
FINAL := blischak-et-al

all: $(FINAL).pdf $(FINAL)-supplement.pdf

$(FINAL).pdf $(FINAL)-supplement.pdf: $(NAME).pdf
	python split-pdf.py $<

$(NAME).pdf: $(NAME).tex $(BIB).bib
	python format-bibtex.py $(BIB).bib
	pdflatex -shell-escape $(NAME)
	bibtex $(NAME)
	pdflatex -shell-escape $(NAME)
	pdflatex -shell-escape $(NAME)

$(NAME).tex: $(NAME).docx
	python build-paper.py $< > $@

clean:
	rm -f $(NAME).aux $(NAME).bbl $(NAME).blg $(NAME).dvi $(NAME).log $(NAME).out $(NAME).pdf $(FINAL)*.pdf *converted-to.pdf ../figure/*converted-to.pdf
