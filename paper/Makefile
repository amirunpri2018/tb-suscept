
NAME := tb-suscept
BIB := references

all: $(NAME).pdf

$(NAME).pdf: $(NAME).tex $(BIB).bib
	python format-bibtex.py $(BIB).bib
	pdflatex $(NAME)
	bibtex $(NAME)
	pdflatex $(NAME)
	pdflatex $(NAME)

$(NAME).tex: $(NAME).docx
	python build-paper.py $< > $@

clean:
	rm -f $(NAME).aux $(NAME).bbl $(NAME).blg $(NAME).dvi $(NAME).log $(NAME).out $(NAME).pdf
