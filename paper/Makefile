
NAME := tb-suscept
BIB := references
FINAL := blischak-et-al

all: $(FINAL).pdf $(FINAL)-supplement.pdf $(FINAL).tex

$(FINAL).pdf: $(FINAL).tex $(FINAL).aux $(FINAL)-supplement.aux $(FINAL)-supplement.pdf
	bibtex $(FINAL)
	pdflatex -shell-escape $(FINAL)

$(FINAL)-supplement.aux: $(NAME)-supplement.aux
	sed "s/$(NAME)/$(FINAL)/" $< > $@

$(FINAL).aux: $(NAME).aux
	sed "s/$(NAME)/$(FINAL)/" $< > $@

$(FINAL).tex: $(NAME).tex
	sed s/'%\\includeonly'/'\\includeonly'/ $< > $@

$(FINAL)-supplement.pdf: $(NAME).pdf
	python split-pdf.py $<
	rm $(FINAL).pdf

%.pdf %.aux %-supplement.aux: %.tex $(BIB).bib
	python format-bibtex.py $(BIB).bib
	pdflatex -shell-escape $(NAME)
	bibtex $(NAME)
	pdflatex -shell-escape $(NAME)
	pdflatex -shell-escape $(NAME)

$(NAME).tex: $(NAME).docx
	python build-paper.py $< > $@

clean:
	rm -f $(NAME)*.aux $(NAME).bbl $(NAME).blg $(NAME).dvi \
              $(NAME).log $(NAME).out $(NAME).pdf \
              $(NAME)-supplement.tex $(FINAL).tex \
              $(FINAL)*.aux $(FINAL).bbl $(FINAL).blg $(FINAL).dvi \
              $(FINAL).log $(FINAL).out $(FINAL)*.pdf \
              $(FINAL)-supplement.tex \
              *converted-to.pdf ../figure/*converted-to.pdf
