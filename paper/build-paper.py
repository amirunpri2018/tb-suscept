#!/usr/bin/env python3

# Convert a Word document to tex for submission to Scientific Reports.
#
# Usage:
#
# python build-paper.py file.docx > file.tex
#
# Resources:
#
# Scientific Reports template: https://www.overleaf.com/latex/templates/template-for-submissions-to-scientific-reports/xyrztqvdccns#.V5Drne3L9z1
# Nature bibliographic style: http://tug.ctan.org/tex-archive/macros/latex/contrib/nature/naturemag.bst

import docx
import os
import sys
import textwrap

doc_class = "\documentclass[fleqn,10pt]{wlscirep}\n"

# List of packages to be loaded in preamble
packages = ["fixltx2e",
            "textcomp", # for tilde: \texttildelow
            "epstopdf", # Convert EPS to PDF
            "filecontents"] # To contain supplement

# LaTeX macro for labeling supplementary tables and figures.
# http://bytesizebio.net/2013/03/11/adding-supplementary-tables-and-figures-in-latex/
label_supp = """
\\newcommand{\\beginsupplement}{%
 \\setcounter{table}{0}
 \\renewcommand{\\thetable}{S\\arabic{table}}%
 \\setcounter{figure}{0}
 \\renewcommand{\\thefigure}{S\\arabic{figure}}%
 }

"""

# Hack to remove supplement from final PDF while retaining
# cross-references.
# https://tex.stackexchange.com/questions/222934/hide-specific-table-keep-cross-references-and-caption-in-listoftables/222948#222948
remove_supp = """
% Uncomment for removing the supplement
%\includeonly{}

"""

# LaTeX to add at start of supplement
start_supp = """
\\clearpage
\\newpage

\\begin{filecontents}{\jobname-supplement}
\\beginsupplement
"""

# LaTeX to add at end of supplement
end_supp = """
\\end{filecontents}

\include{\jobname-supplement}
"""

def convert_style_to_macro(style):
    if style == "Title":
        return "title"
    if style == "Author":
        return "author"
    if style == "Abstract":
        return "abstract"
    return None

def write_title(text):
    return "\\title{%s}\n"%(text)

def write_author(text):
    # Format author and affiliations. Example output below:
    #
    # \author[1,*]{Alice Author}
    # \affil[1]{Affiliation, department, city, postcode, country}
    # \affil[*]{corresponding.author@email.example}
    #
    if text == "":
        return text
    # If it starts with a superscript, it is an affiliation
    elif text[:17] == "\\textsuperscript{":
        num, name = text[17:].split("}")
        return "\\affil[%s]{%s}"%(num, name)
    # Otherwise it is the list of authors with their affiliation
    # numbers following their names.
    else:
        entries = text.split("},")
        result = ""
        for e in entries:
            parts = e.split("\\textsuperscript{")
            author = parts[0].lstrip(" ")
            affiliations = "".join([x.rstrip("}") for x in parts[1:]])
            result = result + "\\author[%s]{%s}\n"%(affiliations,
                                                    author)
        return result

def write_abstract(text):
    wrapped = textwrap.fill(text, break_long_words = False,
                            break_on_hyphens = False)
    return "\\begin{abstract}\n%s\n\\end{abstract}\n"%(wrapped)

def begin_document():
    return "\\begin{document}\n\\flushbottom\n\\maketitle\n\\thispagestyle{empty}\n"

def write_section(text):
    return "\\section*{%s}\n"%(text)

def write_subsection(text):
    return "\\subsection*{%s}\n"%(text)

def write_subsubsection(text):
    return "\\subsubsection*{%s}\n"%(text)

def convert_run(run):
    result = run.text
    if run.font.subscript:
        result = "\\textsubscript{" + result + "}"
    elif run.font.superscript:
        result = "\\textsuperscript{" + result + "}"
    if run.italic:
        result = "\emph{" + result + "}"
    if run.bold:
        result = "\\textbf{" + result + "}"
    if run.underline:
        result = "\\underline{" + result + "}"
    return result

def use_packages(packages):
    # Packages is a list of strings.
    # Add to preamble with \usepackage{}
    result = ""
    for p in packages:
        result = result + "\\usepackage{" + p + "}\n"
    return result

def convert_special_characters(text):
    # Convert special characters to LaTeX format.
    #
    # Currently supports:
    #   %    -> \%
    #   ~    -> \texttidlelow (requires textcomp package)
    #   <    -> \textless \, (the comma adds a space after it)
    #   >    -> \textgreater \,
    #
    # Sources:
    #   http://tex.stackexchange.com/questions/9363/how-does-one-insert-a-backslash-or-a-tilde-into-latex
    #   http://tex.stackexchange.com/questions/10300/include-fewer-than-and-greater-than-inequality-symbols
    # LaTeX macro for creating nice tilde
    # Had difficulty getting this to work properly
    # http://tex.stackexchange.com/a/9372
    #
    result = text.replace("%", "\%")
    result = result.replace("~", "\\texttildelow")
    result = result.replace("<", "\\textless \,")
    result = result.replace(">", "\\textgreater \,")
    result = result.replace("#", "\#")
    result = result.replace("†", "$\dag$")
    result = result.replace("&", "\&")
    # Backslash escape _ only if it is not in a math equation
    if "$" not in result:
        result = result.replace("_", "\_")
    return result

if __name__ == "__main__":
    fname = sys.argv[1]
    assert fname[-4:] == "docx", "Input file is Word document"
    assert os.path.exists(fname), "Input file exists"
    # Add preamble
    sys.stdout.write(doc_class)
    sys.stdout.write(use_packages(packages))
    sys.stdout.write(label_supp)
    sys.stdout.write(remove_supp)
    # Process input Word document
    d = docx.Document(fname)
    for line in d.paragraphs:
        out = ""
        for run in line.runs:
            out = out + convert_run(run)
        style = line.style.name
        if style != "LaTeX":
            out = convert_special_characters(out)
        if style == "Title":
            out = write_title(out)
        elif style == "Author":
            out = write_author(out)
        elif style == "Abstract":
            out = write_abstract(out) + begin_document()
        elif style == "Heading 1":
            out = write_section(out)
            if "Supplementary information" in out:
                out = start_supp + out
#                out = "\\clearpage\\newpage\n\begin{filecontents}{\jobname-support}\\beginsupplement\n" + out
        elif style == "Heading 2":
            out = write_subsection(out)
            if "Supplementary data" in out:
                out = "\\clearpage\\newpage\n" + out
        elif style == "Heading 3":
            out = write_subsubsection(out)
        elif style == "Normal":
            out = textwrap.fill(out, break_long_words = False,
                                break_on_hyphens = False)
        sys.stdout.write(out + "\n")

    sys.stdout.write(end_supp)
    sys.stdout.write("\\end{document}\n")
