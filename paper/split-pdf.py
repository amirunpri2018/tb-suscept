#!/usr/bin/env python3

import os
import sys
import pdb
import PyPDF2

if __name__ == "__main__":
    # Obtain filename
    fname = sys.argv[1]
    assert fname[-3:] == "pdf", "Input file is PDF document"
    assert os.path.exists(fname), "Input file exists"
    fpath = os.path.dirname(fname)
    if fpath == "":
        fpath = "./"
    else:
        fpath = fpath + "/"

    # Read PDF
    handle_full = open(fname, "rb")
    handle_main = open(fpath + "blischak-et-al.pdf", "wb")
    handle_supp = open(fpath + "blischak-et-al-supplement.pdf", "wb")
    pdf_full = PyPDF2.PdfFileReader(handle_full)
    pdf_main = PyPDF2.PdfFileWriter()
    pdf_supp = PyPDF2.PdfFileWriter()
    main = True
    for num in range(pdf_full.numPages):
        page = pdf_full.getPage(num)
        text = page.extractText()
        if "Supplementaryinformation" in text:
            main = False
        if main:
            pdf_main.addPage(page)
        else:
            pdf_supp.addPage(page)

    # Write main and supplement PDFs
    pdf_main.write(handle_main)
    pdf_supp.write(handle_supp)

    # Close file connections
    handle_full.close()
    handle_main.close()
    handle_supp.close()
