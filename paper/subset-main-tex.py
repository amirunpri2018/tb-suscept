#!/usr/bin/env python3

# Scientific Reports wants the article tex file to only include the
# main text.

import os
import sys

if __name__ == "__main__":
    fname = sys.argv[1]
    assert fname[-3:] == "tex", "Input file is tex document"
    assert os.path.exists(fname), "Input file exists"
    handle = open(fname, "r")
    for line in handle:
        sys.stdout.write(line)
        if "bibliography{references}" in line:
            sys.stdout.write("\n\end{document}\n")
            break

    handle.close()
