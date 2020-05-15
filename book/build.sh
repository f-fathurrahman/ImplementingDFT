#!/bin/bash

lualatex -shell-escape ImplementingDFT.tex
bibtex ImplementingDFT
makeindex ImplementingDFT
lualatex -shell-escape ImplementingDFT.tex
