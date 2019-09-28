#!/bin/bash

lualatex -shell-escape ImplementingDFT.tex
makeindex ImplementingDFT
lualatex -shell-escape ImplementingDFT.tex
