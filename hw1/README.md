# Homework 1: Programming and Plotting Basics

## Requirements to Grade
- Python 3 with `numpy` and `matplotlib`
- TeX Live with `xelatex` (required for OpenType font and Unicode math support, should come with a standard MacTex install. Came with mine at least.)

## Build Instructions For All Documents
- Generating Plot: `make plot`
- Exporting ASCII table: `make write TXT=data.txt`
- Read and plot from file: `make read TXT=data.txt`
- Rebuild report: `make report` (or run `xelatex report.tex` directly after generating the plot)