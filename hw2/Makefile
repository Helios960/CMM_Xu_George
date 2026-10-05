CXX = g++ # Replace with your preferred compiler or keep g++
CXXFLAGS = -O3 -Wall -Wextra -std=c++17 
# I wanted to do c++20 but VSC somehow doesn't let me and I dont want to debug with the makefile so c++17 it is.
# I am using 03 to optimize for execution speed.
# I am also using Wall and Wextra to see if I have any issues with the code

LATEX = xelatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error

DOC = HW2
TEX_SRC = $(DOC).tex preamble.tex
PLOTS = carbon_activity.pdf carbon_part_c.pdf golf_trajectories.pdf
CSVS = carbon_dt10.csv carbon_dt100.csv carbon_dt1000.csv golf_theta45.csv golf_theta30.csv golf_theta15.csv golf_theta9.csv

.PHONY: all carbon golf plots report clean distclean

all: carbon golf report

# Compiles the binary if needed using the above compiler and flags, then executes with or without WIDTH
carbon: carbon_bin
	./carbon $(if $(WIDTH),--plot=$(WIDTH),)

carbon_bin: carbon.cpp
	$(CXX) $(CXXFLAGS) -o carbon carbon.cpp

# Compiles the binary if needed, then executes with or without THETA
golf: golf_bin
	./golf $(if $(THETA),--plot=$(THETA),)

golf_bin: golf.cpp
	$(CXX) $(CXXFLAGS) -o golf golf.cpp

# Binaries write CSVs directly to disk via ofstream
$(CSVS): carbon_bin golf_bin
	./carbon
	./golf

# Python script processes the generated CSVs into PDF figures
plots: $(PLOTS)

$(PLOTS): plot_csv.py $(CSVS)
	python3 plot_csv.py

# Compiles the PDF report using XeLaTeX (two passes to resolve cross-references)
report: $(DOC).pdf

$(DOC).pdf: $(TEX_SRC) $(PLOTS)
	$(LATEX) $(LATEXFLAGS) $(DOC).tex
	$(LATEX) $(LATEXFLAGS) $(DOC).tex

clean:
	rm -f carbon golf *.o *.aux *.log *.out *.toc *.synctex.gz *.fls *.fdb_latexmk $(CSVS)

distclean: clean
	rm -f $(DOC).pdf $(PLOTS)