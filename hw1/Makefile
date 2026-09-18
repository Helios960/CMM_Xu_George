PYTHON ?= python3
FXN    ?= cos,sin,sinc # Only useful for plot
TXT    ?= data.txt # Only useful for read and write
FMT    ?= # Defauls to no plot

PRINT_FLAG = $(if $(FMT),--print=$(FMT)) # this line only passes '--print' on if it exists. Otherwise, it does nothing.

.PHONY: plot write read report

plot:
	$(PYTHON) trigonometry.py --function=$(FXN) $(PRINT_FLAG)

write:
	$(PYTHON) trigonometry.py --function=$(FXN) --write=$(TXT) $(PRINT_FLAG)

read:
	$(PYTHON) trigonometry.py --read_from_file=$(TXT) $(PRINT_FLAG)

report:
	$(PYTHON) trigonometry.py --function=cos,sinc,sin --print=jpeg
	xelatex -interaction=nonstopmode -halt-on-error report.tex
	xelatex -interaction=nonstopmode -halt-on-error report.tex