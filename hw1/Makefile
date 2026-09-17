PYTHON ?= /opt/anaconda3/bin/python
FXN    ?= cos,sin,sinc
TXT    ?= data.txt
FMT    ?= 

PRINT_FLAG = $(if $(FMT),--print=$(FMT))

.PHONY: plot write read

plot:
	$(PYTHON) trigonometry.py --function=$(FXN) $(PRINT_FLAG)

write:
	$(PYTHON) trigonometry.py --function=$(FXN) --write=$(TXT) $(PRINT_FLAG)

read:
	$(PYTHON) trigonometry.py --read_from_file=$(TXT) $(PRINT_FLAG)