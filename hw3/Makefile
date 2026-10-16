.PHONY: all clean report oscillator poisson

REPORT_NAME = HW3

# Diagram dependencies
DIAGRAMS = resonance_curve.pdf nonlinear_comparison.pdf equipotential.pdf \
           dipole_far_field.pdf niter_vs_tol.pdf grid_scaling.pdf \
           jacobi_vs_sor.pdf sensitivity.pdf energies.pdf boundary_cube.pdf

all: report

# Primary compilation target
report: $(REPORT_NAME).pdf

$(REPORT_NAME).pdf: $(REPORT_NAME).tex $(DIAGRAMS)
	xelatex -interaction=nonstopmode $(REPORT_NAME).tex
	xelatex -interaction=nonstopmode $(REPORT_NAME).tex

# Python execution rules to generate the diagrams
resonance_curve.pdf: oscillator.py
	python oscillator.py --part 1b

nonlinear_comparison.pdf: oscillator.py
	python oscillator.py --part 1c

equipotential.pdf dipole_far_field.pdf: poisson.py
	python poisson.py --part 2a

niter_vs_tol.pdf: poisson.py
	python poisson.py --part 2b

grid_scaling.pdf: poisson.py
	python poisson.py --part 2c

jacobi_vs_sor.pdf: poisson.py
	python poisson.py --part 2e

sensitivity.pdf: oscillator.py
	python oscillator.py --part bonus1

energies.pdf: oscillator.py
	python oscillator.py --part bonus2

boundary_cube.pdf: poisson.py
	python poisson.py --part bonus3

# Manual execution targets for assignment requirements
oscillator:
	python oscillator.py --part $(PART)

poisson:
	python poisson.py --part $(PART)

clean:
	rm -f *.aux *.log *.out $(DIAGRAMS)