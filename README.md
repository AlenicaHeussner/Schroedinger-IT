# Schroedinger-IT: Schrödinger Equation Solver

This repository provides tools to solve the Schrödinger equation for various scenarios. It can load scenarios from input files, solve the Schrödinger equation for those scenarios, save the results, and visualize the potential, eigenstates, and expected values.

## Modules

### main.py
Responsible for reading the input files, solving the Schrödinger equation for various scenarios, and saving the results.

**Usage:**
\```python
python main.py --dir . --input schrodinger.inp
\```

### solver.py
Provides functionalities to load scenarios from an input file and solve the Schrödinger equation across different scenarios.

### solver_test.py
Contains functions to test the solver's ability to load scenarios from a file and its accuracy in solving the Schrödinger equation for given scenarios.

### visualize.py
Provides functionalities to visualize the potential, eigenstates, and other related data resulting from the solution of the Schrödinger equation.

**Usage for Different Scenarios:**

- **Scenario 1:**
\```python
python visualize.py --dir ./scenario_1 --xlim_left -2 2 --xlim_right 0 1 --ylim 0 4 --scale_factor 0.25 --xtol-start-left -0.2 --xtol-end-left -0.2 --xtol-start-right 0 --xtol-end-right -0.2 --ytol-start -0.2 --ytol-end -0.2
\```
- **Scenario 2:**
\```python
python visualize.py --dir ./scenario_2 --xlim_left -2 2 --xlim_right 0 0.8 --ylim -10 0 --scale_factor 1 --xtol-start-left -0.2 --xtol-end-left -0.2 --xtol-start-right 0 --xtol-end-right -0.07 --ytol-start -0.5 --ytol-end -1.5
\```
- **Scenario 3:**
\```python
python visualize.py --dir ./scenario_3 --xlim_left -5 5 --xlim_right 0 1.5 --ylim 0 2.5 --scale_factor 0.25 --xtol-start-left -0.5 --xtol-end-left -0.5 --xtol-start-right 0 --xtol-end-right -0.07 --ytol-start -0.1 --ytol-end -0.05
\```
- **Scenario 4:**
\```python
python visualize.py --dir ./scenario_4 --xlim_left -20 20 --xlim_right 0 7 --ylim -1.5 2.5 --scale_factor 0.25 --xtol-start-left -2 --xtol-end-left -2 --xtol-start-right 0 --xtol-end-right -0.6 --ytol-start -0.2 --ytol-end -0.1
\```
- **Scenario 5:**
\```python
python visualize.py --dir ./scenario_5 --xlim_left -20 20 --xlim_right 0 7 --ylim -0.5 0.6 --scale_factor 0.16 --xtol-start-left -2 --xtol-end-left -2 --xtol-start-right 0 --xtol-end-right -1.8 --ytol-start -0.13 --ytol-end -0.05
\```
- **Scenario 6:**
\```python
python visualize.py --dir ./scenario_6 --xlim_left 0 20 --xlim_right 0 3 --ylim 0 3.5 --scale_factor 0.27 --xtol-start-left -1 --xtol-end-left -1 --xtol-start-right 0 --xtol-end-right -0.87 --ytol-start -0.2 --ytol-end -0.375
\```

Optional flags include settings for saving plots, scaling wavefunctions, and adjusting visualization settings.

## Requirements
- Python 3.x
- NumPy
- SciPy
- Matplotlib

## Installation

1. Clone this repository:
\```bash
git clone https://github.com/Leni128/Schroedinger-IT
\```

2. Navigate to the repository directory:
\```bash
cd Schroedinger-IT
\```

3. Install the required packages:
\```bash
pip install numpy scipy matplotlib
\```

## Usage
To solve the Schrödinger equation for given scenarios and save the results:
\```python
python main.py --dir . --input schrodinger.inp
\```

To visualize the results:
\```python
python visualize.py --dir ./scenario_x --xlim_left --xlim_right --ylim --scale_factor --xtol-start-left --xtol-end-left --xtol-start-right --xtol-end-right --ytol-start --ytol-end
\```
python visualize.py --dir ./scenario_x --xlim_left --xlim_right --ylim --scale_factor --xtol-start-left --xtol-end-left --xtol-start-right --xtol-end-right --ytol-start --ytol-end --save_pdf
\```

## Testing
To run tests for the solver:
\```python
python3 -m pytest
\```
