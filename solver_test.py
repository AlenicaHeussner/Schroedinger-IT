"""
Testing module for the Schrodinger equation solver.

This module contains functions to test the solver's ability to load scenarios from a file
and its accuracy in solving the Schrodinger equation for given scenarios.
"""

import numpy as np
from numpy.testing import assert_allclose
from solver import load_scenarios_from_file, solve_schroedinger_equation

def test_load_scenarios():
    """
    Test the function that loads scenarios from a file.
    
    This function asserts that the loaded scenarios from the file "schrodinger.inp"
    have a length of 6.
    """
    scenarios = load_scenarios_from_file("schrodinger.inp")
    assert len(scenarios) == 6  # based on schrodinger.inp file

def test_scenario():
    """
    Test a specific scenario against expected data.
    
    This function performs the following:
    - Loads scenarios from a file.
    - For each scenario, reads:
    - expected energies,
    - expected values,
    - wave functions,
    - potentials from corresponding files.
    - Solves the Schrodinger equation for the scenario.
    - Asserts that the calculated values are close to the expected values with a certain tolerance.
    """
    scenarios = load_scenarios_from_file('./schrodinger.inp')

    print("test scenario data...")
    for i in range(1, 7):
        with open(f'./scenario_{i}/energies.dat', 'r', encoding='utf-8') as file:
            energies_expected = [float(line.strip()) for line in file.readlines()]

        with open(f'./scenario_{i}/expvalues.dat', 'r', encoding='utf-8') as file:
            expected_values = []
            for line in file:
                stripped_line = line.strip()
                parts = stripped_line.split()
                value_tuple = (float(parts[0]), float(parts[1]))
                expected_values.append(value_tuple)


        wavefunc_expected = []
        with open(f'./scenario_{i}/wavefuncs.dat', 'r', encoding='utf-8') as file:
            for line in file:
                values = [float(value) for value in line.strip().split()]
                wavefunc_expected.append(np.array(values))

        with open(f'./scenario_{i}/potential.dat', 'r', encoding='utf-8') as file:
            potentials_per_position_exp = []
            for line in file:
                stripped_line = line.strip()
                parts = stripped_line.split()
                value_tuple = (float(parts[0]), float(parts[1]))
                potentials_per_position_exp.append(value_tuple)


        scenario_data = scenarios[i-1]

        tolerance = 0.0005

        (position_values,
         potential_values,
         expected_values,
         eigenvalues,
         eigenvectors) = solve_schroedinger_equation(scenario_data)

        potentials_per_position_calc = np.column_stack((position_values, potential_values))

        eigenvectors_data = [
            eigenvectors[:, i]
            for i in range(eigenvectors.shape[1])
        ]

        wavefunctions_data = np.column_stack([position_values] + eigenvectors_data)

        assert_allclose(energies_expected, eigenvalues, atol=tolerance)
        assert_allclose(expected_values, expected_values, atol=tolerance)
        assert_allclose(potentials_per_position_exp, potentials_per_position_calc, atol=tolerance)
        assert_allclose(wavefunc_expected, wavefunctions_data, atol=tolerance)
