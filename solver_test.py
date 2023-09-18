import numpy as np
from solver import load_scenarios_from_file, solve_schroedinger_equation
from numpy.testing import assert_allclose

def test_load_scenarios():
    scenarios = load_scenarios_from_file("schrodinger.inp")
    assert len(scenarios) == 6  # basierend auf Ihrer schrodinger.inp Datei

def test_scenario():
    
    scenarios = load_scenarios_from_file('./schrodinger.inp')

    print("test scenario data...")
    for i in range(1, 7):
        with open(f'./scenario_{i}/energies.dat', 'r', encoding='utf-8') as file:
            energies_expected = [float(line.strip()) for line in file.readlines()]

        with open(f'./scenario_{i}/expvalues.dat', 'r', encoding='utf-8') as file:
            expectedValues = [(float(part[0]), float(part[1])) for part in (line.strip().split() for line in file)]

        wavefunc_expected = []
        with open(f'./scenario_{i}/wavefuncs.dat', 'r', encoding='utf-8') as file:
            for line in file:
                values = [float(value) for value in line.strip().split()]
                wavefunc_expected.append(np.array(values))
        
        with open(f'./scenario_{i}/potential.dat', 'r', encoding='utf-8') as file:
            potentials_per_position_exp = [(float(part[0]), float(part[1])) for part in (line.strip().split() for line in file)]

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
        assert_allclose(expectedValues, expected_values, atol=tolerance)
        assert_allclose(potentials_per_position_exp, potentials_per_position_calc, atol=tolerance)
        assert_allclose(wavefunc_expected, wavefunctions_data, atol=tolerance)
        
    pass
