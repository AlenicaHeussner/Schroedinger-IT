"""
Modul zum Lösen der Schrödinger-Gleichung für verschiedene Szenarien.
"""

import os
import argparse
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.linalg import eigh_tridiagonal

def load_input_data(filename):
    """Lädt Eingabedaten aus einer Datei und gibt eine Liste von Szenarien zurück."""
    scenarios = []

    with open(filename, 'r', encoding="utf-8") as file:
        while True:
            description = None
            while not description:
                line = file.readline()
                if not line:
                    return scenarios
                description = line.strip().split("#")[0]

            mass_line = file.readline().strip().split("#")[0]
            mass = float(mass_line)

            x_params_line = file.readline().strip().split("#")[0]
            x_min, x_max, n_points = map(float, x_params_line.split())

            eigenvalues_line = file.readline().strip().split("#")[0]
            first_eigenvalue, last_eigenvalue = map(int, eigenvalues_line.split())

            interpolation_type_line = file.readline().strip().split("#")[0]
            interpolation_type = interpolation_type_line.lower().strip()

            num_interpolation_points_line = file.readline().strip().split("#")[0]
            num_interpolation_points = int(num_interpolation_points_line)

            interpolation_data = []
            for _ in range(num_interpolation_points):
                data_line = file.readline().strip().split("#")[0]
                data_point = tuple(map(float, data_line.split()))
                interpolation_data.append(data_point)

            scenario = (
                mass,
                x_min,
                x_max,
                n_points,
                first_eigenvalue,
                last_eigenvalue,
                interpolation_type,
                interpolation_data)
            scenarios.append(scenario)

def generate_potential(x_value, interpolation_type, interpolation_data):
    """Generiert ein Potential basierend auf gegebenen Interpolationsdaten."""
    x_data, y_data = zip(*interpolation_data)

    if interpolation_type == 'linear':
        interpolator = interp1d(x_data, y_data, kind='linear', fill_value="extrapolate")
    elif interpolation_type == 'cspline':
        interpolator = CubicSpline(x_data, y_data, extrapolate=True)
    elif interpolation_type == 'polynomial':
        interpolator = interp1d(x_data, y_data, kind=len(x_data) - 1, fill_value="extrapolate")
    else:
        raise ValueError(f"Unknown interpolation type: {interpolation_type}")

    return interpolator(x_value)

def main_solver_with_expvalues(directory, input_file):
    """Hauptfunktion zum Lösen der Schrödinger-Gleichung."""
    scenarios = load_input_data(input_file)
    results = []

    for index, scenario in enumerate(scenarios, start=1):
        (
        mass,
        x_min,
        x_max,
        n_points,
        first_eigenvalue,
        last_eigenvalue,
        interpolation_type,
        interpolation_data
        ) = scenario

        x_values = np.linspace(x_min, x_max, int(n_points))
        potential = generate_potential(x_values, interpolation_type, interpolation_data)

        delta_x = x_values[1] - x_values[0]
        main_diag = potential + 1.0 / (mass * delta_x**2)
        off_diag = -0.5 / (mass * delta_x**2) * np.ones_like(x_values[:-1])

        eigenvalues, eigenvectors = eigh_tridiagonal(main_diag, off_diag)

        selected_eigenvalues = eigenvalues[first_eigenvalue-1:last_eigenvalue]
        selected_eigenvectors = eigenvectors[:, first_eigenvalue-1:last_eigenvalue]
        for i in range(selected_eigenvectors.shape[1]):
            norm = np.sqrt(np.trapz(selected_eigenvectors[:, i]**2, x_values))
            selected_eigenvectors[:, i] /= norm

        exp_values = []
        for i in range(len(selected_eigenvalues)):
            wave_function = selected_eigenvectors[1:-1, i]
            x_values_interior = x_values[1:-1]
            integrand = wave_function * x_values_interior * wave_function
            expectation_x = np.trapz(integrand, x_values_interior)
            integrand_x2 = wave_function * (x_values_interior**2) * wave_function
            expectation_x2 = np.trapz(integrand_x2, x_values_interior)

            sigma_x = np.sqrt(expectation_x2 - expectation_x**2)
            exp_values.append((expectation_x, sigma_x))

        scenario_dir = os.path.join(directory, f"scenario_{index}")
        os.makedirs(scenario_dir, exist_ok=True)

        # Save the potential values to a file
        file_path = os.path.join(scenario_dir, "potential.dat")
        data_to_save = np.column_stack([x_values, potential])
        np.savetxt(file_path, data_to_save)

        # Save the eigenvalues to a file
        np.savetxt(os.path.join(scenario_dir, "energies.dat"), selected_eigenvalues)

        # Save the eigenvectors to a file
        eigenvectors_list = [
        selected_eigenvectors[:, i]
        for i in range(selected_eigenvectors.shape[1])
        ]
        wavefuncs_data = np.column_stack([x_values] + eigenvectors_list)
        np.savetxt(os.path.join(scenario_dir, "wavefuncs.dat"), wavefuncs_data)

        # Save the exp_values to a file
        expvalues_file_path = os.path.join(scenario_dir, "expvalues.dat")
        np.savetxt(expvalues_file_path, exp_values)

        results.append((x_values, potential, selected_eigenvalues, selected_eigenvectors))

    return results

def main_cli():
    """CLI-Funktion für die Ausführung des Solvers."""
    description_text = "Solve Schrödinger equation for various scenarios."
    parser = argparse.ArgumentParser(description=description_text)
    argument_name = "--dir"
    argument_type = str
    default_value = "."
    help_text = "Directory to save the output. Default is the current directory."

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = "--input"
    argument_type = str
    required_value = True
    help_text = "Input file containing the scenarios."

    parser.add_argument(argument_name,
                    type=argument_type,
                    required=required_value,
                    help=help_text)
    args = parser.parse_args()

    main_solver_with_expvalues(args.dir, args.input)

if __name__ == "__main__":
    main_cli()
