"""
Module for solving the Schrödinger equation across different scenarios.
"""

import os
import argparse
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
from scipy.linalg import eigh_tridiagonal

def load_scenarios_from_file(filename):
    """Loads input data from a file and returns a list of scenarios."""
    loaded_scenarios = []

    with open(filename, 'r', encoding="utf-8") as file:
        while True:
            scenario_description = None
            while not scenario_description:
                line = file.readline()
                if not line:
                    return loaded_scenarios
                scenario_description = line.strip().split("#")[0]

            particle_mass = float(file.readline().strip().split("#")[0])
            line_content = file.readline().strip()
            cleaned_content = line_content.split("#")[0]
            split_values = cleaned_content.split()
            x_range_min, x_range_max, num_points = map(float, split_values)

            line_content = file.readline().strip()
            cleaned_content = line_content.split("#")[0]
            split_values = cleaned_content.split()
            starting_eigenvalue, ending_eigenvalue = map(int, split_values)

            interp_method = file.readline().strip().split("#")[0].lower().strip()
            num_interp_points = int(file.readline().strip().split("#")[0])

            interp_data_points = []
            for _ in range(num_interp_points):
                data_point = tuple(map(float, file.readline().strip().split("#")[0].split()))
                interp_data_points.append(data_point)

            scenario_data = (
                particle_mass,
                x_range_min,
                x_range_max,
                num_points,
                starting_eigenvalue,
                ending_eigenvalue,
                interp_method,
                interp_data_points)
            loaded_scenarios.append(scenario_data)

def create_potential_from_data(x_values, interpolation_method, interpolation_points):
    """Generates a potential based on given interpolation data."""
    x_data_points, y_data_points = zip(*interpolation_points)

    if interpolation_method == 'linear':
        interp_kind = 'linear'
        interp_fill_value = "extrapolate"
        potential_function = interp1d(x_data_points,
                                     y_data_points,
                                     kind=interp_kind,
                                     fill_value=interp_fill_value)

    elif interpolation_method == 'cspline':
        potential_function = CubicSpline(x_data_points, y_data_points, extrapolate=True)
    elif interpolation_method == 'polynomial':
        interp_kind = len(x_data_points) - 1
        interp_fill_value = "extrapolate"

        potential_function = interp1d(
            x_data_points,
            y_data_points,
            kind=interp_kind,
            fill_value=interp_fill_value
        )

    else:
        raise ValueError(f"Unknown interpolation type: {interpolation_method}")

    return potential_function(x_values)

def solve_schrodinger_equation(output_directory, input_filename):
    """Main function for solving the Schrödinger equation."""
    scenarios = load_scenarios_from_file(input_filename)
    computed_results = []

    for index, scenario_data in enumerate(scenarios, start=1):
        (
        particle_mass,
        x_min_range,
        x_max_range,
        num_of_points,
        start_eigenvalue,
        end_eigenvalue,
        interp_method,
        interp_data_points
        ) = scenario_data

        position_values = np.linspace(x_min_range, x_max_range, int(num_of_points))
        potential_values = create_potential_from_data(
            position_values,
            interp_method,
            interp_data_points
        )


        delta_x = position_values[1] - position_values[0]
        main_diag = potential_values + 1.0 / (particle_mass * delta_x**2)
        off_diag_values = -0.5 / (particle_mass * delta_x**2) * np.ones_like(position_values[:-1])

        computed_eigenvalues, computed_eigenvectors = eigh_tridiagonal(main_diag, off_diag_values)

        selected_eigenvalues = computed_eigenvalues[start_eigenvalue-1:end_eigenvalue]
        selected_eigenvectors = computed_eigenvectors[:, start_eigenvalue-1:end_eigenvalue]
        for i in range(selected_eigenvectors.shape[1]):
            norm_factor = np.sqrt(np.trapz(selected_eigenvectors[:, i]**2, position_values))
            selected_eigenvectors[:, i] /= norm_factor

        expected_values = []
        for i in range(len(selected_eigenvalues)):
            eigen_function = selected_eigenvectors[1:-1, i]
            interior_position_values = position_values[1:-1]
            integrand_for_exp_x = eigen_function * interior_position_values * eigen_function
            exp_x = np.trapz(integrand_for_exp_x, interior_position_values)
            integrand_for_exp_x2 = eigen_function * (interior_position_values**2) * eigen_function
            exp_x2 = np.trapz(integrand_for_exp_x2, interior_position_values)

            sigma_x = np.sqrt(exp_x2 - exp_x**2)
            expected_values.append((exp_x, sigma_x))

        scenario_output_dir = os.path.join(output_directory, f"scenario_{index}")
        os.makedirs(scenario_output_dir, exist_ok=True)

        # Save the potential values to a file
        potential_file_path = os.path.join(scenario_output_dir, "potential.dat")
        potential_data_to_save = np.column_stack([position_values, potential_values])
        np.savetxt(potential_file_path, potential_data_to_save)

        # Save the eigenvalues to a file
        np.savetxt(os.path.join(scenario_output_dir, "energies.dat"), selected_eigenvalues)

        # Save the eigenvectors to a file
        eigenvectors_data_list = [
        selected_eigenvectors[:, i]
        for i in range(selected_eigenvectors.shape[1])
        ]
        wavefunctions_combined_data = np.column_stack([position_values] + eigenvectors_data_list)
        np.savetxt(os.path.join(scenario_output_dir, "wavefuncs.dat"), wavefunctions_combined_data)

        # Save the expected values to a file
        expvalues_file_path = os.path.join(scenario_output_dir, "expvalues.dat")
        np.savetxt(expvalues_file_path, expected_values)

        computed_results.append(
            (
                position_values,
                potential_values,
                selected_eigenvalues,
                selected_eigenvectors
            )
        )


    return computed_results

def command_line_interface():
    """CLI function for executing the solver."""
    description_text = "Solve Schrödinger equation for various scenarios."
    parser = argparse.ArgumentParser(description=description_text)
    parser.add_argument("--directory",
                        type=str,
                        default=".",
                        help="Directory to save the output. Default is the current directory.")
    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="Input file containing the scenarios.")
    parsed_arguments = parser.parse_args()

    solve_schrodinger_equation(parsed_arguments.directory, parsed_arguments.input)

if __name__ == "__main__":
    command_line_interface()
