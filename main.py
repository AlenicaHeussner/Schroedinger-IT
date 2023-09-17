"""
Main module for solving the Schrödinger equation for various scenarios.
This module reads input files, solves the equation, and saves the results.
"""

import os
import argparse
import numpy as np
from solver import load_scenarios_from_file, solve_schroedinger_equation

def save_data_to_file(filename, data_array):
    """Saves the data to a file."""
    np.savetxt(filename, data_array)


def process_and_save_scenarios(arguments):
    """Main function that processes scenarios and saves the data."""
    #all_scenarios = process_scenarios(arguments.directory, arguments.input_file)

    scenarios = load_scenarios_from_file(arguments.input_file)

    print(f"Processing {len(scenarios)} scenarios...")
    for index, scenario_data in enumerate(scenarios, start=1):
        scenario_output_dir = os.path.join(arguments.directory, f"scenario_{index}")
        os.makedirs(scenario_output_dir, exist_ok=True)

        (position_values,
         potential_values,
         expected_values,
         eigenvalues,
         eigenvectors) = solve_schroedinger_equation(scenario_data)

        scenario_directory = os.path.join(arguments.directory, f"scenario_{index}")
        os.makedirs(scenario_directory, exist_ok=True)

        file_path = f"{scenario_directory}/potential.dat"
        stacked_data = np.column_stack((position_values, potential_values))
        save_data_to_file(file_path, stacked_data)

        save_data_to_file(f"{scenario_directory}/energies.dat", eigenvalues)

        # Create a single wavefuncs.dat file for all eigenfunctions
        eigenvectors_data = [
            eigenvectors[:, i]
            for i in range(eigenvectors.shape[1])
        ]
        wavefunctions_data = np.column_stack([position_values] + eigenvectors_data)
        save_data_to_file(f"{scenario_directory}/wavefuncs.dat", wavefunctions_data)

        # Save the expected values to a file
        expvalues_file_path = os.path.join(scenario_directory, "expvalues.dat")
        np.savetxt(expvalues_file_path, expected_values)


def main_command_line_interface():
    """CLI handler for the main module."""
    description = "Solves the Schrödinger equation for various scenarios."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--directory",
                        type=str,
                        default=".",
                        help="Directory to save the output. Default is the current directory.")
    parser.add_argument("--input_file",
                        type=str,
                        required=True,
                        help="Input file containing the scenarios.")
    parsed_args = parser.parse_args()

    process_and_save_scenarios(parsed_args)


if __name__ == "__main__":
    main_command_line_interface()
