"""
Hauptmodul für das Lösen der Schrödinger-Gleichung für verschiedene Szenarien.
Dieses Modul liest Eingabedateien, löst die Gleichung und speichert die Ergebnisse.
"""

import os
import argparse
import numpy as np
from solver import main_solver_with_expvalues


def save_to_file(filename, data):
    """Speichert die Daten in einer Datei."""
    np.savetxt(filename, data)


def main(args):
    """Hauptfunktion, die Szenarien bearbeitet und Daten speichert."""
    scenarios = main_solver_with_expvalues(args.dir, args.input)
    print(f"Processing {len(scenarios)} scenarios...")
    for idx, (x_vals,
            potential,
            selected_eigenvalues,
            selected_eigenvectors) in enumerate(scenarios, start=1):
        scenario_dir = os.path.join(args.dir, f"scenario_{idx}")
        os.makedirs(scenario_dir, exist_ok=True)

        save_to_file(f"{scenario_dir}/potential.dat", np.column_stack((x_vals, potential)))
        save_to_file(f"{scenario_dir}/energies.dat", selected_eigenvalues)

        # Create a single wavefuncs.dat file for all eigenfunctions
        eigenvectors_list = [
        selected_eigenvectors[:, i]
        for i in range(selected_eigenvectors.shape[1])
        ]
        wavefuncs_data = np.column_stack([x_vals] + eigenvectors_list)
        save_to_file(f"{scenario_dir}/wavefuncs.dat", wavefuncs_data)


def main_cli():
    """CLI-Handler für das Hauptmodul."""
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

    main(args)


if __name__ == "__main__":
    main_cli()
