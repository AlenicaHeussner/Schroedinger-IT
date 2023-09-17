"""
Visualization module for the results of solving the Schrödinger equation.

This module provides functionalities to visualize the potential, eigenstates, 
and other related data resulting from the solution of the Schrödinger equation.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def load_numerical_data_from_file(file_path):
    """
    Load numerical data from a given file.
    
    Parameters:
    - file_path (str): The path to the file to load.
    
    Returns:
    - ndarray: A numpy array containing the loaded data.
    """
    return np.loadtxt(file_path)

def visualize_schrodinger_results(
    data_directory,
    save_as_pdf,
    wavefunc_scale,
    x_limits_left_graph,
    x_limits_right_graph,
    y_limits
):

    """
    Visualize the potential, eigenstates, and expected values.
    
    The function plots the potential, eigenstates, and expected values based on data
    loaded from the given directory. The visualization can be adjusted using provided arguments.
    
    Parameters:
    - data_directory (str): Directory containing the data files.
    - save_as_pdf (bool): If True, save the plot as a PDF, otherwise display the plot.
    - wavefunc_scale (float): Factor to scale the wavefunctions for visualization.
    - x_limits_left_graph (tuple): Limits for the x-axis of the left graph.
    - x_limits_right_graph (tuple): Limits for the x-axis of the right graph.
    - y_limits (tuple): Limits for the y-axis.
    """
    file_path = f"{data_directory}/potential.dat"
    loaded_data = load_numerical_data_from_file(file_path)
    position_values, potential_values = loaded_data.T

    energy_values = load_numerical_data_from_file(f"{data_directory}/energies.dat")
    wave_functions = load_numerical_data_from_file(f"{data_directory}/wavefuncs.dat")[:, 1:]
    expected_values = load_numerical_data_from_file(f"{data_directory}/expvalues.dat")

    _, (left_graph, right_graph) = plt.subplots(1, 2, figsize=(15, 6))
    plt.subplots_adjust(left=0.3, bottom=0.238, right=0.611, top=0.82, wspace=0.202, hspace=0.192)

    left_graph.plot(position_values, potential_values, label="Potential", color="black")
    colors = ["blue", "red"]  # Alternating colors for wavefunctions
    for index, energy in enumerate(energy_values):
        left_graph.axhline(y=energy, color='lightgray', linestyle='-', linewidth=1.5)
        y_values_for_wavefunc = wave_functions[:, index] * wavefunc_scale + energy
        label_for_wavefunc = f"Eigenstate {index+1}"
        wavefunc_color = colors[index % 2]
        left_graph.plot(position_values, y_values_for_wavefunc,
                 label=label_for_wavefunc,
                 color=wavefunc_color)

        left_graph.plot(expected_values[index, 0], energy, 'gx', markersize=8, markeredgewidth=1)
    left_graph.set_title("Potential, Eigenstates, (x)")
    left_graph.set_xlabel("x [Bohr]")
    left_graph.set_ylabel("Energy [Hartree]")
    if x_limits_left_graph:
        left_graph.set_xlim(x_limits_left_graph)
    if y_limits:
        left_graph.set_ylim(y_limits)

    for index, energy in enumerate(energy_values):
        right_graph.axhline(y=energy, color='lightgray', linestyle='-', linewidth=1.5)
    for index, sigma_value in enumerate(expected_values[:, 1]):
        sigma_marker_color = '#FF00FF'
        sigma_marker_style = '+'
        sigma_marker_size = 12
        sigma_marker_edge_width = 1.2

        right_graph.plot(sigma_value, energy_values[index],
                 color=sigma_marker_color,
                 marker=sigma_marker_style,
                 markersize=sigma_marker_size,
                 markeredgewidth=sigma_marker_edge_width)

    right_graph.set_title("σₓ")
    right_graph.set_xlabel("[Bohr]")
    right_graph.set_ylim(left_graph.get_ylim())
    right_graph.yaxis.set_ticklabels([])
    right_graph.yaxis.set_ticks([])
    if x_limits_right_graph:
        right_graph.set_xlim(x_limits_right_graph)

    if save_as_pdf:
        plt.savefig(f"{data_directory}/revised_visualization.pdf")
    else:
        plt.show()

def command_line_interface():
    """
    Command-line interface (CLI) for the visualization module.
    
    Parses the command-line arguments and calls the visualization function.
    """
    parser_description = "Visualize the results from solving the Schrödinger equation."
    parser = argparse.ArgumentParser(description=parser_description)

    parser.add_argument("--data_dir",
                        default=".",
                        help="Directory where the results are saved.")
    parser.add_argument("--save_as_pdf",
                        action="store_true",
                        help="Save the visualization as a PDF instead of displaying it.")
    parser.add_argument("--wavefunc_scale",
                        type=float,
                        default=1.0,
                        help="Factor to scale the wavefunctions for better visualization.")
    parser.add_argument("--x_limits_left_graph",
                        nargs=2,
                        type=float,
                        default=[None, None],
                        help="Limits for the x-axis of the left graph.")
    parser.add_argument("--x_limits_right_graph",
                        nargs=2,
                        type=float,
                        default=[None, None],
                        help="Limits for the x-axis of the right graph.")
    parser.add_argument("--y_limits",
                        nargs=2,
                        type=float,
                        default=[None, None],
                        help="Limits for the y-axis.")
    # I've kept the tolerance arguments unchanged for brevity
    parser.add_argument('--ytol-start',
                        type=float,
                        default=None,
                        help='Tolerance for the start of the Y-axis.')
    parser.add_argument('--ytol-end',
                        type=float,
                        default=None,
                        help='Tolerance for the end of the Y-axis.')
    parser.add_argument('--xtol-start-left',
                        type=float,
                        default=None,
                        help='Tolerance for the start of the X-axis for the left plot.')
    parser.add_argument('--xtol-end-left',
                        type=float,
                        default=None,
                        help='Tolerance for the end of the X-axis for the left plot.')
    parser.add_argument('--xtol-start-right',
                        type=float,
                        default=None,
                        help='Tolerance for the start of the X-axis for the right plot.')
    parser.add_argument('--xtol-end-right',
                        type=float,
                        default=None,
                        help='Tolerance for the end of the X-axis for the right plot.')

    args = parser.parse_args()

    # Applying the tolerances
    x_limits_left = list(args.x_limits_left_graph)
    x_limits_right = list(args.x_limits_right_graph)
    y_limits = list(args.y_limits)

    if x_limits_left[0] is not None and args.xtol_start_left is not None:
        x_limits_left[0] += args.xtol_start_left
    if x_limits_left[1] is not None and args.xtol_end_left is not None:
        x_limits_left[1] -= args.xtol_end_left
    if x_limits_right[0] is not None and args.xtol_start_right is not None:
        x_limits_right[0] += args.xtol_start_right
    if x_limits_right[1] is not None and args.xtol_end_right is not None:
        x_limits_right[1] -= args.xtol_end_right
    if y_limits[0] is not None and args.ytol_start is not None:
        y_limits[0] += args.ytol_start
    if y_limits[1] is not None and args.ytol_end is not None:
        y_limits[1] -= args.ytol_end

    # Call the visualize_schrodinger_results function
    visualize_schrodinger_results(
        args.data_dir,
        args.save_as_pdf,
        args.wavefunc_scale,
        tuple(x_limits_left),
        tuple(x_limits_right),
        tuple(y_limits)
    )

if __name__ == "__main__":
    command_line_interface()
