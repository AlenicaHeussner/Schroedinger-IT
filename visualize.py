"""
Visualization module for the results of solving the Schrödinger equation.

This module provides functionalities to visualize the potential, eigenstates, 
and other related data resulting from the solution of the Schrödinger equation.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    """
    Load numerical data from a given file.
    
    Parameters:
    - filename (str): The path to the file to load.
    
    Returns:
    - ndarray: A numpy array containing the loaded data.
    """
    return np.loadtxt(filename)

def visualize_data(directory, save_pdf, scale_factor, xlim_left, xlim_right, ylim):
    """
    Visualize the potential, eigenstates, and expected values.
    
    The function plots the potential, eigenstates, and expected values based on data
    loaded from given directory. The visualization can be adjusted using provided arguments.
    
    Parameters:
    - directory (str): Directory containing the data files.
    - save_pdf (bool): If True, save the plot as a PDF, otherwise display the plot.
    - scale_factor (float): Factor to scale the wavefunctions for visualization.
    - xlim_left (tuple): Limits for the x-axis of the left graph.
    - xlim_right (tuple): Limits for the x-axis of the right graph.
    - ylim (tuple): Limits for the y-axis.
    """
    x_values, potential = load_data(f"{directory}/potential.dat").T
    eigenvalues = load_data(f"{directory}/energies.dat")
    wavefuncs = load_data(f"{directory}/wavefuncs.dat")[:, 1:]
    expvalues = load_data(f"{directory}/expvalues.dat")

    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    plt.subplots_adjust(left=0.3, bottom=0.238, right=0.611, top=0.82, wspace=0.202, hspace=0.192)

    ax1.plot(x_values, potential, label="Potential", color="black")
    colors = ["blue", "red"]  # Alternating colors for wavefunctions
    for i, eigenvalue in enumerate(eigenvalues):
        ax1.axhline(y=eigenvalue, color='lightgray', linestyle='-', linewidth=1.5)
        y_values = wavefuncs[:, i] * scale_factor + eigenvalue
        label_text = f"Eigenstate {i+1}"
        color_value = colors[i % 2]
        ax1.plot(x_values, y_values,
                 label=label_text,
                 color=color_value)

        ax1.plot(expvalues[i, 0], eigenvalue, 'gx', markersize=8, markeredgewidth=1)
    ax1.set_title("Potential, Eigenstates, (x)")
    ax1.set_xlabel("x [Bohr]")
    ax1.set_ylabel("Energy [Hartree]")
    if xlim_left:
        ax1.set_xlim(xlim_left)
    if ylim:
        ax1.set_ylim(ylim)

    for i, eigenvalue in enumerate(eigenvalues):
        ax2.axhline(y=eigenvalue, color='lightgray', linestyle='-', linewidth=1.5)
    for i, sigma in enumerate(expvalues[:, 1]):
        color_value = '#FF00FF'
        marker_style = '+'
        marker_size_value = 12
        marker_edge_width = 1.2

        ax2.plot(sigma, eigenvalues[i],
                 color=color_value,
                 marker=marker_style,
                 markersize=marker_size_value,
                 markeredgewidth=marker_edge_width)

    ax2.set_title("σₓ")
    ax2.set_xlabel("[Bohr]")
    ax2.set_ylim(ax1.get_ylim())
    ax2.yaxis.set_ticklabels([])
    ax2.yaxis.set_ticks([])
    if xlim_right:
        ax2.set_xlim(xlim_right)

    if save_pdf:
        plt.savefig(f"{directory}/revised_visualization.pdf")
    else:
        plt.show()

def main():
    """
    Command-line interface (CLI) for the visualization module.
    
    Parses the command-line arguments and calls the visualization function.
    """
    description_text = "Visualize the results from solving the Schrödinger equation."
    parser = argparse.ArgumentParser(description=description_text)

    argument_name = "--dir"
    default_value = "."
    help_text = "Directory where the results are saved."

    parser.add_argument(argument_name,
                    default=default_value,
                    help=help_text)
    argument_name = "--save_pdf"
    action_value = "store_true"
    help_text = "Save the visualization as a PDF instead of displaying it."

    parser.add_argument(argument_name,
                    action=action_value,
                    help=help_text)

    argument_name = "--scale_factor"
    argument_type = float
    default_value = 1.0
    help_text = "Factor to scale the wavefunctions for better visualization."

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = "--xlim_left"
    number_of_arguments = 2
    argument_type = float
    default_value = [None, None]
    help_text = "Limits for the x-axis of the left graph."

    parser.add_argument(argument_name,
                    nargs=number_of_arguments,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = "--xlim_right"
    number_of_arguments = 2
    argument_type = float
    default_value = [None, None]
    help_text = "Limits for the x-axis of the right graph."

    parser.add_argument(argument_name,
                    nargs=number_of_arguments,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = "--ylim"
    number_of_arguments = 2
    argument_type = float
    default_value = [None, None]
    help_text = "Limits for the y-axis."

    parser.add_argument(argument_name,
                    nargs=number_of_arguments,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--ytol-start'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the start of the Y-axis.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--ytol-end'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the end of the Y-axis.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--xtol-start-left'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the start of the X-axis for the left plot.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--xtol-end-left'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the end of the X-axis for the left plot.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--xtol-start-right'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the start of the X-axis for the right plot.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)

    argument_name = '--xtol-end-right'
    argument_type = float
    default_value = None
    help_text = 'Tolerance for the end of the X-axis for the right plot.'

    parser.add_argument(argument_name,
                    type=argument_type,
                    default=default_value,
                    help=help_text)


    args = parser.parse_args()

    # Applying the tolerances
    xlim_left = list(args.xlim_left)
    xlim_right = list(args.xlim_right)
    ylim = list(args.ylim)

    if xlim_left[0] is not None and args.xtol_start_left is not None:
        xlim_left[0] += args.xtol_start_left
    if xlim_left[1] is not None and args.xtol_end_left is not None:
        xlim_left[1] -= args.xtol_end_left
    if xlim_right[0] is not None and args.xtol_start_right is not None:
        xlim_right[0] += args.xtol_start_right
    if xlim_right[1] is not None and args.xtol_end_right is not None:
        xlim_right[1] -= args.xtol_end_right
    if ylim[0] is not None and args.ytol_start is not None:
        ylim[0] += args.ytol_start
    if ylim[1] is not None and args.ytol_end is not None:
        ylim[1] -= args.ytol_end

    # Call the visualize_data function
    visualize_data(
    args.dir,
    args.save_pdf,
    args.scale_factor,
    tuple(xlim_left),
    tuple(xlim_right),
    tuple(ylim)
)


if __name__ == "__main__":
    main()
    