"""Post-processing utilities for wood stiffness homogenization results."""
import json
import os
import itertools

import click
import numpy as np

from .. import constants as cst
from .main import cli


BARS = ['E_R', 'E_T', 'E_L', 'G_TL', 'G_RL', 'G_RT', 'M_cw']
COLOR = {1: '#8c5a2b', 2: '#1f77b4', 3: '#2ca02c'}  # native / infiltrated / delignified + infiltrated
# Results for combination: density=650.0, moisture=12.0, E_poly=4.5, nu_poly=0.4, IF_alpha=0.0, IF_beta=0.0, swelling_pct=10.0, CI=0.58, cellulose=0.41, hemicellulose=0.324, lignin=0.22, extractives=0.032, vessel_pct=25.0, ray_pct=10.5, L_EW=nan, W2_EW=nan, L_LW=nan, W2_LW=nan, cell_aspect_ratio=0.01, superellipse_n=4.0, n_families=20.0, microfibril_slenderness=1e-20, tolerance=8.0
SHORT_TITLE = {
    'density': r'$\rho $ [$kg/m^3$]',
    'moisture': 'M [%]',
    'E_poly': r'$E_{p}$ [GPa]',
    'nu_poly': r'$\nu_{p}$ [-]',
    'IF_alpha': r'$\alpha_{IF}$ [$GPa^{-1}$]',
    'IF_beta': r'$\beta_{IF}$ [$GPa^{-1}$]',
    'swelling_pct': r'$\Delta t$ [%]',
    'CI': 'CI [-]',
    'cellulose': 'C [%]',
    'hemicellulose': 'H [%]',
    'lignin': 'L [%]',
    'extractives': 'X [%]',
    'vessel_pct': 'V [%]',
    'ray_pct': 'R [%]',
    'L_EW': r'$L_{EW}$ [µm]',
    'W2_EW': r'$2W_{EW}$ [µm]',
    'L_LW': r'$L_{LW}$ [µm]',
    'W2_LW': r'$2W_{LW}$ [µm]',
    'cell_aspect_ratio': 'AR [-]',
    'superellipse_n': 'n [-]',
    'n_families': 'N [-]',
    'microfibril_slenderness': 'S [-]',
    'tolerance': 'tol [-]',
}


@cli.group()
def postproc():
    """Post-processing commands for wood stiffness homogenization results."""


@postproc.command()
@click.argument('input_file_csv', type=click.Path(exists=True, dir_okay=False),)
@click.option(
    '-s', '--show-plot',
    is_flag=True,
    default=False,
    help='Show the plot interactively. (Requires a suitable matplotlib backend, e.g., Qt5Agg or TkAgg.)'
)
@click.option(
    '--output-dir',
    type=click.Path(dir_okay=True, writable=True),
    default='output_plots',
    show_default=True,
    help='Path to save the output plot images.'
)
def plot(input_file_csv, show_plot, output_dir):
    """Plot the results of wood stiffness homogenization as a series of bar plots `plot_XXX.png` in the output
    directory. An `indexes.json` file is also created to map the plot files to the corresponding parameter
    combinations."""
    import matplotlib
    import matplotlib.pyplot as plt

    if show_plot:
        messages = []
        for backend in ['Qt5Agg', 'TkAgg']:
            try:
                matplotlib.use(backend)
            except Exception as e:
                messages.append(f"Could not use matplotlib backend '{backend}': {e}")
                continue
            else:
                print(f"Using matplotlib backend '{backend}' for interactive plotting.")
                break
        else:
            for msg in messages:
                click.echo(msg)
            click.echo('Could not set a suitable matplotlib backend to show the plot.')
            show_plot = False

    data = np.genfromtxt(input_file_csv, delimiter=',', names=True, dtype=None, encoding='utf8')

    click.echo(f"Loaded data from {input_file_csv}: {len(data)} rows.")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    lists = []

    for col in cst.INPUTS.keys():
        lists.append([(col, unq) for unq in np.unique(data[col])])

    state_rev_map = {v: k for k, v in cst.STATES.items()}


    combos = list(itertools.product(*lists))

    map_data = {}
    map_file = os.path.join(output_dir, 'indexes.json')
    idx = 0
    for combo in combos:
        row = data
        for col, val in combo:
            typ = row[col].dtype
            if typ == np.float64 and np.isnan(val):
                continue
            row = row[row[col] == val]

        _, ax = plt.subplots()

        width = 0.8 / 3
        for state in row['state']:
            state_idx = state_rev_map[state]
            bar_positions = np.arange(len(BARS)) + (state_idx - 2) * width
            state_row = row[row['state'] == state]
            bar_values = [state_row[bar][0] for bar in BARS]
            ax.bar(
                bar_positions, bar_values, width=width,
                label=f"State: {state}",
                color=COLOR[state_idx]
            )
        ax.set_xticks(range(len(BARS)), BARS)
        ax.set_ylabel('[GPa]')
        ax.grid(alpha=0.3, axis='y')

        title = f"{', '.join(f'{SHORT_TITLE.get(col, col)}={val}' for col, val in combo)}"
        ax.set_title(title, wrap=True)

        ax.legend()

        map_data[f"plot_{idx}"] = dict(combo)

        output_file = os.path.join(output_dir, f"plot_{idx}.png")
        plt.tight_layout()
        plt.savefig(output_file)
        click.echo(f"Saved plot to {output_file}")
        if show_plot:
            plt.show()
        plt.close()
        idx += 1

    with open(map_file, 'w') as f:
        safe = json.dumps(map_data, indent=2).replace('NaN', 'null')
        f.write(safe)


__all__ = ['postproc', 'plot']
