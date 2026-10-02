"""Constants and static input data for the wood stiffness homogenization model."""
import json
from importlib import resources

import numpy as np


def load_resource_file(filename):
    path = resources.files('wood_stiffness_homogenization') / 'data' / filename
    with open(path, encoding='utf8') as f:
        return json.load(f)


# ---------------------------------------------------- input files and batches
# Flat input keys (GUI fields = input-file keys). Numeric values may be a number,
# a list [a, b, c] or a range string 'start:step:stop' -> all combinations are run.
INPUTS = {
    'density': 'Density [kg/m³]',
    'moisture': 'Moisture [% dry mass]',
    'E_poly': 'E polymer [GPa]',
    'nu_poly': 'nu polymer [-]',
    'IF_alpha': 'IF alpha tangential [1/GPa]',
    'IF_beta': 'IF beta normal [1/GPa]',
    'swelling_pct': 'Wall thickness increase [%]',
    'CI': 'Crystallinity CI [-]',
    'cellulose': 'Cellulose [dry mass fr.]',
    'hemicellulose': 'Hemicellulose [dry mass fr.]',
    'lignin': 'Lignin [dry mass fr.]',
    'extractives': 'Extractives [dry mass fr.]',
    
    'vessel_pct': 'Vessels [%]',
    'ray_pct': 'Rays [%]',
    'L_EW': 'L EW [µm]',
    'W2_EW': '2W EW [µm]',
    'L_LW': 'L LW [µm]',
    'W2_LW': '2W LW [µm]',
    'cell_aspect_ratio': 'Cell aspect ratio',
    'superellipse_n': 'Superellipse n',
    'n_families': 'Fibril orientation families',
    'microfibril_slenderness': 'Microfibril slenderness',
    'tolerance': 'Tolerance (1-16)'
}
STATES = {1: 'native', 2: 'infiltrated', 3: 'delignified + infiltrated'}
OUTPUTS = {
    'E_R': 'E_R [GPa]',
    'E_T': 'E_T [GPa]',
    'E_L': 'E_L [GPa]',
    'G_TL': 'G_TL [GPa]',
    'G_RL': 'G_RL [GPa]',
    'G_RT': 'G_RT [GPa]',
    'nu_TR': 'nu_TR',
    'nu_LR': 'nu_LR',
    'nu_TL': 'nu_TL',
    'nu_LT': 'nu_LT',
    'M_cw': 'M_cellwall [GPa]'
}
POLY = ('E_poly', 'nu_poly')                      # irrelevant for native wood

PH = load_resource_file('phases.json')                         # cell-wall constituents
WOODS = load_resource_file('woods.json')                       # species data (from data_Wood.xlsx)
POLYMERS = {k: (v['E_GPa'], v['nu']) for k, v in load_resource_file('polymers.json').items() if k[0] != '_'}
CRYCEL = {k: np.array(v['C_GPa']) for k, v in PH['crystalline_cellulose'].items() if isinstance(v, dict)}
RHO = dict(crycel=PH['crystalline_cellulose']['density'], amcel=PH['amorphous_cellulose']['density'],
           hemcel=PH['hemicellulose']['density'], lignin=PH['lignin']['density'],
           extr=PH['extractives']['density'])
