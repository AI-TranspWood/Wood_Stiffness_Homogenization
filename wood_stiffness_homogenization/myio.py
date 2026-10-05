"""Input/Output functions for the wood stiffness homogenization model."""
import numpy as np
import csv

from .constants import INPUTS, OUTPUTS, STATES

def write_csv(batch, filename):
    """One line per row of run_batch: all inputs, outputs, upper triangles of C."""
    iu, f = np.triu_indices(6), batch['inputs']
    with open(filename, 'w', newline='', encoding='utf8') as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ['wood', 'polymer', 'cellulose_material', 'swelling', 'state'] + list(INPUTS) +
            ['MFA'] + list(OUTPUTS) + [f'Cclw_{a + 1}{c + 1}' for a, c in zip(*iu)] +
            [f'Ccellwall_{a + 1}{c + 1}' for a, c in zip(*iu)]
        )
        for d, r, i in batch['rows']:
            writer.writerow(
                [f['wood'], f['polymer'], f['cellulose_material'], f['swelling'], STATES[d['state']]] +
                [d[k] for k in INPUTS] + [d['MFA']] + [r[k][i] for k in OUTPUTS] +
                list(r['C_clearwood'][i][iu]) + list(r['C_cellwall'][i][iu])
            )
