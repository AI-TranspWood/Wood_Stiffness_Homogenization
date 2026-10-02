"""Run the wood stiffness homogenization."""
import json
import os

from ..homogenization_tw import run_batch, write_csv
from .main import cli, click


@cli.command()
@click.argument('input_json', type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option(
    '--output', '-o',
    default=None,
    type=click.Path(dir_okay=False, writable=True), help='Output file name (CSV).'
)
def run(input_json, output):
    """Run the wood stiffness homogenization."""
    click.echo(f'Running wood stiffness homogenization with input file: {input_json}')
    with open(input_json, encoding='utf8') as f:
        input_data = json.load(f)
    if output is None:
        input_dir = os.path.dirname(input_json)
        output_name = input_data.get('output')
        if output_name is None:
            output_name = os.path.splitext(os.path.basename(input_json))[0] + '_results.csv'
        output = os.path.join(input_dir, output_name)

    result = run_batch(input_data)
    write_csv(result, output)

    click.echo(f'Output written to: {output}')
    for err in result.get('errors', []):
        click.echo(f'Error: {err}')

__all__ = ['run']
