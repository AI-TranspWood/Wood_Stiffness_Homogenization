"""Run the wood stiffness homogenization."""
import json
import os

from ..homogenization_tw import run_batch
from ..myio import write_csv
from ..params import BatchHomegenizationParams
from .main import cli, click


@cli.command()
@click.pass_context
@BatchHomegenizationParams.to_click_options
def run(ctx, config_file):
    """Run the wood stiffness homogenization."""
    click.echo(f'Running wood stiffness homogenization with input file: {config_file}')

    data = {}
    if config_file:
         with open(config_file, encoding='utf8') as f:
            data = json.load(f)

    overrides: dict = ctx.obj.get('override_params', {})

    output = overrides.pop('output', None)
    if output is None:
        if config_file is None:
            input_dir = os.getcwd()
            output_name = 'homogenization_results.csv'
        else:
            input_dir = os.path.dirname(config_file)
            output_name = os.path.splitext(os.path.basename(config_file))[0] + '_results.csv'
        output = os.path.join(input_dir, output_name)

    if overrides:
        data.update(overrides)

    result = run_batch(data)
    write_csv(result, output)

    click.echo(f'Output written to: {output}')
    for err in result.get('errors', []):
        click.echo(f'Error: {err}')

__all__ = ['run']
