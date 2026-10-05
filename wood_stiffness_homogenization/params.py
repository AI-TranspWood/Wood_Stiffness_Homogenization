"""Input paramemters"""
import re
from collections.abc import Callable
from copy import copy
from dataclasses import dataclass, field, fields, MISSING, make_dataclass
import json
from typing import ClassVar, Self

import numpy as np
import rich_click as click

from .constants import CRYCEL, POLYMERS, WOODS, INPUTS


class DelimitedList(click.ParamType):
    """A custom Click parameter type that parses a comma-separated list of values."""

    name = 'list'

    def __init__(self, delimiter=',', subtype=str, exact_length=None):
        super().__init__()
        self.delimiter = delimiter
        self.subtype = subtype
        self.exact_length = exact_length

        name = str(subtype.__name__)
        self.name = f'[{name}[{self.delimiter}{name}]]'

    def convert(self, value, param, ctx):
        print(f'----Converting value: {value} for param: {param.name}')
        if isinstance(value, list):
            res = [self.subtype(item) for item in value]
        elif not isinstance(value, str):
            self.fail(f"Expected a string or list, got {type(value).__name__}", param, ctx)

        mch = re.match(r'^([\d\.\+\-eE]+):([\d\.\+\-eE]+):([\d\.\+\-eE]+)$', value)
        if mch:
            try:
                start, stop, step = map(float, mch.groups())
            except Exception as e:
                self.fail(f"Could not parse range `{value}`: {e}", param, ctx)
            # Make range right-inclusive by adding step times a small eps to the stop value
            range = np.round(np.arange(start, stop + step * 1e-9, step), decimals=12)
            res = [self.subtype(x) for x in range]
        else:
            try:
                items = value.split(self.delimiter)
                res = [self.subtype(item.strip()) for item in filter(None, items)]

                if not res:
                    res = None
            except Exception as e:
                self.fail(f"Could not parse list `{value}`: {e}", param, ctx)

        if self.exact_length is not None and len(res) != self.exact_length:
            self.fail(f"Expected exactly {self.exact_length} items, got {len(res)}", param, ctx)

        return res


class JsonParams:
    params_map: ClassVar[dict[str, str]] = {}
    post_set: ClassVar[list[str]] = []

    def _to_json(self, data: dict) -> dict:
        """Convert the parameters to a JSON serializable dictionary"""
        return data

    def to_json(self, json_file: str):
        """Save the parameters to a JSON file"""
        data = {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
        data = self._to_json(data)

        with open(json_file, 'w') as f:
            json.dump(data, f, indent=4)

    @classmethod
    def from_json(cls, json_file: str) -> list[Self]:
        """Create an instance from a JSON file"""
        with open(json_file, 'r') as f:
            data = json.load(f)
        res = []
        if isinstance(data, dict):
            res = [cls.from_dict(data)]
        elif isinstance(data, list):
            res = [cls.from_dict(item) for item in data]
        else:
            raise ValueError('Invalid data format in JSON file')
        return res

    @classmethod
    def from_dict(cls, data: dict) -> Self:
        """Create an instance from a JSON file"""
        post = {}
        data = {cls.params_map.get(k, k): v for k, v in data.items()}
        for k in cls.post_set:
            if k in data:
                post[k] = data.pop(k)
        res = cls(**data)
        for k, v in post.items():
            setattr(res, k, v)
        return res

    @classmethod
    def to_click_options(cls, func: Callable) -> Callable:
        """Decorator to add click options for the parameters to a click command"""
        name_map = {}

        OVERRIDE_GROUP = 'Override Input Parameters'

        groups = {'Options'}

        def callback(ctx: click.Context, param: click.Parameter, value):
            ctx.ensure_object(dict)
            # source = ctx.get_parameter_source(param.name)
            # if source == ParameterSource.DEFAULT:
            #     return
            # print(f'Callback for {param.name:<25s} with value: {value}')
            name = name_map.get(param.name)
            overrides = ctx.obj.setdefault('override_params', {})
            if value is None:
                # print(f'    No value provided for {name}, skipping override.')
                return
            print(f"    Overriding parameter {name:<25s} with value: {value}")
            overrides[name] = value

        for fld in fields(cls)[::-1]:
            extra_help = ''
            metadata = getattr(fld, 'metadata', {})
            if not fld.init:
                continue
            if fld.name.startswith('_'):
                continue
            if fld.type not in (
                    int, float, str, bool,
                    list[str], tuple[int, int, int], list[int], list[float]
                ):
                continue

            default = fld.default
            if fld.default_factory is not MISSING:
                default = fld.default_factory()
            elif default is not MISSING:
                default = fld.default
            else:
                default = None

            typ = fld.type
            # print(f'Processing field: {fld.name} with type: {typ} and default: {default}')
            if typ in (int, float):
                min_val = metadata.get('min', None)
                max_val = metadata.get('max', None)
                if min_val is not None or max_val is not None:
                    cls_typ = click.FloatRange if typ == float else click.IntRange
                    typ = cls_typ(min=min_val, max=max_val)
            elif typ == str:
                if metadata.get('file', False):
                    typ = click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True)
                elif metadata.get('dir', False):
                    typ = click.Path(exists=True, file_okay=False, readable=True, resolve_path=True)
                elif metadata.get('choices', None) is not None:
                    typ = click.Choice(metadata['choices'])
            elif typ == list[str]:
                delimiter = metadata.get('delimiter', ',')
                default = delimiter.join(default) if default else None

                typ = DelimitedList(delimiter=delimiter, subtype=str)
                extra_help = ' (comma-separated list)'
            elif typ == list[int]:
                delimiter = metadata.get('delimiter', ',')
                default = delimiter.join(str(x) for x in default) if default else None

                typ = DelimitedList(delimiter=metadata.get('delimiter', ','), subtype=int)
                extra_help = ' (comma-separated list of integers)'
            elif typ == list[float]:
                delimiter = metadata.get('delimiter', ',')
                default = delimiter.join(str(x) for x in default) if default else None

                typ = DelimitedList(delimiter=metadata.get('delimiter', ','), subtype=float)
                extra_help = ' (comma-separated list of floats)'
            elif typ == tuple[int, int, int]:
                typ = DelimitedList(subtype=int, exact_length=3)
                extra_help = ' (comma-separated list of 3 integers)'


            expose = metadata.get('expose_value', False)
            # prefix = '--param-' if not expose else '--'
            prefix = '--'

            decl = f'{prefix}{fld.name}'
            if fld.type == bool:
                decl = f'{prefix}{fld.name}/{prefix}no-{fld.name}'
            name_map[fld.name.lower()] = fld.name

            required = metadata.get('required', False)

            group = metadata.get('group', OVERRIDE_GROUP)
            groups.add(group)

            help_str = metadata.get('help', None)
            if help_str is not None and extra_help:
                help_str += extra_help

            kwargs = {
                'type': typ,
                'default': default,
                'is_flag': fld.type == bool,
                'required': required,
                'expose_value': expose,
                'callback': callback,
                'help': help_str,
                'panel': group,
                'show_default': True,
            }

            if 'flag_value' in metadata:
                kwargs['flag_value'] = metadata['flag_value']

            func = click.option(decl, **kwargs)(func)

        func = click.option(
            '--config-file', type=click.Path(exists=True), help='Path to file with parameters'
        )(func)

        for group in sorted(groups, key=lambda x: (x != 'Options', x))[::-1]:
            func = click.option_panel(group)(func)

        return func

@dataclass(kw_only=True)
class HomegenizationParams(JsonParams):
    """Define the parameters for wood stiffness homogenization"""
    output: str = field(
        default=None,
        metadata={
            'help': 'CSV file path to save the homogenization results',
            'group': 'Options',
        }
    )

    ######################################################################
    # Wood properties
    wood: str = field(
        default='Birch',
        metadata={
            'help': 'Wood species name as defined in woods.json',
            'group': '1. Wood Properties',
            'choices': list(WOODS.keys()),
        }
    )
    density: str = field(
        metadata={
            'help': 'Native wood density [kg/m³] as a range or single value',
            'group': '1. Wood Properties',
        }
    )
    moisture: float = field(
        metadata={
            'help': 'Moisture content value [% dry mass]',
            'group': '1. Wood Properties',
        }
    )

    ######################################################################
    # Polymers and interface
    polymer: str = field(
        default='HEMA',
        metadata={
            'help': 'Polymer name as defined in polymers.json (sets E_poly, nu_poly unless these are given)',
            'group': '2. Polymers and interface',
            'choices': list(POLYMERS.keys()),
        }
    )
    E_poly: float = field(
        metadata={
            'help': 'Polymer Young\'s modulus values [GPa]',
            'group': '2. Polymers and interface',
        }
    )
    nu_poly: float = field(
        metadata={
            'help': 'Polymer Poisson\'s ratio',
            'group': '2. Polymers and interface',
            'min': 0.0, 'max': 0.5,
        }
    )
    IF_alpha: float = field(
        metadata={
            'help': 'Fibril–matrix interface compliance, tangential [1/GPa]',
            'group': '2. Polymers and interface',
        }
    )
    IF_beta: float = field(
        metadata={
            'help': 'Fibril–matrix interface compliance, normal [1/GPa]',
            'group': '2. Polymers and interface',
        }
    )


    ######################################################################
    # Swelling
    swelling: bool = field(
        # default=True,
        metadata={
            'help': 'Wall swelling on/off',
            'group': '3. Swelling Properties',
        }
    )
    swelling_pct: float = field(
        # default_factory=lambda: [10.0],
        metadata={
            'help': 'Increase of the wall thickness [%]',
            'group': '3. Swelling Properties',
            'min': 0.0, 'max': 100.0,
        }
    )

    ######################################################################
    # Advanced parameters
    cellulose_material: str = field(
        default='Dri2014_TI',
        metadata={
            'help': 'Cellulose material model: Dri2014_TI, IMWS or Dri',
            'group': '4. Advanced Parameters',
            'choices': list(CRYCEL.keys()),
        }
    )
    CI: float = field(
        default=None,
        metadata={
            'help': 'Mass-based crystallinity of cellulose [%]',
            'group': '4. Advanced Parameters',
        }
    )

    cellulose: float = field(
        metadata={
            'help': 'Dry-mass fraction of cellulose in the cell wall [%]',
            'group': '4. Advanced Parameters',
        }
    )
    hemicellulose: float = field(
        metadata={
            'help': 'Dry-mass fraction of hemicellulose in the cell wall [%]',
            'group': '4. Advanced Parameters',
        }
    )
    lignin: float = field(
        metadata={
            'help': 'Dry-mass fraction of lignin in the cell wall [%]',
            'group': '4. Advanced Parameters',
        }
    )
    extractives: float = field(
        metadata={
            'help': 'Dry-mass fraction of extractives in the cell wall [%]',
            'group': '4. Advanced Parameters',
        }
    )
    vessel_pct: float = field(
        metadata={
            'help': 'Vessel volume fraction [%]',
            'group': '4. Advanced Parameters',
        }
    )
    ray_pct: float = field(
        metadata={
            'help': 'Ray volume fraction [%]',
            'group': '4. Advanced Parameters',
        }
    )
    L_EW: float = field(
        metadata={
            'help': 'Earlywood lumen width [µm]; null for hardwoods',
            'group': '4. Advanced Parameters',
        }
    )
    W2_EW: float = field(
        metadata={
            'help': 'Earlywood double wall thickness [µm]; null for hardwoods',
            'group': '4. Advanced Parameters',
        }
    )
    L_LW: float = field(
        metadata={
            'help': 'Latewood lumen width [µm]; null for hardwoods',
            'group': '4. Advanced Parameters',
        }
    )
    W2_LW: float = field(
        metadata={
            'help': 'Latewood double wall thickness [µm]; null for hardwoods',
            'group': '4. Advanced Parameters',
        }
    )

    cell_aspect_ratio: float = field(
        metadata={
            'help': 'Fibre-lumen slenderness (the ray-cell lumen is fixed at 1/5)',
            'group': '4. Advanced Parameters',
            'min': 0.0, 'max': 1.0,
        }
    )
    superellipse_n: float = field(
        metadata={
            'help': 'Lumen shape exponent (the ray-cell lumen is fixed at 1/5)',
            'group': '4. Advanced Parameters',
            'min': 1.0, 'max': 10.0,
        }
    )
    n_families: int = field(
        metadata={
            'help': 'Number of microfibril families',
            'group': '4. Advanced Parameters',
            'min': 1, 'max': 100,
        }
    )
    microfibril_slenderness: float = field(
        metadata={
            'help': 'Microfibril slenderness',
            'group': '4. Advanced Parameters',
            'min': 1e-30, 'max': 1e-10,
        }
    )
    tolerance: float = field(
        metadata={
            'help': 'Tolerance for numerical calculations',
            'group': '4. Advanced Parameters',
            'min': 1e-10,
        }
    )

def list_fields(cls: type, *field_names: str) -> type:
    """Create a new dataclass with the same fields as cls, but with the specified fields converted to lists."""
    selected = set(field_names)
    new_fields = []

    for fld in fields(cls):
        new_fld = copy(fld)

        type_ = fld.type
        if fld.name in selected:
            type_ = list[fld.type]
            new_fld.type = type_

        new_fields.append((fld.name, type_, new_fld))

    return make_dataclass(
        cls.__name__ + "List", new_fields,
        bases=(cls,), namespace=dict(cls.__dict__),
        kw_only=True
    )

BatchHomegenizationParams = list_fields(HomegenizationParams, *INPUTS.keys())
