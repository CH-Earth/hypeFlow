# hypeFlow
A model specific workflow for HYPE

## Installation

Install the package using pip:

```bash
pip install -e .
```

## Usage

### Command Line Interface (CLI)

The hypeFlow CLI allows you to generate HYPE model setup files from a JSON configuration file.

#### Basic Usage

```bash
hypeflow --json config.json
```

#### Override Output Path

```bash
hypeflow --json config.json --output-path /path/to/custom/output
```

#### Verbose Output

```bash
hypeflow --json config.json --verbose
```

#### Help

```bash
hypeflow --help
```

### Configuration File

Create a JSON configuration file with the following structure (see `example_config.json`):

```json
{
    "easymore_output": "/path/to/easymore/outputs",
    "output_path": "/path/to/hype/output",
    "timeshift": -6,
    "forcing_units": {
        "temperature": {
            "in_varname": "CaSR_v3.1_P_TT_09975",
            "in_units": "celsius",
            "out_units": "celsius"
        },
        "precipitation": {
            "in_varname": "CaSR_v3.1_A_PR0_SFC",
            "in_units": "m/hr",
            "out_units": "mm/day"
        }
    },
    "geofabric_mapping": {
        "basinID": {
            "in_varname": "COMID"
        },
        "nextDownID": {
            "in_varname": "NextDownID"
        },
        "area": {
            "in_varname": "unitarea",
            "in_units": "km^2",
            "out_units": "m^2"
        },
        "rivlen": {
            "in_varname": "new_len_km",
            "in_units": "km",
            "out_units": "m"
        }
    },
    "gistool_outputs": {
        "soil": "/path/to/soil/stats.csv",
        "landcover": "/path/to/landcover/stats.csv",
        "elevation": "/path/to/elevation/stats.csv"
    },
    "subbasins_shapefile": "/path/to/subbasins.shp",
    "rivers_shapefile": "/path/to/rivers.shp",
    "frac_threshold": 0.01,
    "spinup_days": 274
}
```

### Python API

You can also use the functions directly in Python:

```python
import hypeflow as hf

# Configuration parameters
easymore_output = '/path/to/easymore/outputs'
output_path = '/path/to/hype/output'
timeshift = -6
# ... other parameters

# Run the workflow
hf.write_hype_forcing(easymore_output, timeshift, forcing_units, geofabric_mapping, output_path)
hf.write_hype_geo_files(gistool_outputs, subbasins_shapefile, rivers_shapefile, frac_threshold, geofabric_mapping, output_path)
hf.write_hype_par_file(output_path)
hf.write_hype_info_filedir_files(output_path, spinup_days)
```

## CLI Options

- `--json`: Path to the JSON configuration file (required)
- `--output-path`: Override the output path specified in the JSON file (optional)
- `--verbose`, `-v`: Enable verbose output (optional)

The CLI will execute the following steps in sequence:
1. Write HYPE forcing files
2. Write HYPE geometry files
3. Write HYPE parameter file
4. Write HYPE info and filedir files

# License
Copyright (c) 2023 Mohamed Moghairib.
