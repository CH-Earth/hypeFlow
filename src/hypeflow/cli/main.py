#!/usr/bin/env python
"""
HypeFlow CLI - Command Line Interface for HYPE model setup

This CLI accepts input arguments from a JSON file and runs the HYPE setup functions in sequence.
"""

import json
import os
import sys
from pathlib import Path

import click

import hypeflow as hf


def validate_json_file(ctx, param, value):
    """Validate that the JSON file exists and is readable."""
    if value is None:
        return None
    
    json_path = Path(value)
    if not json_path.exists():
        raise click.BadParameter(f"JSON file '{value}' does not exist.")
    
    if not json_path.is_file():
        raise click.BadParameter(f"'{value}' is not a file.")
    
    try:
        with open(json_path, 'r') as f:
            json.load(f)
    except json.JSONDecodeError as e:
        raise click.BadParameter(f"Invalid JSON file: {e}")
    except Exception as e:
        raise click.BadParameter(f"Error reading JSON file: {e}")
    
    return value


def load_config(json_file):
    """Load configuration from JSON file."""
    with open(json_file, 'r') as f:
        config = json.load(f)
    
    # Validate required fields
    required_fields = [
        'easymore_output',
        'timeshift',
        'forcing_units',
        'geofabric_mapping',
        'gistool_outputs',
        'subbasins_shapefile',
        'rivers_shapefile',
        'frac_threshold',
        'spinup_days'
    ]
    
    missing_fields = [field for field in required_fields if field not in config]
    if missing_fields:
        raise click.ClickException(f"Missing required fields in JSON file: {', '.join(missing_fields)}")
    
    return config


@click.command()
@click.option(
    '--json',
    'json_file',
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    callback=validate_json_file,
    help='Path to the JSON file containing input configuration.'
)
@click.option(
    '--output-path',
    'output_path_override',
    type=click.Path(file_okay=False, dir_okay=True),
    help='Override the output path specified in the JSON file.'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose output.'
)
def cli(json_file, output_path_override, verbose):
    """
    HypeFlow CLI - Generate HYPE model setup files from configuration.
    
    This command reads configuration from a JSON file and runs the HYPE setup
    functions in sequence to generate forcing files, geometry files, parameter
    files, and info files.
    
    Example usage:
        hypeflow --json config.json --output-path /path/to/output
    """
    try:
        # Load configuration
        if verbose:
            click.echo(f"Loading configuration from: {json_file}")
        
        config = load_config(json_file)
        
        # Override output path if provided
        if output_path_override:
            config['output_path'] = output_path_override
            if verbose:
                click.echo(f"Output path overridden to: {output_path_override}")
        
        # Ensure output_path is set
        if 'output_path' not in config and output_path_override is None:
            raise click.ClickException("Output path must be specified either in JSON file or via --output-path option")
        
        output_path = config.get('output_path', output_path_override)
        
        # Create output directory if it doesn't exist
        os.makedirs(output_path, exist_ok=True)
        
        if verbose:
            click.echo("Starting HYPE model setup...")
        
        # Step 1: Write HYPE forcing files
        click.echo("Step 1/4: Writing HYPE forcing files...")
        if verbose:
            click.echo(f"  - Easymore output: {config['easymore_output']}")
            click.echo(f"  - Time shift: {config['timeshift']} hours")
        
        hf.write_hype_forcing(
            easymore_output=config['easymore_output'],
            timeshift=config['timeshift'],
            forcing_units=config['forcing_units'],
            geofabric_mapping=config['geofabric_mapping'],
            path_to_save=output_path
        )
        click.echo("✓ HYPE forcing files written successfully")
        
        # Step 2: Write HYPE geometry files
        click.echo("Step 2/4: Writing HYPE geometry files...")
        if verbose:
            click.echo(f"  - Subbasins shapefile: {config['subbasins_shapefile']}")
            click.echo(f"  - Rivers shapefile: {config['rivers_shapefile']}")
            click.echo(f"  - Fraction threshold: {config['frac_threshold']}")
        
        hf.write_hype_geo_files(
            gistool_outputs=config['gistool_outputs'],
            subbasins_shapefile=config['subbasins_shapefile'],
            rivers_shapefile=config['rivers_shapefile'],
            frac_threshold=config['frac_threshold'],
            geofabric_mapping=config['geofabric_mapping'],
            path_to_save=output_path
        )
        click.echo("✓ HYPE geometry files written successfully")
        
        # Step 3: Write HYPE parameter file
        click.echo("Step 3/4: Writing HYPE parameter file...")
        hf.write_hype_par_file(path_to_save=output_path)
        click.echo("✓ HYPE parameter file written successfully")
        
        # Step 4: Write HYPE info and filedir files
        click.echo("Step 4/4: Writing HYPE info and filedir files...")
        if verbose:
            click.echo(f"  - Spinup days: {config['spinup_days']}")
        
        hf.write_hype_info_filedir_files(
            path_to_save=output_path,
            spinup_days=config['spinup_days']
        )
        click.echo("✓ HYPE info and filedir files written successfully")
        
        click.echo("\n🎉 HYPE model setup completed successfully!")
        click.echo(f"Output files written to: {output_path}")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    cli()