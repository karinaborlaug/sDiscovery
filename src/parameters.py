#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions to extract parameters from config.yaml file."""

import yaml
import os

try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


def load_config(config_file):
    """Load config file with yaml.

    Args:
        config_file: config_file

    Returns:
           config file

    """
    with open(config_file, 'r') as stream:
        config = yaml.load(stream, Loader=Loader)

    return config


def get_strict_filter_parameters(config_file):
    """Get strict filter from config file.

    Args:
        config_file: config_file

    Returns:
           strict filter

    """
    config = load_config(config_file)
    return config['strict_filter']


def get_loose_filter_parameters(config_file):
    """Get loose filter from config file.

    Args:
        config_file: config_file

    Returns:
           loose filter

    """
    config = load_config(config_file)
    return config['loose_filter']


def get_sample_ids(config_file):
    """Get sample ids from config file.

    Args:
        config_file: config_file

    Returns:
           sample ids

    """
    config = load_config(config_file)

    return config['samples'].keys()


def get_result_dir(config_file):
    """Get path to results directory from config file.

    Args:
        config_file: config_file

    Returns:
           results directory

    """
    config = load_config(config_file)

    return os.path.join(config['working_dir'], config['results_dir'])


def get_vcf_files(config_file):
    """Get path to vcf files from config file.

    Args:
        config_file: config_file

    Returns:
           path to vcf files

    """
    vcf_files = {}
    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        if os.path.exists(config['samples'][sample_id]['vcf_file']):
            vcf_files[sample_id] = config['samples'][sample_id]['vcf_file']
        else:
            vcf_files[sample_id] = os.path.join(config['working_dir'], config['samples'][sample_id]['vcf_file'])

    return vcf_files


def get_maf_files(config_file):
    """Get path to maf files from config file.

    Args:
        config_file: config_file

    Returns:
           path to maf files

    """
    maf_files = {}
    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        if os.path.exists(config['samples'][sample_id]['maf_file']):
            maf_files[sample_id] = config['samples'][sample_id]['maf_file']
        else:
            maf_files[sample_id] = os.path.join(config['working_dir'], config['samples'][sample_id]['maf_file'])

    return maf_files


def get_type(config_file):
    """Get sampel types from config file.

    Args:
        config_file: config_file

    Returns:
           sample types

    """
    types = {}
    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        types[sample_id] = config['samples'][sample_id]['type']

    return types
