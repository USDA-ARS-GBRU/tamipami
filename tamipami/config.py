# -*- coding: utf-8 -*-
import os
import yaml
import logging

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_config(config_file=None):
    """
    Loads a YAML configuration file.

    If no configuration file path is provided, it defaults to 'assets/config.yaml'
    within the root directory. Handles errors related to file not found, YAML parsing,
    and other exceptions, logging them appropriately.

    Parameters:
        config_file (str, optional): The path to the configuration file. Defaults to None.

    Returns:
        dict: The loaded configuration data.

    Raises:
        Exception: If an unexpected error occurs during file loading.
    """
    if config_file is None:
        config_file = os.path.join(ROOT_DIR, "assets", "config.yaml")
    try:
        with open(config_file, "r") as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: The configuration file {config_file} was not found.")
    except yaml.YAMLError as e:
        print(f"Error: Failed to parse YAML file {config_file}: {e}")
    except Exception as e:
        logging.error(f"Failed to load configuration file {config_file}: {e}")
        raise


config = load_config()
