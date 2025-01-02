# -*- coding: utf-8 -*-
# #!/usr/bin/local python

import yaml
import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
def load_config(config_file=os.path.join(ROOT_DIR, 'assets/config.yaml')):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

config = load_config()