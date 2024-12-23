# -*- coding: utf-8 -*-
# #!/usr/bin/local python

import yaml

def load_config(config_file='assets/config.yaml'):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

config = load_config()