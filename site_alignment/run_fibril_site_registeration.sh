#!/bin/bash

# Activate "fibrilsite" environment

python ./fibril_site_registeration.py --info_file "./2025-04-07_sites_parsed_info/2025-04-07_all_sites_input_feats.csv" \
                                                   --sites_folder "../sel_fibril_sites/"  \
                                                   --output_folder "."