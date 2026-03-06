#!/usr/bin/env python

import os
import sys
import json
from genome3danalysis.structfeat import SfFile


def main():
    """ Run the Structural Feature (SF) pipeline.
    
    This function gets the configuration json file name from the command line.
    
    Then it loads the configuration, initializes the SF file and runs it.
    """
    
    # Check that the command line has the correct number of arguments
    if len(sys.argv) != 2:
        raise ValueError('Error: structfeat-run requires one argument')

    # Get the configuration json file name from the command line
    config_name = sys.argv[1]

    # Check that configuration file exists and is a json file
    if not os.path.isfile(config_name):
        raise FileNotFoundError(f'Error: {config_name} not found')
    if not config_name.endswith('.json'):
        raise ValueError(f'Error: {config_name} is not a json file')

    # Convert the configuration file name to an absolute path
    config_name = os.path.abspath(config_name)

    # Load the configuration json file
    with open(config_name, 'r') as f:
        config = json.load(f)

    # Get the hss file name from the configuration
    hss_name = config['hss_name']

    # Resolve output SF file name
    # Priority:
    # 1) config['sf_name'] (explicit full path)
    # 2) config['output_dir'] + basename(hss_name).replace('.hss', '.sf.h5')
    # 3) default: hss_name.replace('.hss', '.sf.h5')
    if 'sf_name' in config:
        sf_name = config['sf_name']
    else:
        sf_basename = os.path.basename(hss_name).replace('.hss', '.sf.h5')
        if 'output_dir' in config:
            sf_name = os.path.join(config['output_dir'], sf_basename)
        else:
            sf_name = hss_name.replace('.hss', '.sf.h5')

    sf_name = os.path.abspath(sf_name)
    os.makedirs(os.path.dirname(sf_name), exist_ok=True)

    sf = SfFile(sf_name, 'w')

    # Run the SF file
    sf.run(config)

    # Close the SF file
    sf.close()


if __name__ == '__main__':
    main()
