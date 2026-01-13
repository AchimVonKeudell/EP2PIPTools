import argparse
from pathlib import Path
from generate_campaign_config import list_command, configure_command
import ftir_data_processing
"""The project can be run by calling this script. 
There are two modes of operation: 
    --config: following the settings of a single `config.yaml` file (which is given with the command)
    --batch: reads all config.yaml files in a measurement_campaign directory
        \output\`measurement campaign`
            \data_set_1\config.yaml
            \data_set_2\config.yaml
        these are previously configured using `generate_campaign_config.py`
    --configure: configures .yaml files following config_selection_rules.yaml
        
Steijn Vervloedt
25/11-2025, Bochum
"""



### running the project done below:
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        usage="example: python run.py [option] [filepath or directory]",
        description='IR Spectrum Analysis'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--batch', type=Path, help='Process all configs in directory')
    group.add_argument('--config', type=Path, help='Single config file')
    group.add_argument('--configure', type=Path,
                       help="Configures the given directory with the config.yaml files following the rules in "
                            "[config_selection_rules.yaml]. The rules follow the name of either the measurement campaign "
                            "or the respective dataset.")
    group.add_argument('--list', type=Path, help='Listing all the directories containing Bruker OPUS files')

    args = parser.parse_args()
    print(f'{args=}')

    if args.list:
        list_command(args.list)
    elif args.configure:
        configure_command(args.configure)
    elif args.config:
        print(args.config)
        ftir_data_processing.run(args.config)
    elif args.batch:
        print(f'{args.batch=}')
        # Find all config.yaml files in subdirectories
        config_files = sorted(args.batch.glob('**/config.yaml'))

        if not config_files:
            print(f"No config.yaml files found in subdirectories of {args.batch}")
        else:
            print(f"Found {len(config_files)} config files:")

            for config_file in config_files:
                print(f"  - {config_file}")

            print('Continue with analysis:\n')

            for config_file in config_files:
                ftir_data_processing.run(config_file)