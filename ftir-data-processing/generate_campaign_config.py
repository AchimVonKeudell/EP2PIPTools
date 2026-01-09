"""This script will select the correct .yaml files for each data set and initialise the \output-director.

It should be called as: python generate_campaign_config.py \raw\`measurement_campaign`.
    1. it will read the different data sets
    2. it will create the corresponding directory structure under \output\`measurement_campaign`
    3. it will select the correct .yaml files for each data set
        following the defined rules -> use the naming '{index}_{label}'
        where index is used to identify the order of the experiment,
        `label` refers to which .yaml file should be used

Afterwards, the individual setting can be tweaked for each data_set. This should give sufficient freedom to
accurately process all types of measured spectra.

It should return Errors/Warnings when config.yaml files already exist, since the corresponding .yaml files
may already have been personalised.
"""
import os

import yaml
from pathlib import Path

BRUKER_OPUS_HEADER_ID = b'\xfe\xfe\x00\x00\x00\x00'
CONFIG_RULES_FILEPATH = 'config_selection_rules.yaml'

### base functions
def _is_bruker_files_in_directory(directory: Path) -> bool:
    if not directory.is_dir():
        return False

    for child in directory.iterdir():
        if child.is_file():
            with child.open('rb') as f_data:
                _ = f_data.read(2)
                header = f_data.read(6)
            if header == BRUKER_OPUS_HEADER_ID:
                return True  # found a OPUS file
    return False  # did not find a OPUS file


def _search_directory(directory: Path) -> Path:
    if directory.is_dir():
        for _sub_directory in directory.iterdir():
            if _is_bruker_files_in_directory(_sub_directory):
                return _sub_directory
            _search_directory(_sub_directory)
    return


class KeywordConfigSelector:
    def __init__(self, rules_file=CONFIG_RULES_FILEPATH):
        with open(rules_file, 'r') as f:
            self._config = yaml.safe_load(f)

        self.path_levels = self._config['path_structure']['levels']

        self.rules = self._config['rules']
        for rule in self.rules:
            if 'name' not in rule:
                raise ValueError(f'add name to {rule=}')

            if 'keywords' not in rule:
                raise ValueError(f'add keywords to {rule=}')

            if 'config' not in rule:
                raise ValueError(f'add config to {rule=}')

        self.default_config_file = Path(self._config['default_config'])

    def _get_path_components(self, path: Path) -> dict:
        dataset_path_parts = path.parts

        components = {}
        for label, location in self.path_levels.items():
            components[label] = dataset_path_parts[location]
        return components


    def select_config(self, dataset_path: Path) -> Path:
        """Select config based on keywords in directory path.

        :param dataset_path: Path to dataset (e.g. D:/data/raw/campaign/dataset/FTIR)
        :return: Path to appropriate config file
        """
        # Extract path components
        dataset_components = self._get_path_components(dataset_path)

        # Check rules in priority order
        for rule in self.rules:
            matches = True

            for label, keyword in rule['keywords'].items():
                if keyword:
                    if not (keyword in dataset_components[label]):
                        matches = False

            # found a match - using the corresponding given config file
            if matches:
                return Path(rule['config'])

        # No match - use default
        return self.default_config_file

    def adjust_data_path_location(self, config_file: Path, dataset_path: Path):

        dataset_components = self._get_path_components(dataset_path)

        output_directory = '\\'.join(dataset_path.parts[:self.path_levels['campaign']-1])
        output_directory += '\\output'

        sorted_order_output_director = dict(
            sorted(self._config['path_structure']['output'].items(), key=lambda item: item[1])
        )

        for label in sorted_order_output_director.keys():
            output_directory += '\\' + dataset_components[label]

        print('-', output_directory)
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        with open(config_file, 'r') as fd:
            metadata = yaml.safe_load(fd)

        metadata['input data location'] = str(dataset_path)
        metadata['output data location'] = output_directory

        config_file_output = f'{output_directory}\\config.yaml'
        print('-', config_file_output)
        if os.path.exists(config_file_output):
            print('\talready exists, so not remaking it')
        with open(config_file_output, "w", encoding="utf-8") as f:
            yaml.safe_dump(metadata, f, sort_keys=False, allow_unicode=True)
        return

### commands
def list_command(argument: Path):
    print('Found:')
    for directory in argument.iterdir():
        ftir_dir = _search_directory(directory)
        if ftir_dir:
            print('-', ftir_dir)
    return


def configure_command(argument: Path):

    selector = KeywordConfigSelector()

    print('Found:')
    for directory in argument.iterdir():
        ftir_dir = _search_directory(directory)
        if ftir_dir:
            config_file = selector.select_config(ftir_dir)
            print(f'- {ftir_dir}'
                  f'\n\t{config_file}')
            selector.adjust_data_path_location(config_file, ftir_dir)
            print()
    return