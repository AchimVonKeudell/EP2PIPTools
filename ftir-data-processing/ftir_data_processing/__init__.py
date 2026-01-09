import os
import shutil
import importlib
from . import yaml_handler


VALID_ROUTINES = {
    'baseline_correction_and_air_correction',
    'fit_background_spectra',
    'fit_equilibrium_spectra',
    'fit_non_equilibrium_spectra',
    'fit_reflection_spectra',
}


def run(metadata_file_path: str):

    ## 0. get metadata
    print(f'{metadata_file_path=}')
    metadata = yaml_handler.read_yaml(metadata_file_path)

    ### 0.1. check if the data exists
    raw_data_path = metadata['input data location']
    if not os.path.exists(raw_data_path):
        raise FileExistsError(f'{raw_data_path=} not found')

    ## 1. save the settings
    output_data_path = metadata['output data location']
    if not os.path.exists(output_data_path):
        os.makedirs(output_data_path)

    metadata_file_output = output_data_path + '\\' + 'config.yaml'
    try:
        shutil.copy(metadata_file_path, metadata_file_output)
    except shutil.SameFileError:
        pass

    ## 2. acquire the timing values for the spectra
    try:
        from .process_spectra.collect_results import create_timing_dat
        create_timing_dat(raw_data_path, output_data_path, **metadata)
    except IndexError:
        pass # If timing extraction fails (e.g. background file missing), just continue

    ## 3. follow the routines set by each specie
    for processing_step, metadata_specie in metadata['processing_steps'].items():
        method = metadata_specie['method']
        print(f'{processing_step=}, {method=}')

        if method in VALID_ROUTINES:
            module = importlib.import_module(f'.{method}', package=__package__)
            module.process_spectra(metadata, processing_step)
        else:
            raise ValueError(f'{method=} not known/implemented\n'
                             f'try one of {VALID_ROUTINES=}')

    return