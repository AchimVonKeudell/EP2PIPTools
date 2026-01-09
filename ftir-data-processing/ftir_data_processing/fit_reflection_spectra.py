import os
from .file_extension import *
from .process_spectra import reference_spectra, collect_results
from .vibrationalspectrum.fit import fit_reflection_spectrum


def process_spectra(metadata: dict, processing_step: str):

    metadata_processing_step = metadata['processing_steps'][processing_step]
    metadata_conditions = metadata['conditions']

    ### 0.0. configure input & output directory path
    input_file_path = metadata['input data location']
    output_file_path = metadata['output data location']
    output_directory_reference_spectra = output_file_path + '\\reference_spectra'
    output_directory_method = output_file_path + '\\' + processing_step

    # check if any OTHER steps involve `baseline_correction_and_air_correction`
    for key, value in metadata['processing_steps'].items():
        # if found, then changing directory of the input data accordingly
        if value['method'] == 'baseline_correction_and_air_correction':
            input_file_path = output_file_path + '\\' + key

    if not os.path.exists(input_file_path): # ensure that directory exists
        raise FileNotFoundError(f'{input_file_path=} not found')

    if not os.path.exists(output_directory_reference_spectra): # ensure that directory exists
        os.makedirs(output_directory_reference_spectra)

    if not os.path.exists(output_directory_method): # ensure that directory exists
        os.makedirs(output_directory_method)

    ### 1. acquiring the reference gaseous absorption features
    reference_spectra_gas_molecule = {}
    for gas_molecule_name, gas_molecule_metadata in metadata_processing_step['contributions_gas_phase'].items():
        reference_data_fn = f'{output_directory_reference_spectra}\\{gas_molecule_name}.{DATA_EXTENSION}'
        reference_metadata_fn = f'{output_directory_reference_spectra}\\{gas_molecule_name}.{METADATA_EXTENSION}'
        if os.path.exists(reference_data_fn):
            reference_spectra_gas_molecule[gas_molecule_name] = reference_spectra.ReferenceSpectrum.read_from_file(
                reference_data_fn, reference_metadata_fn
            )
            continue  # program is lazy, not going to redo work if it is not necessary
        reference_spectra_gas_molecule[gas_molecule_name] = reference_spectra.get_synthetic_spectrum(
            molecule_name=gas_molecule_name, **{**metadata_conditions, **gas_molecule_metadata}
        )
        reference_spectra_gas_molecule[gas_molecule_name].save_spectrum(
            file_names=[f'{output_directory_reference_spectra}\\{gas_molecule_name}.{extension}' for extension in EXTENSIONS]
        )

    ### 2. fitting the spectra
    def fit_spectrum(input_file_name: str):

        full_input_fn = input_file_path + '\\' + input_file_name

        file_name, _ = os.path.splitext(input_file_name)
        output_files = [f"{output_directory_method}\\{file_name}.{extension}" for extension in EXTENSIONS]

        if os.path.exists(output_files[0]):
            print(f'skip {output_files[0]}')
            return  # already done, not redoing work ;p

        fit_reflection_spectrum(
            spectrum_file=full_input_fn,
            output_files=output_files,
            reference_molecular_spectra=reference_spectra_gas_molecule,
            **metadata_conditions,
            **metadata_processing_step
        )
        print(f'done with {output_files[0]=}')
        return

    for file in os.listdir(input_file_path):
        if file[-4:] == '.dat':
            fit_spectrum(file)

    ### 3. collect outcome
    parameters = [f'fraction_{molecule_name}' for molecule_name in metadata_processing_step['contributions_gas_phase']]
    parameters += [f'amplitude_{band_name}' for band_name in metadata_processing_step['contributions_surface']]  # add surface contribution intensities
    for parameter in parameters:
        collect_results(
            directory_output_data=output_directory_method,
            output_directory=output_file_path + '\\',
            molecule_label=processing_step,
            parameter_name=parameter
        )
    return