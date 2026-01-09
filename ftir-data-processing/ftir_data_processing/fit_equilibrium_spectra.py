import os
from .file_extension import *
from .process_spectra import obtain_and_correct_transmittance_spectra, actual_molecule_name, collect_results
from .vibrationalspectrum.fit import fit_equilibrium_spectrum


def process_spectra(metadata: dict, processing_step: str):

    metadata_processing_step = metadata['processing_steps'][processing_step]
    metadata_conditions = metadata['conditions']

    ### 0. configure input & output directory path
    input_file_path = metadata['input data location']
    output_file_path = metadata['output data location']

    output_directory_spectra = f'{output_file_path}\\{processing_step}\\spectra'
    output_directory_fits = f'{output_file_path}\\{processing_step}\\fits'

    if not os.path.exists(output_directory_spectra): # ensure that directory exists
        os.makedirs(output_directory_spectra)
    if not os.path.exists(output_directory_fits): # ensure that directory exists
        os.makedirs(output_directory_fits)


    ### 1. collect the spectra
    def read_spectrum(input_file_name: str):

        full_input_fn = f'{input_file_path}/{input_file_name}'

        output_files = [f'{output_directory_spectra}\\{input_file_name}.{file_type}' for file_type in [DATA_EXTENSION, FIGURE_EXTENSION]]
        if os.path.exists(output_files[0]):
            print(f'skip {input_file_name}')
            return  # already done, not redoing work ;p
        try:
            obtain_and_correct_transmittance_spectra(
                input_data_file=full_input_fn,
                output_data_file_name=output_files[0],
                output_figure_file_name=output_files[1],
                **metadata_processing_step
            )
        except PermissionError:
            print(f'PermissionError: {input_file_name=} is still opened in OPUS')
        except FileExistsError as e:
            print(e)
        return

    for file in os.listdir(input_file_path):
        read_spectrum(file)


    ### 2. fit the spectra
    def fit_spectrum(input_file_name: str):

        full_input_fn = f'{output_directory_spectra}\\{input_file_name}'

        fn_index, _ = os.path.splitext(input_file_name)
        output_files = [f"{output_directory_fits}\\{fn_index}.{extension}" for extension in EXTENSIONS]

        if os.path.exists(output_files[0]):
            print(f'skip {output_files[0]}')
            return  # already done, not redoing work ;p

        fit_equilibrium_spectrum(
            spectrum_file=full_input_fn,
            output_files=output_files,
            molecule_name=actual_molecule_name(processing_step),
            **metadata_processing_step,
            **metadata_conditions  # parameters of metadata_conditions will overwrite that of metadata_processing_step
        )
        return

    for file in os.listdir(output_directory_spectra):
        if file[-4:] == '.' + DATA_EXTENSION:
            fit_spectrum(file)


    ### 3. collect outcome
    for parameter in [f'fraction_{actual_molecule_name(processing_step).lower()}',
                      'temperature_rotation_kelvin'
                      ]:
        if parameter + '_range' in metadata_processing_step:
            collect_results(
                directory_output_data=output_directory_fits,
                output_directory=output_file_path + '\\',
                molecule_label=processing_step,
                parameter_name=parameter
            )
    return