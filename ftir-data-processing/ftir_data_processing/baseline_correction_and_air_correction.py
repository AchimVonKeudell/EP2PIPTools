import os
from .file_extension import *
from .process_spectra import actual_molecule_name, obtain_and_correct_transmittance_spectra, reference_spectra


def get_air_contributions(output_directory: str, **metadata) -> dict:
    """

    :param output_directory:
    :param metadata:    sub-dictionary 'pre_treatment'/'air_contributions'
    :return:
    """

    air_contributions = {}

    for molecule_label, contribution_metadata in metadata.items():
        output_data_fn = f'{output_directory}\\{molecule_label}.{DATA_EXTENSION}'
        output_metadata_fn = f'{output_directory}\\{molecule_label}.{METADATA_EXTENSION}'
        output_figure_fn = f'{output_directory}\\{molecule_label}.{FIGURE_EXTENSION}'

        if os.path.exists(output_data_fn):
            air_contributions[molecule_label] = reference_spectra.ReferenceSpectrum.read_from_file(
                output_data_fn, output_metadata_fn
            )
            print(f'read {output_data_fn=}')
        else:
            method = contribution_metadata['method']
            if method == 'synthetic_spectrum':
                spectrum_data = reference_spectra.get_synthetic_spectrum(
                    molecule_name=actual_molecule_name(molecule_label),
                    **contribution_metadata
                )
            elif method == 'reference_spectrum':
                spectrum_data = reference_spectra.get_measured_spectrum(
                    **contribution_metadata
                )
            else:
                raise ValueError(f'{method=} not known/implemented')

            air_contributions[molecule_label] = spectrum_data
            spectrum_data.save_spectrum([output_data_fn, output_metadata_fn, output_figure_fn])


            print(f'acquired reference spectrum {molecule_label=} following {method=}')
    return air_contributions


def process_spectra(metadata: dict, processing_step: str):

    metadata_processing_step = metadata['processing_steps'][processing_step]
    metadata_conditions = metadata['conditions']

    ### 0.0. configure input & output directory path
    input_file_path = metadata['input data location']
    output_file_path = metadata['output data location']

    output_directory_reference_spectra = output_file_path + '\\reference_spectra'
    output_directory_method = output_file_path + '\\' + processing_step


    if not os.path.exists(input_file_path): # ensure that directory exists
        raise FileNotFoundError(f'{input_file_path=} not found')

    if not os.path.exists(output_directory_reference_spectra):  # ensure that directory exists
        os.makedirs(output_directory_reference_spectra)

    if not os.path.exists(output_directory_method):  # ensure that directory exists
        os.makedirs(output_directory_method)

    if not os.path.exists(output_directory_reference_spectra): # ensure that directory exists
        os.makedirs(output_directory_reference_spectra)

    print(f'Preprocessing:\n'
          f'\t{input_file_path=}\n'
          f'\t{output_directory_method=}\n'
          f'\t{output_directory_reference_spectra=}\n')


    ### 1. acquiring the reference gaseous absorption features
    air_contribution = get_air_contributions(
        output_directory=output_directory_reference_spectra,
        **metadata_processing_step['contributions']
    )

    def read_spectrum(input_file_name: str):

        full_input_file_name = f'{input_file_path}\\{input_file_name}'
        output_spectrum_file_names = [f'{output_directory_method}\\{input_file_name}.{ext}' for ext in (DATA_EXTENSION, FIGURE_EXTENSION)]

        if all([os.path.exists(fn) for fn in output_spectrum_file_names]):
            return # already done, not redoing work ;p

        try:
            obtain_and_correct_transmittance_spectra(
                input_data_file=full_input_file_name,
                output_data_file_name=output_spectrum_file_names[0],
                output_figure_file_name=output_spectrum_file_names[1],
                air_contribution=air_contribution,
                **metadata_processing_step
            )
        except PermissionError:
            print(f'PermissionError: {full_input_file_name=} is still opened in OPUS')
        except FileExistsError as e:
            print(e)
        return


    # 1.2) read the spectra and save them under ..\spectra
    for file in os.listdir(input_file_path):
        read_spectrum(file)
