import os

from .process_spectra import obtain_and_correct_transmittance_spectra, actual_molecule_name
from .vibrationalspectrum.fit import fit_equilibrium_spectrum


def process_spectra(metadata: dict, molecule_name: str):

    input_file_path = metadata['input data location']
    output_file_path = metadata['output data location']
    molecule_metadata = metadata['routines'][molecule_name]
    conditions_metadata = metadata['conditions']

    output_paths = [f'{output_file_path}\\{molecule_name}\\{key}' for key in ['spectra', 'fits']]

    if not os.path.exists(output_paths[0]):
        os.makedirs(output_paths[0])
    if not os.path.exists(output_paths[1]):
        os.makedirs(output_paths[1])


    def read_spectrum(input_file_name: str):

        full_input_fn = f'{input_file_path}/{file}'

        output_file_name = file[len(file) - file[::-1].find('.'):]

        output_files = [f'{output_paths[0]}/{output_file_name}.{file_type}' for file_type in ['txt', 'png']]
        if os.path.exists(output_files[0]):
            print(f'skip {file}')
            return  # already done, not redoing work ;p
        try:
            obtain_and_correct_transmittance_spectra(
                full_input_fn, output_files,
                molecule_metadata['fit wavenumber range'],
                wavenumber_range_baseline_correction=molecule_metadata['baseline wavenumber range'],
                baseline_polynomial_order=molecule_metadata['baseline polynomial order'],
                datetype='ScRf'
            )
        except PermissionError:
            print(f'PermissionError: {input_file_name=} is still opened in OPUS')
        except FileExistsError as e:
            print(e)
        return


    def fit_spectrum(input_file_name: str):

        full_input_fn = f'{output_paths[0]}/{input_file_name}'

        fn_index = input_file_name.split('.')[0]
        output_files = [f"{output_paths[1]}\\{fn_index}.{file_type}" for file_type in ['metadata', 'txt', 'png']]

        if os.path.exists(output_files[0]):
            print(f'skip {output_files[0]}')
            return  # already done, not redoing work ;p

        kwargs = {'p0': 0., **molecule_metadata}

        if 'instrumental broadening' in molecule_metadata:
            kwargs['instrumental_broadening_sigma'] = molecule_metadata['instrumental broadening']
        if 'wavenumber offset' in molecule_metadata:
            kwargs['temperature_rotation_kelvin'] = molecule_metadata['wavenumber offset']

        if 'gas temperature [K]' in conditions_metadata:
            kwargs['temperature_rotation_kelvin'] = conditions_metadata['gas temperature [K]']


        fit_equilibrium_spectrum(
            spectrum_file=full_input_fn,
            output_files=output_files,
            database_file=molecule_metadata['database filepath'],
            pressure_atm=conditions_metadata['pressure [atm]'],
            path_length_cm=conditions_metadata['path length [cm]'],
            molecule_name=actual_molecule_name(molecule_name),
            **kwargs
        )
        return


    for file in os.listdir(input_file_path):
        read_spectrum(file)

    for file in os.listdir(output_paths[0]):
        if file[-4:] == '.txt':
            fit_spectrum(file)

