import os
import numpy as np
import brukeropusreader
from . import plot, reference_spectra

OPUS_REFERENCE_SPECTRUM_KEY = 'ScRf'
OPUS_SIGNAL_SPECTRUM_KEY = 'ScSm'
OPUS_TRANSMITTANCE_KEY = 'AB'



def get_snip(data_array: np.ndarray, lower_x_limit, upper_x_limit):
    """Function that selects a given part of the data_array

    :param data_array: array with size [2,N];
    :param lower_x_limit:
    :param upper_x_limit:
    :return:
    """
    mask = (data_array[:, 0] > lower_x_limit) & (data_array[:, 0] < upper_x_limit)
    return data_array[mask].copy()


def get_measured_spectrum(filepath: str, **kwargs) -> np.ndarray:
    """Reads the measured spectrum from the given filepath.

    :param filepath:
    :param kwargs:
    :return:
    """
    def read_from_opus(raw_opus_data, data_key):
        _x = raw_opus_data.get_range(data_key)
        _y = raw_opus_data[data_key]

        if _x[0] > _x[-1]: # reverse if order is incorrect
            _x ,_y = _x[::-1], _y[::-1]
        return np.stack([_x, _y], axis=1)

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{filepath=} not found")

    if filepath[-4:] == '.dat': # filepath refers to a .dat file that already contains the spectrum
        return np.loadtxt(filepath)[:, :2]

    ## acquiring the spectrum from OPUS file:
    opus_data = brukeropusreader.read_file(filepath)

    if OPUS_TRANSMITTANCE_KEY not in opus_data:
        raise FileExistsError(f'{filepath=} does not contain {OPUS_TRANSMITTANCE_KEY=}, so skipping it.'
                              f'\n\tIt is probably a reference spectrum.')

    if 'background_spectrum_file_name' in kwargs:
        signal_spectrum = read_from_opus(opus_data, OPUS_SIGNAL_SPECTRUM_KEY)

        # read background spectrum
        filepath_background_spectrum = kwargs['background_spectrum_file_name']
        if filepath_background_spectrum is None:
            raise ValueError(f'{filepath_background_spectrum=}, why do you do that to me? '
                             f'Either define a file, or remove the parameter')
        background_opus_data = brukeropusreader.read_file(filepath_background_spectrum)

        background_spectrum = read_from_opus(
            background_opus_data,
            OPUS_SIGNAL_SPECTRUM_KEY if OPUS_SIGNAL_SPECTRUM_KEY in background_opus_data else OPUS_REFERENCE_SPECTRUM_KEY
        )
        signal_spectrum[:, 1] /= np.interp(signal_spectrum[:, 0], background_spectrum[:, 0], background_spectrum[:, 1])
        return signal_spectrum

    transmittance_spectrum = read_from_opus(opus_data, OPUS_TRANSMITTANCE_KEY)[:-1]
    return transmittance_spectrum


def get_baseline_coefficients(spectrum_data: np.ndarray, baseline_polynomial_order: int, **kwargs):
    """Estimates the baseline from the given parameters:
        if wavenumber_range_baseline_correction is given:
            the baseline is estimated from the given wavenumber range
        if baseline_x_points is given:
            the baseline is estimated from the given waveumber points

    :param spectrum_data:
    :param baseline_polynomial_order:
    :param kwargs:
    :return:
    """

    if 'wavenumber_range_baseline_correction' in kwargs:
        baseline_spectrum_data = get_snip(spectrum_data, *kwargs['wavenumber_range_baseline_correction'])
        baseline_x_points = baseline_spectrum_data[:, 0]
        baseline_y_points = baseline_spectrum_data[:, 1]
    elif 'baseline_x_points' in kwargs:
        baseline_x_points = kwargs['baseline_x_points']
        baseline_y_points = np.interp(baseline_x_points, *spectrum_data.T)
    else:
        print('No baseline correction')
        return [1]

    baseline_coefficients = np.polyfit(baseline_x_points, baseline_y_points, baseline_polynomial_order)
    return baseline_coefficients


def obtain_and_correct_transmittance_spectra(
        input_data_file,
        output_data_file_name,
        output_figure_file_name,
        baseline_polynomial_order,
        **kwargs
):
    """Processes a sample FTIR spectrum by optionally:
        * applying background (baseline) correction,
        * correcting for air contribution
    It saves the corrected spectrum in a data file and a figure.

    :param input_data_file: Path to the .dat or OPUS file, containing spectrum.
    :param output_data_file_name: full file-name to save the corrected spectrum as a .txt or .dat file,
    :param output_figure_file_name: full file-name to save the diagnostic plot as a .png or .pdf.
    :param baseline_polynomial_order: an integer
    :param kwargs: May contain the following parameters:
        - wavenumber_range_spectrum
        - background_spectrum_file_name
        - wavenumber_range_baseline_correction or baseline_x_points
        - air_contribution
    :return:
    """

    # 1) obtaining the data
    experimental_data = get_measured_spectrum(input_data_file, **kwargs)

    # 2. determine baseline from a known 'emtpy' part of the spectrum
    baseline_coefficients = get_baseline_coefficients(experimental_data, baseline_polynomial_order, **kwargs)

    # 3. selecting only the part that we want (if range is given)
    if 'wavenumber_range_spectrum' in kwargs:
        wavenumber_range_spectrum = kwargs['wavenumber_range_spectrum']
        experimental_data =  get_snip(experimental_data, *wavenumber_range_spectrum)

    # 4. apply baseline correction
    baseline_spectrum = np.polyval(baseline_coefficients, experimental_data[:, 0])
    corrected_experimental_data = experimental_data.copy()
    corrected_experimental_data[:, 1] /= baseline_spectrum

    # 3. remove contributions from an open air path (if `air_contributions` is given)
    if 'air_contribution' in kwargs:
        best_fit_air_contributions = {}
        air_contribution = kwargs['air_contribution']
        for label, artificial_spectrum in air_contribution.items():
            best_fit_air_contributions[label] = reference_spectra.fit_contribution(
                measured_spectrum=corrected_experimental_data,
                reference_spectrum=artificial_spectrum
            )
            corrected_experimental_data[:, 1] /= best_fit_air_contributions[label]

        kwargs['best_fit_air_contributions'] = best_fit_air_contributions

    # 4. save the spectrum
    np.savetxt(output_data_file_name, corrected_experimental_data,
               delimiter='\t', header=f"{baseline_coefficients=}\nwavenumber/cm-1\tdata")

    # 5. plot the measured spectrum with the baseline correction
    plot.plot_transmittance_spectrum(
        filepath_figure=output_figure_file_name,
        wavenumbers=experimental_data[:, 0],
        original_spectrum=experimental_data[:, 1],
        corrected_spectrum=corrected_experimental_data[:, 1],
        baseline=baseline_spectrum,
        **kwargs
    )
    return





