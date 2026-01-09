from dataclasses import dataclass
import os.path
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from .. import yaml_handler
from ..vibrationalspectrum import GetSimulator
from ..vibrationalspectrum.spectrumsimulator import spectrum, distributions, conditions, constants


@dataclass
class ReferenceSpectrum:

    path_length_cm: float
    pressure_mbar: float
    temperature_kelvin: float
    molar_fraction: float
    instrumental_broadening_cm: float

    file_name: str

    wavenumber_array: np.ndarray
    data_array: np.ndarray

    def __repr__(self):
        return f'{self.__class__.__name__}({self.file_name})'

    def __add__(self, other):
        if not isinstance(other, ReferenceSpectrum):
            raise ValueError(f'{other=} should be {self.__class__.__name__} object')

        parameters = {}
        for key, value in self.__dict__.items():
            if '_array' not in key:
                parameters[key] = value

        return ReferenceSpectrum(
            **parameters,
            wavenumber_array=np.append(self.wavenumber_array, other.wavenumber_array),
            data_array=np.append(self.data_array, other.data_array)
        )

    @staticmethod
    def read_from_file(data_array_file_name, reference_spectrum_metadata_file_name):
        """

        :param data_array_file_name: with DATA_EXTENSION
        :param reference_spectrum_metadata_file_name: with METADATA_EXTENSION
        :return:
        """
        _data_array = np.loadtxt(data_array_file_name)

        metadata = yaml_handler.read_yaml(reference_spectrum_metadata_file_name)
        pressure_mbar = metadata['pressure_mbar']
        temperature_kelvin = metadata['temperature_kelvin']
        path_length_cm = metadata['path_length_cm']
        molar_fraction = metadata['molar_fraction']
        instrumental_broadening_cm = metadata['instrumental_broadening_cm']


        return ReferenceSpectrum(
            pressure_mbar=pressure_mbar,
            temperature_kelvin=temperature_kelvin,
            path_length_cm=path_length_cm,
            molar_fraction=molar_fraction,
            instrumental_broadening_cm=instrumental_broadening_cm,
            file_name=data_array_file_name,
            wavenumber_array=_data_array[:, 0],
            data_array=_data_array[:, 1],
        )

    @property
    def reference_number_density_m3(self):
        return self.molar_fraction * (self.pressure_mbar * 1e2) / (constants.k_B * self.temperature_kelvin)

    @property
    def reference_number_density_cm3(self):
        return self.reference_number_density_m3 * 1e-6

    def interpolate(self, x, is_absorbance=False):
        if is_absorbance:
            return -np.log(np.interp(x, self.wavenumber_array, self.data_array, left=1, right=1))
        return np.interp(x, self.wavenumber_array, self.data_array, left=1, right=1)


    def save_spectrum(self, file_names: list, **other_parameters):
        """

        :param file_names: a list the file name ending the extensions following file_extension.EXTENSIONS
        :param other_parameters:
        :return:
        """
        self.file_name = file_names[0]

        np.savetxt(
            fname=self.file_name,
            X=np.stack([self.wavenumber_array, self.data_array], axis=1),
            delimiter='\t',
            header=f'metadata file name = {file_names}'
        )

        yaml_handler.write_yaml(
            filepath=file_names[1],
            metadata={'pressure_mbar': self.pressure_mbar,
                      'temperature_kelvin': self.temperature_kelvin,
                      'path_length_cm': self.path_length_cm,
                      'molar_fraction': self.molar_fraction,
                      'instrumental_broadening_cm': self.instrumental_broadening_cm,
                      'data file name': file_names[0],
                      **other_parameters}
        )

        plt.rc('axes', labelweight='bold')
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set(xlabel='wavenumber / cm$^{-1}$',
               ylabel='transmittance / -')
        ax.plot(self.wavenumber_array, self.data_array)
        fig.tight_layout()
        fig.savefig(file_names[2])
        plt.close()
        return


def get_synthetic_spectrum(
        molecule_name: str,
        hitran_database_file: str,
        pressure_atm: float,
        molar_fraction_reference: float,
        gas_temperature_kelvin: float,
        path_length_cm: float,
        instrumental_broadening_sigma: float,
        absorbance_wavenumber_bands: list,
        wavenumber_resolution: float,
        **kwargs
) -> ReferenceSpectrum:

    if absorbance_wavenumber_bands and isinstance(absorbance_wavenumber_bands[0], (list, tuple)):
        # this part deals with multiple absorbance bands, the function can be recycled, so another function does not
        # need to be written :)
        output = None
        for wavenumber_band in absorbance_wavenumber_bands:
            reference_spectrum = get_synthetic_spectrum(
                molecule_name=molecule_name,
                hitran_database_file=hitran_database_file,
                pressure_atm=pressure_atm,
                molar_fraction_reference=molar_fraction_reference,
                gas_temperature_kelvin=gas_temperature_kelvin,
                path_length_cm=path_length_cm,
                instrumental_broadening_sigma=instrumental_broadening_sigma,
                absorbance_wavenumber_bands=wavenumber_band,  # Pass the individual band
                wavenumber_resolution=wavenumber_resolution,
                **kwargs
            )
            if output is None:
                output = reference_spectrum
            else:
                output += reference_spectrum
        return output

    if absorbance_wavenumber_bands[0] > absorbance_wavenumber_bands[1]:
        raise ValueError(f'wavenumber range: {absorbance_wavenumber_bands=} does not make sense, [lower_limit, upper_limit]')

    wavenumber_array = np.arange(absorbance_wavenumber_bands[0], absorbance_wavenumber_bands[1] + wavenumber_resolution,
                                   wavenumber_resolution)
    spectrum_data = spectrum.Spectrum(wavenumber_array=wavenumber_array, data_array=np.ones_like(wavenumber_array))


    molecule_simulator = getattr(GetSimulator, f'get_{molecule_name.upper()}_simulator')(
        hitran_database_file=hitran_database_file,
        wavenumber_start=min(absorbance_wavenumber_bands),
        wavenumber_end=max(absorbance_wavenumber_bands)
    )

    distribution = distributions.EquilibriumDistribution()
    distribution.set_temperature_K(gas_temperature_kelvin)
    molecule_simulator.set_distribution(distribution)

    molecule_simulator.set_instrumental_broadening_coefficient(instrumental_broadening_sigma)

    condition = conditions.TransmissionMeasurementGasConditions(
        pressure_atm=pressure_atm,
        temperature_K=gas_temperature_kelvin,
        path_length_cm=path_length_cm
    )

    spectrum_data.data_array *= molecule_simulator.add_transmission_to_spectrum(
        spectrum=spectrum_data,
        conditions=condition.get_species_conditions_from_mixing_ratio(molar_fraction_reference),
        inclusion_strength_threshold=1e-25, simulation_transmission_threshold=1-1e-6
    )
    molecule_simulator.apply_instrumental_broadening(spectrum_data)
    return ReferenceSpectrum(
        pressure_mbar=pressure_atm * 1013.25,
        temperature_kelvin=gas_temperature_kelvin,
        path_length_cm=path_length_cm,
        molar_fraction=molar_fraction_reference,
        instrumental_broadening_cm=instrumental_broadening_sigma,
        file_name='',
        wavenumber_array=spectrum_data.wavenumber_array,
        data_array=spectrum_data.data_array,
    )


def get_measured_spectrum(
        reference_spectrum_file_name: str,
        absorbance_wavenumber_bands: list,
        **reference_conditions
) -> ReferenceSpectrum:
    if not os.path.exists(reference_spectrum_file_name):
        raise ValueError(f'{reference_spectrum_file_name=} does not exist')

    reference_spectrum = np.loadtxt(reference_spectrum_file_name)
    mask = (reference_spectrum[:, 0] > absorbance_wavenumber_bands[0]) & (
            reference_spectrum[:, 0] < absorbance_wavenumber_bands[1])
    reference_spectrum = reference_spectrum[mask]
    return ReferenceSpectrum(
        path_length_cm=None,
        pressure_mbar=None,
        temperature_kelvin=None,
        molar_fraction=1.0,
        instrumental_broadening_cm=None,
        file_name='',
        wavenumber_array=reference_spectrum[:, 0],
        data_array=reference_spectrum[:, 1],
    )


def fit_contribution(measured_spectrum: np.ndarray, reference_spectrum: ReferenceSpectrum):
    """Fitting the contribution by fitting the absorbance instead of the transmittance.

    :param measured_spectrum:
    :param reference_spectrum:
    :return:
    """
    x_data = measured_spectrum[:, 0]
    y_data = measured_spectrum[:, 1]

    mask = (x_data > reference_spectrum.wavenumber_array.min()) & (
            x_data < reference_spectrum.wavenumber_array.max())

    output_fit = curve_fit(
        lambda x, ratio: np.diff(reference_spectrum.interpolate(x) ** ratio),
        x_data[mask], np.diff(y_data[mask]), p0=[1.], bounds=[[-2], [2]]
    )
    best_fit_ratio = float(output_fit[0])
    best_fit_spectrum = reference_spectrum.interpolate(measured_spectrum[:, 0]) ** best_fit_ratio
    return best_fit_spectrum


