from abc import ABC
from . import constants as const


class TransmissionMeasurementConditions(ABC):
    """
    Abstract class for representation of the measurement conditions for a transmission measurement. This includes the
    temperature and pressure of the measured gas mixture and the absorption path length
    """

    def __init__(self, pressure_atm, temperature_K, path_length_cm):
        self.pressure_atm = pressure_atm
        self.temperature_K = temperature_K
        self.path_length_cm = path_length_cm

    def get_number_density_cm3(self):
        """
        Get the number density in molecules per cm3. This assumes the pressure to be set in atmosphere and the 
        temperature in K
        :return: 
        """
        return self.pressure_atm * 101325 / (const.k_B * self.temperature_K) / 1e6


class TransmissionMeasurementSpeciesConditions(TransmissionMeasurementConditions):
    """
    Representation of the measurement conditions of a single species for a transmission measurement. In addition to the
    full gas mixture, this saves the total pressure, as well as the partial pressure of the species
    """
    def __init__(self, pressure_atm_total, pressure_atm_partial, temperature_K, path_length_cm):
        super().__init__(pressure_atm_partial, temperature_K, path_length_cm)
        self.pressure_atm_total = pressure_atm_total


class TransmissionMeasurementGasConditions(TransmissionMeasurementConditions):
    """
    Representation of the measurement conditions of the whole gas for a transmission measurement.
    """
    def get_species_conditions_from_mixing_ratio(self, mixing_ratio):
        """
        Get the measurement conditions for a species with a given mixing ratio. The returned object includes the
        temperature and path length, as well as total and partial pressures.
        :param mixing_ratio: The mixing ratio of a species, where 1 means that the total gas consists of molecules of 
        the species, and 0 means no molecules of the given species are present.
        :return: TransmissionMeasurementSpeciesConditions object
        """
        return TransmissionMeasurementSpeciesConditions(
            pressure_atm_total=self.pressure_atm,
            pressure_atm_partial=self.pressure_atm * mixing_ratio,
            temperature_K=self.temperature_K,
            path_length_cm=self.path_length_cm
        )
