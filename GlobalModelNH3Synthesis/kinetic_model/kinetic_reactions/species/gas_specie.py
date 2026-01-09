from dataclasses import dataclass
from ...constants import avogadros_constant


@dataclass
class VolumeSpecie:
    mass_g_per_mol: float
    diffusion_volume: float

    @property
    def mass_kg(self):
        return 1e-3 * self.mass_g_per_mol / avogadros_constant

    def diffusion_constant_m2_s(self, specie_b, gas_temperature_kelvin=300., pressure_atm=1.0):
        """Using the function described in Fuller et al (1966) https://doi.org/10.1021/ie50677a007
        When considering a binary gas phase system of particle `self` and specie_b.

        The factor 1e-7 originates from the 1e-3 found in the reference and 1e-4 to convert the unit from cm2/s to m2/s.

        :param specie_b:
        :param gas_temperature_kelvin:
        :param pressure_atm:
        :return:
        """
        _mass_factor = (1 / self.mass_g_per_mol + 1 / specie_b.mass_g_per_mol) ** 0.5
        _collision_factor = (self.diffusion_volume ** 0.333 + specie_b.diffusion_volume ** 0.333) ** 2
        return 1e-7 * gas_temperature_kelvin ** 1.75 * _mass_factor / pressure_atm / _collision_factor

