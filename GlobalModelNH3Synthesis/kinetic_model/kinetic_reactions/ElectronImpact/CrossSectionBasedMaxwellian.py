from dataclasses import dataclass
import numpy as np
import os.path
from ...constants import electron_volt_to_speed_squared, NORM_FACTOR_MAXWELLIAN_EEDF, cross_sections_path
from .. import basic_reaction, ElectronImpact


@dataclass
class CrossSectionSet:

    energies_ev: np.ndarray
    cross_sections_m2: np.ndarray

    energy_steps_ev: np.ndarray = None

    def __post_init__(self):
        """Better approximating the energy step for each point in the EEDF.
        To create a histogram where the rectangles are centered at the `energy_ev` values,
            and the width is midway between the values around it.
        For instance, for point i is centered at E_i and
            the energy step on the left is [E_i-E_(i-1)] / 2
            the energy step on the right is [E_(i+1)-E_i] / 2
        """
        diff_energies = np.diff(self.energies_ev)
        self.energy_steps_ev = np.append(0, diff_energies) + np.append(diff_energies, 0)
        re_correct_for_double_delta_energy_values = np.ones_like(self.energies_ev)
        re_correct_for_double_delta_energy_values[1:-1] *= 0.5
        self.energy_steps_ev *= re_correct_for_double_delta_energy_values

    @staticmethod
    def from_file(file_name: str):
        _data = np.loadtxt(file_name, skiprows=1)
        return CrossSectionSet(energies_ev=_data[:, 0], cross_sections_m2=_data[:, 1])


class ElectronImpactMaxwellianEEDF(ElectronImpact.ElectronImpactReaction):
    """
        k = int{cross_section * v_e * f_EEPF(E)}
        where,
            v_e = sqrt{2 * e * E / m_e}
            f_EEPF(E) = 2 * (E/pi) ^ 0.5 * f_EEDF(E)
            f_EEDF(E) = E ^ -1.5 * exp(-E/T_e)
        with,
            v_e:    electron speed [m/s]
            E:      electron energy [eV]
            T_e:    electron temperature [eV]
            e:      electron charge [C]
            m_e:    mass electron [kg]
    """

    normalisation_factor = NORM_FACTOR_MAXWELLIAN_EEDF
    cross_section_set: CrossSectionSet

    def set_constants(self, *args):
        """
            Only one variable: cross-section file
        """
        if len(args) != 1:
            raise ValueError(f'Either not enough or too many constants in {args=}!\n')

        self.rate_constant = basic_reaction.Constant(0, 'cm3/s')

        full_file_path = f"{cross_sections_path}\\{args[0][0]}"
        if not os.path.exists(full_file_path):
            raise FileNotFoundError(f"{full_file_path=} not found")
        self.cross_section_set = CrossSectionSet.from_file(full_file_path)

    def electron_energy_density_function(self, electron_energy_ev: np.ndarray):
        return electron_energy_ev ** -1.5 * np.exp(
            -electron_energy_ev / self.conditions['electron_temperature_ev']
        )

    def electron_energy_probability_function(self, electron_energy_ev: np.ndarray):
        return (self.normalisation_factor * np.sqrt(electron_energy_ev)
                * self.electron_energy_density_function(electron_energy_ev))

    @property
    def rate_equation(self) -> float:
        electron_speed_m_s = np.sqrt(electron_volt_to_speed_squared * self.cross_section_set.energies_ev)
        self.rate_constant.value = 1e6 * np.sum(
            self.cross_section_set.cross_sections_m2 * electron_speed_m_s
            * self.electron_energy_probability_function(self.cross_section_set.energies_ev)
            * self.cross_section_set.energy_steps_ev
        )
        return super().rate_equation
