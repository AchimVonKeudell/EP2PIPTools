import numpy as np
from ..constants import constant_kJmol_to_kelvin, constant_eV_to_kelvin
from . import basic_reaction


class ArrheniusReaction(basic_reaction.Reaction):
    """
        k = C0 * exp(-C1 / T_gas),

        C0: pre-exponential term, in [1/torr/s], [cm3/s], [cm6/s]
        C1: energy barrier, in [K]
    """

    activation_energy_barrier: basic_reaction.Constant

    def set_constants(self, *args):
        if len(args) != 2:
            raise ValueError(f'Either not enough or too many constants in {args=}!')
        super().set_constants(args[0])
        self.activation_energy_barrier = basic_reaction.Constant(*args[1])

    @property
    def _exponent_term(self) -> float:
        return self.activation_energy_barrier.value / self.conditions['gas_temperature_kelvin']

    @property
    def rate_equation(self) -> float:
        rate = super().rate_equation

        if self.activation_energy_barrier.value == 0.:
            return rate        # exp(0) = 1, so no calculations required
        rate *= np.exp(-self._exponent_term)
        return rate


class ArrheniusSurfaceReaction(ArrheniusReaction):
    """
        k = C0 * exp(-C1 / T_surface),

        C0: pre-exponential term, in [1/torr/s], ([cm3/s], [cm6/s] is weird but okay)
        C1: energy barrier, in [kJ/mol]
    """

    @property
    def _exponent_term(self) -> float:
        if self.activation_energy_barrier.unit == 'kJ/mol':
            energy_term = constant_kJmol_to_kelvin * self.activation_energy_barrier.value
        elif self.activation_energy_barrier.unit == 'eV':
            energy_term = constant_eV_to_kelvin * self.activation_energy_barrier.value
        else:
            raise NotImplementedError(f'{self.activation_energy_barrier.value=} not implemented')

        return energy_term / self.conditions['surface_temperature_kelvin']


class ExtendedArrheniusReaction(ArrheniusReaction):
    """
        k = C0 * (T / C1) ** C2 * exp(-C3/T)

        C0: pre-exponential term, in [1/torr/s], [cm3/s], [cm6/s]
        C1: reference temperature, in [K]
        C2: temperature exponent, in [-]
        C3: energy barrier, in [K]
    """

    reference_temperature: basic_reaction.Constant
    temperature_exponent: basic_reaction.Constant

    def set_constants(self, *args):
        if len(args) != 4:
            raise ValueError(f'Either not enough or too many constants {args=}')
        super().set_constants(args[0], args[3])

        self.reference_temperature = basic_reaction.Constant(*args[1])
        if self.reference_temperature.unit != 'K':
            raise ValueError(f'Unit of {self.reference_temperature=} does not make sense to me :(')
        self.temperature_exponent = basic_reaction.Constant(*args[2])

    @property
    def rate_equation(self):
        rate = super().rate_equation
        rate *= (self.conditions['gas_temperature_kelvin'] / self.reference_temperature.value
                 ) ** self.temperature_exponent.value
        return rate


class LangmuirHinshelwoodReaction(ArrheniusSurfaceReaction):
    """From [Hong et al. (2017) 10.1088/1361-6463/aa6229]:
        k = nu / 4 * exp[-(E_d + E_a) / (kB * T_surface)]
        with,
            nu          :   surface diffusional jump frequency, which is 4.3e12 [1/s] for N atoms on Fe(001)
            E_d         :   diffusion energy barrier    -> surface dependent
            E_a         :   activation energy barrier   -> reaction dependent
            kB          :   boltzmann's constant ~ 1.38e-23 [J/K]
            T_surface   :   surface temperature [K]
    """

    def _pre_factor_rate(self, number_of_species: int):
        factors = super()._pre_factor_rate(number_of_species)

        for species_index, species_label in zip(self.reactants.indexes + self.products.indexes,
                                                         self.reactants.labels + self.products.labels):
            if 's' in species_label:
                pass
            else:
                factors[species_index] /= self.volume_to_surface_coverage_ratio

        return factors

    def set_constants(self, *args):
        if len(args) != 3:
            raise ValueError(f'Either not enough or too many constants in {args=}!')
        surface_diffusional_jump_frequency = basic_reaction.Constant(*args[0])
        if surface_diffusional_jump_frequency.unit == '1/s':
            self.rate_constant = basic_reaction.Constant(
                value=surface_diffusional_jump_frequency.value / 4,
                unit=surface_diffusional_jump_frequency.unit
            )
        else:
            raise NotImplementedError(f'{surface_diffusional_jump_frequency.unit=} not implemented')
        self.activation_energy_barrier = basic_reaction.Constant(*args[1])
        diffusion_energy_barrier = basic_reaction.Constant(*args[2])
        self.activation_energy_barrier.value += diffusion_energy_barrier.value
        return
