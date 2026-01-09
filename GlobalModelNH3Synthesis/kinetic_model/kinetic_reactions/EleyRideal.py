from ..constants import boltzmann_constant, pi
from . import basic_reaction, species


class EleyRidealReaction(basic_reaction.Reaction):
    """Using Cantry's formula to calculate the rate for (dissociative) adsorption and Eley-Rideal kinetic_reactions.

        k = [Lambda^2 / D + V/A * 2 * (2 - gamma) / (v_av * gamma)] ^ -1, in [1/s]
            * diffusion term, only dominant when gamma = 1
            * sticking term, is dominant for gamma << 1
        where,
            Lambda  :   characteristic diffusion length    -> related to V/A, is d/2 for planar geometries
            D       :   diffusion constant      -> gas specie dependent
            V/A     :   volume to surface ratio
            gamma   :   sticking coefficient    -> only reaction dependent variable
            v_av    :   average speed of the particle towards the wall  -> gas specie dependent

        Adopted from [Carrasco et al. (2011) 10.1039/C1CP22284H, Shah et al. (2018) 10.1021/acsaem.8b00898].
    """

    sticking_coefficient: basic_reaction.Constant
    _gaseous_reactant_species = ''
    _reactant_specie: species.VolumeSpecie
    _background_gas: species.VolumeSpecie

    def set_constants(self, *args):
        if len(args) != 2:
            raise ValueError(f'Either not enough or too many constants in {args=}!\n')

        self.sticking_coefficient = basic_reaction.Constant(*args[0])
        for reactant_specie in self.reactants.labels:
            if 's' not in reactant_specie:
                self._gaseous_reactant_species = reactant_specie
                break
        if hasattr(species, self._gaseous_reactant_species):
            self._reactant_specie = getattr(species, self._gaseous_reactant_species)
        else:
            raise AttributeError(f'{self._gaseous_reactant_species=} not known/implemented')

        colliding_specie = args[1][0]
        if hasattr(species, colliding_specie):
            self._background_gas = getattr(species, colliding_specie)
        else:
            raise AttributeError(f'{colliding_specie=} not known/implemented')

    def _pre_factor_rate(self, number_of_species: int):
        factors = super()._pre_factor_rate(number_of_species)

        for species_index, species_label in zip(self.reactants.indexes + self.products.indexes,
                                                         self.reactants.labels + self.products.labels):
            if 's' in species_label:
                factors[species_index] *= self.volume_to_surface_coverage_ratio

        return factors

    @property
    def rate_equation(self) -> float:
        """Returns the adsorption rate, in [1/s]"""
        if self.sticking_coefficient.value == .0:
            return 0.
        average_velocity_m_s = (
            8 / pi * boltzmann_constant
            * self.conditions['gas_temperature_kelvin']
            / self._reactant_specie.mass_kg
        ) ** 0.5
        sticking_term = (
            2e-2 * self.conditions['volume_to_surface_area_ratio_cm']
            * (2 - self.sticking_coefficient.value)
            / (average_velocity_m_s * self.sticking_coefficient.value)
        )

        diffusion_m2_s = self._reactant_specie.diffusion_constant_m2_s(
            specie_b=self._background_gas,
            gas_temperature_kelvin=self.conditions['gas_temperature_kelvin'],
            pressure_atm=self.conditions['pressure_mbar'] / 1.01325e3
        )
        diffusion_term = 1e-4 * self.conditions['characteristic_diffusion_length_cm'] ** 2 / diffusion_m2_s
        return 1 / (sticking_term + diffusion_term)

