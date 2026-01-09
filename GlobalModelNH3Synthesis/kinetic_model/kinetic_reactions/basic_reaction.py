from dataclasses import dataclass
import numpy as np

boltzmann_constant = 1.38e-23       # J/K
torr_to_pa = 101325 / 760           # Pa / torr


@dataclass
class SpecieList:
    labels: list
    indexes: list
    coefficients: list

    def is_colliding_with_non_reacting_specie(self, reaction_unit: str) -> bool:
        """For gas-phase reactions, the unit is given in cm3/s or cm6/s.
        When the reaction is occurs with a 'third' party specie M,
        then the reaction rate is modified with the total density.

        :param reaction_unit:
        :return:
        """
        order = 0
        for i, specie_name in enumerate(self.labels):
            if 's' not in specie_name:
                order += self.coefficients[i]       # e.g. '2 H' counts as 2

        if ((reaction_unit == 'cm3/s') & (order == 1)) | (
                (reaction_unit == 'cm6/s') & (order == 2)):
            return True
        return False

    @staticmethod
    def set(species_reaction, species):
        return SpecieList(
            labels=list(species_reaction),
            indexes=[species.index(specie_name) for specie_name in species_reaction],
            coefficients=[coefficient for coefficient in species_reaction.values()],
        )


@dataclass
class Constant:
    value: float
    unit: str

    def __post_init__(self):
        self.value = float(self.value)


class Reaction:
    """
        k = C0

        C0: rate constants, in [1/torr/s], [cm3/s], [cm6/s]
    """

    rate_constant: Constant

    def __init__(self, raw_reaction: str, species: list, reactants: dict, products: dict, conditions: dict):
        self._raw_reaction = raw_reaction
        self.reactants = SpecieList.set(reactants, species)
        self.products = SpecieList.set(products, species)
        self.conditions = conditions

    def __repr__(self):
        return f'{self.__class__.__name__}({self._raw_reaction})'

    def set_constants(self, *args):
        if len(args) != 1:
            raise ValueError(f'Either not enough or too many constants in {args=}!\n')
        self.rate_constant = Constant(*args[0])
        return

    @property
    def volume_to_surface_coverage_ratio(self):
        """Returns volume / (surface area * surface site density) in cm^3.
        """
        return self.conditions['volume_to_surface_area_ratio_cm'] / self.conditions['surface_site_density_cm2']

    def _pre_factor_rate(self, number_of_species: int):
        """
        :param number_of_species:
        :return: array with all pre-factors for each species in the model
        """
        return np.ones((number_of_species, 1))

    @property
    def rate_equation(self) -> float:
        """Returns the rate of the reaction in [cm3/s] or [1/s].
        :return:
        """
        if self.rate_constant.unit == '1/s':
            return self.rate_constant.value
        elif self.reactants.is_colliding_with_non_reacting_specie(self.rate_constant.unit):
            return self.rate_constant.value * self.conditions['density_cm3']
        elif self.rate_constant.unit in ['cm3/s', 'cm6/s']:
            return self.rate_constant.value
        elif self.rate_constant.unit == '1/torr/s':
            conversion_to_cm3 = boltzmann_constant / torr_to_pa * self.conditions['gas_temperature_kelvin']
            return conversion_to_cm3 * self.rate_constant.value
        raise NotImplementedError(f'Unit of {self.rate_constant=} not implemented')

    def rate(self, specie_concentrations: list):

        # 1. determine if specie is lost or gained, and by how much, e.g. A + A -> A2 consumes two A species.
        phase = np.zeros((len(specie_concentrations), 1))
        for reactant_index, reactant_coefficients in zip(self.reactants.indexes, self.reactants.coefficients):
            phase[reactant_index] -= reactant_coefficients
        for product_index, product_coefficients in zip(self.products.indexes, self.products.coefficients):
            phase[product_index] += product_coefficients

        # 2. apply pre-factor, is needed when considering volume and surface reactions
        pre_factor = self._pre_factor_rate(len(specie_concentrations))

        # 3. calculate the reaction rate
        rate = self.rate_equation
        for reactant_index, reactant_coefficient in zip(self.reactants.indexes, self.reactants.coefficients):
            if reactant_coefficient == 1:
                rate *= specie_concentrations[reactant_index]
            else:
                rate *= specie_concentrations[reactant_index] ** reactant_coefficient
        return phase * pre_factor * rate
