import numpy as np

from .parse import acquire_species_from_reaction_list, parse_reaction
from .configure_reactions import configure_reaction


class KineticModel:

    def __init__(self, reactions_raw_format: list, conditions: dict):
        """
        :param reactions_raw_format:
        :param conditions:
        return;
        """
        reactions_data = []
        for _reaction in reactions_raw_format:
            if '#' in _reaction:
                continue
            if len(_reaction.split(',')) > 0:
                reactions_data.append(_reaction.split(','))

        # collect the species
        self.species = acquire_species_from_reaction_list(reactions_data)
        print(f'Found species: {self.species}')

        # configure the reactions
        self.reactions = [configure_reaction(_reaction, self.species, conditions) for _reaction in reactions_data]

        # # incorporate the volume-surface correction factors
        # if 'surface_to_volume_ratio' in kwargs:
        #     _surface_to_volume_ratio = kwargs['surface_to_volume_ratio']
        # else:
        #     _surface_to_volume_ratio = 1
        # self.surface_volume_correction_factors = self.get_correction_factor(_surface_to_volume_ratio)

    def get_correction_factor(self, surface_to_volume_ratio: float) -> np.ndarray:
        return np.ones((len(self.species), 1))

    def ode(self, t, y) -> np.ndarray:
        dy_dt = np.zeros((len(self.species), 1))
        for reaction in self.reactions:
            dy_dt += reaction.rate(y)
        # dy_dt *= self.surface_volume_correction_factors
        return dy_dt.T[0]


class SurfaceKineticModel(KineticModel):
    """
        If gas-surface reactions are in
            * [1/torr/s], then use molar fractions [-],
            * [cm3/s], then use concentrations [1/cm3];

    """

    def get_correction_factor(self, surface_to_volume_ratio: float) -> np.ndarray:
        factors = np.ones((len(self.species), 1))
        for index, specie in enumerate(self.species):
            if 's' not in specie:
                factors[index] *= surface_to_volume_ratio
        return factors


class PlasmaKineticModel(KineticModel):
    """
        Unit of gaseous species [1/cm3];
    """
    pass


__all__ = ['SurfaceKineticModel', 'PlasmaKineticModel']
