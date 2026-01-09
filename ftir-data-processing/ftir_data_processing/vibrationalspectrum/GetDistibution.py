from .spectrumsimulator.distributions import EquilibriumDistribution, \
    TreanorTriatomicLinearFermiResonantFractionalDistribution, TreanorDiatomicFractionalDistribution, \
    NonEquilibriumBoltzmannPyramidalTetratomicFractionDistribution
from .spectrumsimulator.constants import *


def get_equilibrium_distributions(number_of_distributions=1):
    """Returns a list of the object EquilibriumDistribution, one for each isotopologue.

    :param number_of_distributions:
    :return:
    """
    if number_of_distributions == 1:
        return EquilibriumDistribution()
    return [EquilibriumDistribution()] * number_of_distributions


def get_CO2_treanor_distributions():
    """
    Get a list of Treanor distribution objects for the first three isotopologues of CO2.
    :return:
    """
    distributions = []

    # Set a Treanor distribution for the isotopologues
    for index in range(3):
        distribution = TreanorTriatomicLinearFermiResonantFractionalDistribution(
            rotational_statistical_weight=CO2_statistical_weights[index],
            molecular_constants=CO2_molecular_properties[index],
            vibrational_constants=CO2_vibrational_constants[index],
            rotational_constants=CO2_rotational_constants[index],
            fermi_resonance_constants=CO2_resonance_constants[index],
            v1_max=4,
            v2_max=8,
            v3_max=6,
            J_max=60
        )
        distributions.append(distribution)

    return distributions


def get_CO_treanor_distributions():
    """Get a list of Treanor distribution objects for the first three isotopologues of CO.

    :return:
    """
    distributions = []

    # Set a Treanor distribution for the isotopologues
    for index in range(3):
        distribution = TreanorDiatomicFractionalDistribution(
            rotational_statistical_weight=CO_statistical_weights[index],
            dunham_matrix=CO_dunham_matrices[index],
            v_max=8,
            J_max=60
        )
        distributions.append(distribution)

    return distributions


def get_NH3_distributions():
    """
    Get a list of Boltzmann distributions for 14NH3 and 15NH3 - the isotopologues included in the HITRAN database.
    :return:
    """
    # TODO: implement 15NH3, look at statistical weights + vibrational energy levels
    return [NonEquilibriumBoltzmannPyramidalTetratomicFractionDistribution(
        statistical_weight=NH3_statistical_weights[index],
        rotational_constants=NH3_rotational_constants[0],
        vibrational_constants=NH3_vibrational_constants[0]
    ) for index in range(2)]
