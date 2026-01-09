from . import kinetic_reactions
from .parse import parse_reaction


def configure_reaction(reaction_raw_data: list, species: list, conditions: dict) -> kinetic_reactions.Reaction:
    """This function conforms the raw string format of the reaction to the correct Reaction object.

    Example:
        1f,N2 + * -> N2*,ArrheniusReaction(3.85E+03 1/torr/s;0.0 kJ/mol)

        * `ArrheniusReaction` refers to the corresponding reaction formula, in this case: k = A * exp(-E/R/T)

    :param reaction_raw_data:
    :param species:
    :param conditions:
    :return:
    """
    reactants, products = parse_reaction(reaction_raw_data[1])

    reaction_full_function = reaction_raw_data[2].split('(')
    reaction_function_name = reaction_full_function[0]

    if not hasattr(kinetic_reactions, reaction_function_name):
        raise AttributeError(f'{reaction_function_name=} not known/implemented, check for spelling(?)')
    reaction_object: kinetic_reactions.Reaction = getattr(kinetic_reactions, reaction_function_name)(
        raw_reaction=reaction_raw_data[1],
        species=species,
        reactants=reactants,
        products=products,
        conditions=conditions
    )
    reaction_constants = reaction_full_function[1][:reaction_full_function[1].find(')')]
    reaction_object.set_constants(*[reaction_constant.split(' ') for reaction_constant in reaction_constants.split(';')])
    return reaction_object
