from kinetic_model.parse.reaction import parse_reaction


def acquire_species_from_reaction_list(reactions: list) -> list:
    reactant_species, product_species = [], []
    for reaction in reactions:
        reactants_reaction, products_reaction = parse_reaction(reaction[1])
        reactant_species += list(reactants_reaction)
        product_species += list(products_reaction)

    # re-order the list
    all_species = sorted({*reactant_species, *product_species})
    # reactant_species = sorted({*reactants_species})
    # product_species = sorted({*products_species})

    # # make sure that all products are included as well
    # for product_specie in product_species:
    #     raise ValueError(f'{reactants_species=} not equal to \n'
    #                      f'{products_species=}')
    return all_species
