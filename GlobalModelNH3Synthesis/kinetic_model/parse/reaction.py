
def parse_species(species_str: str, standard_coefficient=1) -> dict:
    species_dict = {}

    coefficient = standard_coefficient
    for string_part in species_str.split(' + '):
        for specie in string_part.split():
            if specie.isdigit():
                coefficient = int(specie)
            else:
                species_dict[specie.strip()] = coefficient
                coefficient = standard_coefficient
    return species_dict


def parse_reaction(reaction_str: str) -> tuple:
    if '->' not in reaction_str:
        raise ValueError(f"reaction={reaction_str} does not contain '->'...")
    reactants_str, products_str = reaction_str.split('->')

    reactants = parse_species(reactants_str)
    products = parse_species(products_str)
    return reactants, products


if __name__ == "__main__":
    reaction = '2 H + N -> NH + H'
    print(parse_reaction(reaction))
