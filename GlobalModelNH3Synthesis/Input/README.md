The chemistry is defined in the .csv files. The lines must have the following syntax:

'reaction number','chemical formula','rate coefficient function'

* 'reaction number': is only for aesthetics, not used by the python scripts
* 'chemical formula' of the reaction: species involved in the chemistry are defined here
    * use same syntax for the same molecule/atom
    * delimiter = space
    * '->' indicates direction of reaction; only use this syntax. If you want to consider the reverse reaction, then
    write this down separately (i.e. treat it as another reaction)
    * 'Ms' is an adsorbed specie
    * 's' is a free surface site
* 'rate coefficient function': refers to the reactions defined under ./model/kinetic_reactions/. The variables are
defined between the brackets
    * Reaction(rate_constant): constant reaction rate, defined in m3/s or cm3/s
    * ArrheniusReaction(C0, C1):
        * k = C0 * exp(-C1 / T_gas)
    * ArrheniusSurfaceReaction(C0, C1):
        * k = C0 * exp(-C1 / T_surface)
    * ExtendedArrheniusReaction(C0, C1, C2, C3):
        * k = C0 * (T / C1) ** C2 * exp(-C3/T)
    * ElectronImpactMaxwellianEEDF(cross_section_set): 'cross_section_set' is a filename where set is saved
    * ElectronImpactTemperatureReaction(C0, C1, C2, C3, C4):
        * k = C0 * (T_electron / C1) ^ C2 * (T_gas / C3) ^ C4 * electron_density [1/s]
        * C3, and C4 are optional
        * set C2=0 if you want to omit the linear, quadratic, .. term
    * ElectronImpactReaction: not finished yet
    * EleyRidealReaction(sticking coefficient, background gas):
        * 'sticking coefficient' of the reaction, i.e. sticking probability
        * 'background gas' is used to calculate the diffusion constant
    * LangmuirHinshelwoodReaction(nu, E_a, E_d)
        * k = nu / 4 * exp[-(E_d + E_a) / (kB * T_surface)]
        * E_d and E_a are diffusion and activation energy barriers