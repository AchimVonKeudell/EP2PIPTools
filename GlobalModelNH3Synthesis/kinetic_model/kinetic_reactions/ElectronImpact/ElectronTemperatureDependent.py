from .. import basic_reaction


class ElectronImpactReaction(basic_reaction.Reaction):
    """
        k = A * electron_density [1/s]
    """

    @property
    def rate_equation(self) -> float:
        if self.rate_constant.unit == 'cm3/s':
            return self.rate_constant.value * self.conditions['electron_density_cm3']
        raise NotImplementedError(f'Unit of {self.rate_constant=} not implemented')


class ElectronImpactTemperatureReaction(ElectronImpactReaction):
    """
        k = C0 * (T_electron / C1) ^ C2 * (T_gas / C3) ^ C4 * electron_density [1/s]

        C0: rate constants, in [1/torr/s], [cm3/s], [cm6/s]
        C1: reference electron temperature, in [eV]
        C2: electron temperature exponent, in [-]
        C3: reference gas temperature, in [K]; default = 300 K
        C4: gas temperature exponent, in [-]; default = 0.
    """

    reference_electron_temperature: basic_reaction.Constant
    exponent_electron_temperature: basic_reaction.Constant

    reference_gas_temperature: basic_reaction.Constant
    exponent_gas_temperature: basic_reaction.Constant

    def set_constants(self, *args):
        if len(args) < 2:
            raise ValueError(f'Either not enough or too many constants in {args=}!\n')
        super().set_constants(args[0])
        self.reference_electron_temperature = basic_reaction.Constant(*args[1])
        self.exponent_electron_temperature = basic_reaction.Constant(*args[2])

        if len(args) == 4:
            self.reference_gas_temperature = basic_reaction.Constant(*args[3])
            self.exponent_gas_temperature = basic_reaction.Constant(*args[4])
        else:
            self.reference_gas_temperature = basic_reaction.Constant(300, 'K')
            self.exponent_gas_temperature = basic_reaction.Constant(0., '-')

    @property
    def rate_equation(self) -> float:
        rate = super().rate_equation
        rate *= (self.conditions['electron_temperature_ev'] / self.reference_electron_temperature.value
                 ) ** self.exponent_electron_temperature.value
        if self.exponent_gas_temperature.value > 0:
            rate *= (self.conditions['gas_temperature_kelvin'] / self.reference_gas_temperature.value
                     ) ** self.exponent_gas_temperature.value
        return rate
