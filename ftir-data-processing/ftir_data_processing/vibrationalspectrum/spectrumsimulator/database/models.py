from dataclasses import dataclass
from . import diatomic_molecule_ids, triatomic_molecule_ids, triatomic_linear_molecule_ids, \
    triatomic_linear_fermi_resonant_molecule_ids, pyramidal_tetratomic_molecule_ids, pentatomic_molecule_ids


class Molecule:
    """
    This object represents a molecule as stored in the HITRAN database. A molecule is defined by a unique id and the
    name of the molecule
    """
    molecule_id = 0
    molecule_name = ""

    def __init__(self, molecule_id, molecule_name):
        self.molecule_id = molecule_id
        self.molecule_name = molecule_name


@dataclass
class Isotopologue:
    """
    This object represents an isotopologue of a molecule as stored in the HITRAN database. An isotopologue is defined by
    the id of the molecule to which it belongs, an order number, where 1 is the most abundant isotopologue for a given
    molecule id, the old Air Force Geophysics Laborator shorthand notation (afgl_code), the isotopologues abundance 
    fraction and the isotopologues molar mass.
    """
    molecule_id: int
    isotopologue_order_num: int
    afgl_code: str
    abundance: float
    molar_mass: float


class Structure:
    """
    This object represents an absorption structure, consisting of multiple absorption lines. The structure is defined by
    a global upper and lower quantum, and is linked to a specific isotopologue. Also, the maximum line strength is
    stored per structure to estimate relative strengths of structures.
    """
    structure_id = 0
    molecule_id = 0
    isotopologue_order_num = 0
    global_quanta_upper = None
    global_quanta_lower = None
    max_line_strength_296K = 0

    def __init__(self, molecule_id, isotopologue_order_num, global_quanta_upper, global_quanta_lower, max_line_strength_296K=0, structure_id=-1):
        self.structure_id = structure_id
        self.molecule_id = molecule_id
        self.isotopologue_order_num = isotopologue_order_num
        self.global_quanta_upper = global_quanta_upper
        self.global_quanta_lower = global_quanta_lower
        self.max_line_strength_296K = max_line_strength_296K

    @staticmethod
    def get_global_quanta_model_type_by_molecule_id(molecule_id):
        """
        Return the model type that should be used for a given molecule id. For some molecules, it is possible to extract
        more information on the global upper and lower quantum, but to store this, a different global quanta model has
        to be used.
        :param molecule_id:
        :return:
        """
        if molecule_id in diatomic_molecule_ids:
            return DiatomicGlobalQuanta
        elif molecule_id in triatomic_linear_molecule_ids:
            return TriatomicLinearGlobalQuanta
        elif molecule_id in triatomic_linear_fermi_resonant_molecule_ids:
            return TriatomicLinearFermiResonantGlobalQuanta
        elif molecule_id in pyramidal_tetratomic_molecule_ids:
            return PyramidalTetratomicGlobalQuanta
        elif molecule_id in pentatomic_molecule_ids:
            return PentatomicGlobalQuanta
        return GlobalQuanta


class GlobalQuanta:
    """
    This object represents the most general form of a global quanta in the HITRAN database structure. If no additional
    information is known on the type of molecule to which the global quanta belong, then the only known property is the
    raw string representation of the quanta in the HITRAN format.
    """
    hitran_raw_quanta = ""

    def __init__(self, hitran_raw_quanta):
        self.hitran_raw_quanta = hitran_raw_quanta

    def __repr__(self):
        return f"{self.__class__.__name__}({self.hitran_raw_quanta})"


class DiatomicGlobalQuanta(GlobalQuanta):
    """
    This object represents the global quanta for a diatomic molecule, which only has a single vibrational quantum.
    """
    vibrational_quantum = 0

    def __init__(self, vibrational_quantum):
        self.vibrational_quantum = vibrational_quantum

    @property
    def v(self):
        return self.vibrational_quantum

    @property
    def hitran_raw_quanta(self):
        return "%13s%2i" % ("", self.vibrational_quantum)


class TriatomicGlobalQuanta(GlobalQuanta):
    """
    This object represents the global quanta for a (non-linear) triatomic molecule, which has three vibrational
    modes.
    """
    vibrational_quantum_1 = 0
    vibrational_quantum_2 = 0
    vibrational_quantum_3 = 0

    def __init__(self, vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3):
        self.vibrational_quantum_1 = vibrational_quantum_1
        self.vibrational_quantum_2 = vibrational_quantum_2
        self.vibrational_quantum_3 = vibrational_quantum_3

    @property
    def v(self):
        return self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3

    @property
    def v1(self):
        return self.vibrational_quantum_1

    @property
    def v2(self):
        return self.vibrational_quantum_2

    @property
    def v3(self):
        return self.vibrational_quantum_3

    @property
    def hitran_raw_quanta(self):
        return "%7s%2i%2i%2i" % ("", self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3)


class TriatomicLinearGlobalQuanta(TriatomicGlobalQuanta):
    """
    This object represents the global quanta for a linear triatomic molecule, which has three vibrational
    modes, one of which being degenerate. Therefore, the vibrational angular momentum also has to be stored.
    """
    vibrational_angular_momentum = 0

    def __init__(self, vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3,
                 vibrational_angular_momentum):
        super().__init__(vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3)
        self.vibrational_angular_momentum = vibrational_angular_momentum

    @property
    def l(self):
        return self.vibrational_angular_momentum

    @property
    def hitran_raw_quanta(self):
        return "%7s%2i%2i%2i%2i" % (
            "",
            self.vibrational_quantum_1,
            self.vibrational_quantum_2,
            self.vibrational_angular_momentum,
            self.vibrational_quantum_3
        )


class TriatomicLinearFermiResonantGlobalQuanta(TriatomicLinearGlobalQuanta):
    """
    This object represents the global quanta for a linear triatomic molecule with a strong fermi resonance, which has
    three vibrational modes, one of which being degenerate and resonance between at least two vibrational modes.
    Therefore, the vibrational ranking index is stored additionally, which is unity for the highest level among a Fermi
    resonant pair of levels.
    """
    vibrational_ranking_index = 0

    def __init__(self, vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3,
                 vibrational_angular_momentum, vibrational_ranking_index):
        super().__init__(vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3, vibrational_angular_momentum)
        self.vibrational_ranking_index = vibrational_ranking_index

    @property
    def hitran_raw_quanta(self):
        return "%6s%2i%2i%2i%2i%1i" % (
            "",
            self.vibrational_quantum_1,
            self.vibrational_quantum_2,
            self.vibrational_angular_momentum,
            self.vibrational_quantum_3,
            self.vibrational_ranking_index
        )


class PyramidalTetratomicGlobalQuanta(TriatomicGlobalQuanta):
    """
    This object describes the global quanta of pyramidal tetratomic molecules, such as NH3 and PH3.
    """

    vibrational_quantum_4 = 0
    vibrational_angular_momentum_3 = 0
    vibrational_angular_momentum_4 = 0
    vibrational_angular_momentum = 0
    vibrational_symmetry = ''

    def __init__(self, vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3, vibrational_quantum_4,
                 vibrational_angular_momentum_3, vibrational_angular_momentum_4, vibrational_angular_momentum,
                 vibrational_symmetry):
        super().__init__(vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3)
        self.vibrational_quantum_4 = vibrational_quantum_4
        self.vibrational_angular_momentum_3 = vibrational_angular_momentum_3
        self.vibrational_angular_momentum_4 = vibrational_angular_momentum_4
        self.vibrational_angular_momentum = vibrational_angular_momentum
        self.vibrational_symmetry = vibrational_symmetry

    def __eq__(self, other):   # TODO move this to a parent object?
        if not isinstance(other, PyramidalTetratomicGlobalQuanta):
            return False
        for key, value in self.__dict__.items():
            if (key != 'vibrational_angular_momentum') & (value != other.__dict__[key]):
                return False
        return True

    @property
    def v(self):
        return self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3, \
               self.vibrational_quantum_4

    @property
    def l3(self):
        return self.vibrational_angular_momentum_3

    @property
    def l4(self):
        return self.vibrational_angular_momentum_4

    @property
    def v4(self):
        return self.vibrational_quantum_4

    @property
    def hitran_raw_quanta(self):
        if self.vibrational_symmetry in ['s', 'a']:
            return "%2i%2i%2i%2i%2s" % (
                self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3,
                self.vibrational_quantum_4, self.vibrational_symmetry)
        if self.vibrational_angular_momentum is None:
            return "%s %s%4s" % (
                f'{self.vibrational_quantum_1:1d}{self.vibrational_quantum_2:1d}'
                f'{self.vibrational_quantum_3:1d}{self.vibrational_quantum_4:1d}',
                f'{self.vibrational_angular_momentum_3:1d}' \
                f'{self.vibrational_angular_momentum_4:1d}  ',
                self.vibrational_symmetry
            )
        return "%s %s%4s" % (
            f'{self.vibrational_quantum_1:1d}{self.vibrational_quantum_2:1d}'
            f'{self.vibrational_quantum_3:1d}{self.vibrational_quantum_4:1d}',
            f'{self.vibrational_angular_momentum_3:1d}' \
            f'{self.vibrational_angular_momentum_4:1d} {self.vibrational_angular_momentum:1d}',
            self.vibrational_symmetry
        )


class PentatomicGlobalQuanta(TriatomicGlobalQuanta):
    """
    This object describes the global quanta for pentatomic and greater polyatomic molecules,
        Rothman et al Journal of Quantitative Spectroscopy & Radiative Transfer 96 (2005) 139–204
    """
    vibrational_quantum_4 = 0
    multiplicity_index = 0
    vibrational_symmetry = ''

    def __init__(self, vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3, vibrational_quantum_4,
                 multiplicity_index, vibrational_symmetry):
        super().__init__(vibrational_quantum_1, vibrational_quantum_2, vibrational_quantum_3)
        self.vibrational_quantum_4 = vibrational_quantum_4
        self.multiplicity_index = multiplicity_index
        self.vibrational_symmetry = vibrational_symmetry

    @property
    def v(self):
        return self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3, \
               self.vibrational_quantum_4

    @property
    def v4(self):
        return self.vibrational_quantum_4

    @property
    def n(self):
        return self.multiplicity_index

    @property
    def hitran_raw_quanta(self):
        return "%2i%2i%2i%2i%2i%-2s" % (
            self.vibrational_quantum_1, self.vibrational_quantum_2, self.vibrational_quantum_3,
            self.vibrational_quantum_4, self.n, self.vibrational_symmetry)


class Line:
    """
    This object represents a single absorption line, and holds multiple parameters of this line, allowing calculation of
    line profiles.
    """
    line_id = 0
    structure_id = 0
    wavenumber_vacuum = 0
    line_strength_296K = 0
    einstein_A = 0
    broadening_hw_air = 0
    broadening_hw_self = 0
    energy_lower_state = 0
    broadening_temp_coefficient = 0
    pressure_shift = 0
    local_quanta_upper = ""
    local_quanta_lower = ""
    rotational_statistical_weight_upper = 0
    rotational_statistical_weight_lower = 0

    def __init__(self, structure_id, wavenumber_vacuum, line_strength_296K, einstein_A, broadening_hw_air,
                 broadening_hw_self, energy_lower_state, broadening_temp_coefficient, pressure_shift,
                 local_quanta_upper, local_quanta_lower, rotational_statistical_weight_upper,
                 rotational_statistical_weight_lower, line_id = -1):
        self.line_id = line_id
        self.structure_id = structure_id
        self.wavenumber_vacuum = wavenumber_vacuum
        self.line_strength_296K = line_strength_296K
        self.einstein_A = einstein_A
        self.broadening_hw_air = broadening_hw_air
        self.broadening_hw_self = broadening_hw_self
        self.energy_lower_state = energy_lower_state
        self.broadening_temp_coefficient = broadening_temp_coefficient
        self.pressure_shift = pressure_shift
        self.local_quanta_upper = local_quanta_upper
        self.local_quanta_lower = local_quanta_lower
        self.rotational_statistical_weight_upper = rotational_statistical_weight_upper
        self.rotational_statistical_weight_lower = rotational_statistical_weight_lower

    @staticmethod
    def get_local_quanta_model_type_by_molecule_id(molecule_id):
        """
        Return the local quanta model type that should be used for a given molecule id. For some molecules, it is
        possible to extract more information on the local upper and lower quantum, but to store this, a different
        local quanta model has to be used.
        :param molecule_id:
        :return:
        """
        if molecule_id in [2, 5]:
            return DiatomicOrLinearLocalQuanta
        elif molecule_id == 3:
            return TriatomicNonLinearLocalQuanta
        elif molecule_id == 11:
            return PyramidalTetratomicLocalQuanta
        return LocalQuanta


class LocalQuanta:
    """
    This object represents the most general form of a local quanta in the HITRAN database structure. If no additional
    information is known on the type of molecule to which the local quanta belong, then the only known property is the
    raw string representation of the quanta in the HITRAN format.
    """
    hitran_raw_quanta = ""

    def __init__(self, hitran_raw_quanta):
        self.hitran_raw_quanta = hitran_raw_quanta

    def __repr__(self):
        return f"{self.__class__.__name__}({self.hitran_raw_quanta})"


class DiatomicOrLinearLocalQuanta(LocalQuanta):
    """
    This object represents the global quanta for a diatomic molecule, which only has a single vibrational quantum.
    """
    J = 0
    F = 0
    symmetry = ""

    def __init__(self, J, symmetry, F):
        self.J = J
        self.symmetry = symmetry
        self.F = F

    @property
    def hitran_raw_quanta(self):
        """
        Tries to construct the HITRAN database raw local quanta property. Since the HITRAN notation also includes the
        branch property, this function does not correspond exactly with HITRAN, since the lower local quanta does not
        know about the upper local quanta and thus the branch is unknown. This property is this always filled with
        spaces
        :return:
        """
        return "%6s%3i%1s%5s" % (
            "",
            self.J,
            self.symmetry,
            self.F
        )


class TriatomicNonLinearLocalQuanta(LocalQuanta):
    J = 0
    K = 0

    def __init__(self, J, K):
        self.J = J
        self.K = K

    @property
    def hitran_raw_quanta(self):
        return "%3i%3i" % (self.J, self.K)


class PyramidalTetratomicLocalQuanta(LocalQuanta):

    J = 0
    K = 0
    inversion_symmetry = ''
    symmetry_rotational = ''
    symmetry_total = ''

    def __init__(self, J, K, inversion_symmetry, symmetry_rotational, symmetry_total):
        self.J = J
        self.K = K
        self.inversion_symmetry = inversion_symmetry
        self.symmetry_rotational = symmetry_rotational
        self.symmetry_total = symmetry_total

    def __eq__(self, other):
        if isinstance(other, PyramidalTetratomicLocalQuanta):
            if self.__dict__ == other.__dict__:
                return True
        return False

    @property
    def hitran_raw_quanta(self):
        return "%2i%3i%2s %3s%3s " % (self.J, self.K, self.inversion_symmetry,
                                      self.symmetry_rotational, self.symmetry_total)

