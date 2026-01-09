import functools

import numpy as np

from .database import models, exomol
from . import constants
from . import weights
from . import cache


class Distribution:
    def set_temperature_K(self, *args, **kwargs):
        raise NotImplementedError


class EquilibriumDistribution(Distribution):
    def set_temperature_K(self, temperature_K):
        self.temperature_rot_K = temperature_K

    def get_partition_sum(self):
        raise NotImplementedError


class NonEquilibriumDistribution(Distribution):

    def get_rotational_partition_sum(self):
        raise NotImplementedError

    def get_vibrational_partition_sum(self):
        raise NotImplementedError

    def get_partition_sum(self):
        return self.get_vibrational_partition_sum() * self.get_rotational_partition_sum()

    def get_fractional_population(self, global_quanta: models.GlobalQuanta, local_quanta: models.LocalQuanta):
        raise NotImplementedError


class NonCharacterisableEnergyLevels:

    def get_fractional_population(self, global_quanta: models.PyramidalTetratomicGlobalQuanta,
                                  local_quanta: models.PyramidalTetratomicLocalQuanta,
                                  total_energy: float) -> float:
        raise NotImplementedError


class BoltzmannDiatomicFractionalDistribution(EquilibriumDistribution):
    """
    Object to represent the fractional distribution of a molecule among its states. This object does this for
    the specific case of a diatomic, based on Dunham coefficients, as described in Mantz, A. W., Maillard, J.-P., Roh,
    W. B., & Narahari Rao, K. (1975). Ground state molecular constants of 12C16O. Journal of Molecular Spectroscopy,
    57(1), 155–159. https://doi.org/10.1016/0022-2852(75)90049-1
    """
    def __init__(self, rotational_statistical_weight: weights.RotationalStatisticalWeight, dunham_matrix, J_max = 100, v_max = 10):
        """
        Initialize a distribution indicating the fraction of diatomic molecules in a given state
        :param rotational_statistical_weight: An object indicating how the state dependent and state independent
        rotational statistical weights are for the relevant isotopologue.
        :param dunham_matrix: The dunham matrix for a diatomic
        :param J_max: The maximum rotational level to include for calculation of the rotational partition sum
        :param v_max: The maximum vibrational level to include for calculation of the vibrational partition sum
        """
        self._rotational_statistical_weight = rotational_statistical_weight
        self._dunham_matrix = np.matrix(dunham_matrix)
        self._l_max, self._m_max = self._dunham_matrix.shape

        self._J_max = J_max
        self._v_max = v_max

        self.temperature_rot_K = 0

    @cache.reset(check_parameter_change=True)
    def set_temperature_K(self, temperature_rot_K):
        """
        Set the rotational temperature to use for calculation of densities.
        :param temperature_rot_K:
        :return:
        """
        self.temperature_rot_K = temperature_rot_K

    def get_rotational_energy(self, global_quanta: models.DiatomicGlobalQuanta,
                              local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Calculate the rotational energy for a level determined by the global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        energy = 0

        for l in range(self._l_max):
            for m in range(1, self._m_max):
                energy += self._dunham_matrix.item(l, m) * (global_quanta.vibrational_quantum + 0.5)**l * \
                    local_quanta.J**m * (local_quanta.J + 1)**m

        return energy

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def get_vibrational_energy(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Calculate the vibrational energy for a level determined by the global quanta. The ground-state energy
        is subtracted, since most processes are defined relative to the ground-state energy.
        :param global_quanta:
        :return:
        """
        energy = 0
        energy_ground_state = 0

        for l in range(self._l_max):
            energy += self._dunham_matrix.item(l, 0) * (global_quanta.vibrational_quantum + 0.5)**l
            energy_ground_state += self._dunham_matrix.item(l, 0) * 0.5**l

        return energy - energy_ground_state

    def get_total_energy(self, global_quanta: models.DiatomicGlobalQuanta,
                         local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the total energy (rotational+vibrational) for a level determined by the global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self.get_vibrational_energy(global_quanta) + self.get_rotational_energy(global_quanta, local_quanta)

    def _get_unnormalized_rotational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta,
                                                           local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the unnormalized rotational fraction population. This fractional population is not normalized because it has
        not been divided by the rotational partition sum yet.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self._rotational_statistical_weight.get_rotational_weight(global_quanta, local_quanta) * \
            np.exp(-constants.c2 * self.get_rotational_energy(global_quanta, local_quanta) / self.temperature_rot_K)

    @cache.cache_by_global_quanta_decorator
    def get_rotational_partition_sum(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Retrieve the rotational partition sum by summing over all unnormalized rotational fraction populations in a
        given global quanta.
        :param global_quanta:
        :return:
        """
        return np.sum(
            [self._get_unnormalized_rotational_fractional_population(
                global_quanta,
                models.DiatomicOrLinearLocalQuanta(
                    J=J,
                    symmetry="",
                    F=""
                )
            ) for J in range(self._J_max)])

    def get_rotational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta,
                                             local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the (normalized) rotational fractional population for a given state described by given global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self._get_unnormalized_rotational_fractional_population(global_quanta, local_quanta) / \
            self.get_rotational_partition_sum(global_quanta)

    @cache.cache_by_global_quanta_decorator
    def _get_unnormalized_vibrational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Get the unnormalized vibrational fraction population. This fractional population is not normalized because it has
        not been divided by the vibrational partition sum yet.
        :param global_quanta:
        :return:
        """
        return np.exp(-constants.c2 * self.get_vibrational_energy(global_quanta) / self.temperature_rot_K)

    def get_vibrational_partition_sum(self):
        """
        Retrieve the vibrational partition sum by summing over all unnormalized vibrational fraction populations.
        :return:
        """
        return np.sum(
            [self._get_unnormalized_vibrational_fractional_population(
                models.DiatomicGlobalQuanta(vibrational_quantum=v)
            ) for v in range(self._v_max)])

    @cache.cache_by_global_quanta_decorator
    def get_vibrational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Get the (normalized) vibrational fractional population for a given state described by given global quanta.
        :param global_quanta:
        :return:
        """
        return self._get_unnormalized_vibrational_fractional_population(global_quanta) / \
            self.get_vibrational_partition_sum()

    def get_fractional_population(self, global_quanta: models.GlobalQuanta, local_quanta: models.LocalQuanta):
        """
        Calculate the fractional population in a given level determined by the local and global quanta. By calculating
        the fractional population, it is possible to easily find the 'real' number density in the given level by
        multiplying the fractional population by the total number density for the molecule under consideration.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self.get_vibrational_fractional_population(global_quanta) * \
            self.get_rotational_fractional_population(global_quanta, local_quanta)


class TreanorDiatomicFractionalDistribution(BoltzmannDiatomicFractionalDistribution, NonEquilibriumDistribution):
    """
    Object representing a Treanor distribution for a diatomic molecule where a seperate temperature describes the
    vibrational distribution among the vibrational states.
    """
    def __init__(self, rotational_statistical_weight, dunham_matrix, J_max=100, v_max=10):
        """
        Initialize a Treanor distribution indicating the fraction of diatomic molecules in a given state
        :param rotational_statistical_weight: An object indicating how the state dependent and state independent
        rotational statistical weights are for the relevant isotopologue.
        :param dunham_matrix: The dunham matrix for a diatomic
        :param J_max: The maximum rotational level to include for calculation of the rotational partition sum
        :param v_max: The maximum vibrational level to include for calculation of the vibrational partition sum
        """
        super().__init__(rotational_statistical_weight, dunham_matrix, J_max, v_max)

        self._v_dunham_max, _ = self._dunham_matrix.shape
        self._G01 = self._calculate_G01()

    def _calculate_G01(self):
        """
        Calculate the spacing between v=0,J=0 and v=1,J=0 as a zeroth order approximation of the vibrational spacing.
        :return:
        """
        return \
            self.get_vibrational_energy(models.DiatomicGlobalQuanta(vibrational_quantum=1)) - \
            self.get_vibrational_energy(models.DiatomicGlobalQuanta(vibrational_quantum=0))

    @cache.reset(check_parameter_change=True)
    def set_temperature_K(self, temperature_rot_K, temperature_vib_K):
        """
        Set the vibrational and rotational temperatures to use for calculation of densities.
        :param temperature_rot_K:
        :param temperature_vib_K:
        :return:
        """
        self.temperature_rot_K = temperature_rot_K
        self.temperature_vib_K = temperature_vib_K

    def _get_unnormalized_vibrational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Get the unnormalized vibrational fraction population. This fractional population is not normalized because it has
        not been divided by the vibrational partition sum yet.
        :param global_quanta:
        :return:
        """
        return np.exp(-constants.c2 * (
                global_quanta.vibrational_quantum * self._G01 / self.temperature_vib_K +
                global_quanta.vibrational_quantum * (global_quanta.vibrational_quantum - 1) * self._dunham_matrix.item(2, 0) / self.temperature_rot_K
        ))

    def get_vibrational_partition_sum(self):
        """
        Retrieve the vibrational partition sum by summing over all unnormalized vibrational fraction populations.
        :return:
        """
        return np.sum(
            [self._get_unnormalized_vibrational_fractional_population(
                models.DiatomicGlobalQuanta(vibrational_quantum=v)
            ) for v in range(self._v_max)])

    @cache.cache_by_global_quanta_decorator
    def get_vibrational_fractional_population(self, global_quanta: models.DiatomicGlobalQuanta):
        """
        Get the (normalized) vibrational fractional population for a given state described by given global quanta.
        :param global_quanta:
        :return:
        """
        return self._get_unnormalized_vibrational_fractional_population(global_quanta) / \
            self.get_vibrational_partition_sum()

    def get_fractional_population(self, global_quanta: models.GlobalQuanta, local_quanta: models.LocalQuanta):
        """
        Calculate the fractional population in a given level determined by the local and global quanta. By calculating
        the fractional population, it is possible to easily find the 'real' number density in the given level by
        multiplying the fractional population by the total number density for the molecule under consideration.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self.get_vibrational_fractional_population(global_quanta) * \
            self.get_rotational_fractional_population(global_quanta, local_quanta)


class BoltzmannTriatomicLinearFermiResonantFractionalDistribution(EquilibriumDistribution):
    """
    Object to represent a Boltzmann distribution for a triatomic linear molecule with a strong Fermi resonance. 
    
    The vibrational modes are split according to the following convention, which follows HITRAN
    * All symmetric stretch vibrations and its Fermi resonant states are grouped under the symmetric stretch mode. This
      includes vibrations in the bending mode with v2 != l2
    * All bending mode vibrations which do not have a Fermi resonant with any other states (e.g. v2=l2)
    * All asymmetric stretch modes are grouped under the asymmetric stretch vibration
    """
    def __init__(self, rotational_statistical_weight: weights.RotationalStatisticalWeight, molecular_constants,
                 vibrational_constants, fermi_resonance_constants, rotational_constants, J_max=100, v1_max=10,
                 v2_max=10, v3_max=10):
        """
        Calculation of a Boltzmann distribution for a triatomic molecule with strong Fermi resonance. These calculations
        are based on the following paper:
        Stephenson, D. A., & Blint, R. J. (1979). Theoretical Fitting of Computer Processed Laser Raman Spectra from 
        Methane- and Propane-Air Flames. Applied Spectroscopy, 33(1), 41–45. https://doi.org/10.1366/0003702794926218
        :param rotational_statistical_weight: An object indicating how the state dependent and state independent
        rotational statistical weights are for the relevant isotopologue.
        :param molecular_constants: An array containing some molecular properties such as the classical symmetry
        number, kappa and the nuclear spin degeneracy.
        :param vibrational_constants: A multi dimensional array containing the constants to calculate the unperturbed 
        vibrational levels according to Stephenson et al. Index 0 should contain the first order coefficients (3 total),
        labeled as w_oi in the paper. Index 1 contains the second order coefficients (an upper triangular 3x3 matrix),
        index 2 contains the third order coefficients (an upper triangular 3x3x3 matrix), and index 4 contains the
        angular momentum coupling constant xl2l2.
        :param fermi_resonance_constants: An array containing the constants for resonance between Fermi levels. The
        value with index 0 is defined by Stephenson et al. as W0, where index 1, 2 and 3 are the constants associated
        with vibrational mode v1, v2 and v3 respectively. Index 4 is the constant associated to the vibrational angular
        momentum.
        :param rotational_constants: The rotational constants involved to calculate the energies of rotational levels.
        This is an array that includes multiple constants, including the rotational constant B at index 0, first order
        vibrational correction factors (3x1 array) at index 1, and second order vibrational correction factors 
        (3x3 matrix) at index 2. Rotational D and H constants are stored at index 3 and 4 respectively.
        :param J_max: 
        :param v1_max:
        :param v2_max:
        :param v3_max:
        """
        super().__init__()

        self._rotational_statistical_weight = rotational_statistical_weight
        self._molecular_constants = molecular_constants
        self._vibrational_constants = vibrational_constants
        self._fermi_resonance_constants = fermi_resonance_constants
        self._rotational_constants = rotational_constants

        self._J_max = J_max
        self._v1_max = v1_max
        self._v2_max = v2_max
        self._v3_max = v3_max

        self.temperature_rot_K = 0

    @cache.reset(check_parameter_change=True)
    def set_temperature_K(self, temperature_rot_K):
        """
        Set the rotational temperature to use for calculation of densities.
        :param temperature_rot_K:
        :return:
        """
        self.temperature_rot_K = temperature_rot_K

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_absolute_unperturbed_vibrational_energy(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the absolute energy of a state based on the vibrational quanta for unperturbed levels (e.g. no Fermi
        resonance). Since this is the absolute energy, this value returns a non-zero energy for vibrational level
        (v1=0, v2=0, v3=0, l2=0).
        :param global_quanta:
        :return:
        """
        # define the vibrational degeneracies of each of the vibrational modes
        degeneracy = [1, 2, 1]

        return np.sum([self._vibrational_constants[0][i] * (global_quanta.v[i] + degeneracy[i] / 2) for i in range(3)]) + \
               np.sum([self._vibrational_constants[1][i, j] * (global_quanta.v[i] + degeneracy[i] / 2) * (global_quanta.v[j] + degeneracy[j] / 2) for i in range(3) for j in range(3)]) + \
               self._vibrational_constants[2] * global_quanta.l ** 2

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_unperturbed_vibrational_energy(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the energy of a state, relative to the vibrational ground state, based on the vibrational quanta for
        unperturbed levels (e.g. no Fermi resonance)
        :param global_quanta: 
        :return: 
        """
        ground_state_global_quanta = models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=0,
            vibrational_quantum_2=0,
            vibrational_quantum_3=0,
            vibrational_angular_momentum=0,
            vibrational_ranking_index=1
        )

        return self._get_absolute_unperturbed_vibrational_energy(global_quanta) - \
               self._get_absolute_unperturbed_vibrational_energy(ground_state_global_quanta)

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_perturbed_vibrational_coupling(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Calculate the vibrational coupling between the level (v1, v2, v3, l2) and (v1-1, v2+2, v3, l2) to be used as 
        off-axis elements in the Hamiltonian matrix
        :param global_quanta: Defines the level (v1, v2, v3, l2)
        :return: 
        """
        We = self._fermi_resonance_constants[0]
        We -= sum([self._fermi_resonance_constants[i+1]*global_quanta.v[i] for i in range(3)])

        return 0.5 * We * ((global_quanta.v2 + 2) ** 2 - global_quanta.l ** 2) ** 0.5 * global_quanta.v1 ** 0.5

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_resonant_states(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        For a given vibrational level, described by global_quanta, find all Fermi-resonant states and return the global
        quanta of these resonant levels. NOTE: The level for which resonant states are requested is always included in 
        the result.
        :param global_quanta: 
        :return: 
        """
        resonant_states = []
        resonant_state_count = int(global_quanta.v1 + (global_quanta.v2 - global_quanta.l) / 2 + 1)

        for i in range(resonant_state_count):
            # decrease v1 stepwise by i, while increasing v2 stepwise by 2i
            v1 = resonant_state_count - (i + 1)
            v2 = global_quanta.l + 2 * ((i + 1) - 1)

            resonant_states.append(models.TriatomicLinearFermiResonantGlobalQuanta(
                vibrational_quantum_1=v1,
                vibrational_quantum_2=v2,
                vibrational_quantum_3=global_quanta.v3,
                vibrational_angular_momentum=global_quanta.vibrational_angular_momentum,
                vibrational_ranking_index=0  # as long as these quanta are used within this class, this is irrelevant
            ))

        return resonant_states

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_hamiltonian(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve the Hamiltonian for a given vibrational level described by the global quanta. The size of the 
        Hamiltonian depends on the amount of resonant states for the given vibrational state. If there are 3 resonant
        states, the hamiltonian is 3x3 in size.
        :param global_quanta: 
        :return: 
        """
        resonant_states = self._get_resonant_states(global_quanta)
        resonant_state_count = len(resonant_states)
        hamiltonian = np.zeros([resonant_state_count, resonant_state_count])

        for i, resonant_state in enumerate(resonant_states):
            for j in range(resonant_state_count):
                if i == j:
                    hamiltonian[i, j] = self._get_unperturbed_vibrational_energy(resonant_state)
                elif j == i + 1:
                    Wij = self._get_perturbed_vibrational_coupling(resonant_state)
                    hamiltonian[i, j] = Wij
                    hamiltonian[j, i] = Wij
        return hamiltonian

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_eigenstates(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrive the eigenvalues and eigen vectors of the hamiltonian for a given vibrational level defined by its
        global quanta and return these as a tuple.
        :param global_quanta: 
        :return: 
        """
        eigen_energies, eigen_vectors = np.linalg.eig(self._get_hamiltonian(global_quanta))

        return eigen_energies, eigen_vectors

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_resonance_energies(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve the energies of a Fermi resonant vibrational level, sorted from low to high energy. 
        NOTE: This also works for a non-resonant level, since an array with a single energy will then be returned
        :param global_quanta: 
        :return: 
        """
        resonance_energies, _ = self._get_eigenstates(global_quanta)
        return resonance_energies

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_resonance_weights(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        The resonance weights, indicating how much each of the resonant levels contributes to the final state is
        calculated and returned by this function by squaring the eigen-vectors element-wise.
        :param global_quanta: 
        :return: 
        """
        _, resonance_weights = self._get_eigenstates(global_quanta)
        return np.square(resonance_weights)

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_resonance_sorting_indexes_desc(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve the indexes that sort the resonance energies from high to low energy. This is used to correspond 
        energies to their HITRAN ranking index, which distinguishes a Fermi resonant group by their energy sorting.
        :param global_quanta: 
        :return: 
        """
        resonance_energies, _ = self._get_eigenstates(global_quanta)
        return np.argsort(resonance_energies)[::-1]

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def get_vibrational_energy(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Calculate the vibrational energy for a level determined by the global quanta.
        :param global_quanta:
        :return:
        """
        # get the resonance energies for a given level, and invert it, so it goes from high energy to low energy
        resonance_sorting_indexes = self._get_resonance_sorting_indexes_desc(global_quanta)
        resonance_energies = self._get_resonance_energies(global_quanta)[resonance_sorting_indexes]

        # the ranking index gives which state belongs to which resonant energy, where a ranking index of 1 corresponds
        # to the highest energy of a Fermi resonant group, with increasing numbers corresponding to lower energies
        return resonance_energies[global_quanta.vibrational_ranking_index-1]

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_unperturbed_rotational_constant(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Calculate the unperturbed (e.g. without Fermi resonance) rotational B constant for a given vibrational level.
        :param global_quanta: 
        :return: 
        """
        # define the vibrational degeneracies of each of the vibrational modes
        degeneracy = [1, 2, 1]

        return self._rotational_constants[0] - \
               np.sum([self._rotational_constants[1][i]*(global_quanta.v[i] + degeneracy[i]/2) for i in range(3)])

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_rotational_constant(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Calculate the rotational constant for a given vibrational level, while taking Fermi resonance into account. This
        effectively causes mixing of the levels by a given weight per level. The effective rotational B constant is then
        the weighted average of the unperturbed rotational B constants.
        :param global_quanta: 
        :return: 
        """
        sorting_indexes = self._get_resonance_sorting_indexes_desc(global_quanta)
        resonant_states = np.array(self._get_resonant_states(global_quanta))[sorting_indexes]
        resonance_weights = self._get_resonance_weights(global_quanta)[sorting_indexes]

        B_eff = 0

        for i, resonant_state in enumerate(resonant_states):
            B_eff += resonance_weights[global_quanta.vibrational_ranking_index-1][i] * \
                 self._get_unperturbed_rotational_constant(resonant_state)

        return B_eff

    def get_rotational_energy(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta,
                              local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Calculate the rotational energy for a level determined by the global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self._get_rotational_constant(global_quanta) * local_quanta.J * (local_quanta.J + 1) \
               - self._rotational_constants[2] * local_quanta.J**2 * (local_quanta.J + 1)**2 \
               + self._rotational_constants[3] * local_quanta.J**3 * (local_quanta.J + 1)**3

    def get_total_energy(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta,
                         local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the total energy (rotational+vibrational) for a level determined by the global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self.get_vibrational_energy(global_quanta) + self.get_rotational_energy(global_quanta, local_quanta)

    def _get_unnormalized_rotational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta,
                                                           local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the unnormalized rotational fraction population. This fractional population is not normalized because it has
        not been divided by the rotational partition sum yet.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self._rotational_statistical_weight.get_rotational_weight(global_quanta, local_quanta) * \
               np.exp(-constants.c2 * self.get_rotational_energy(global_quanta, local_quanta) / self.temperature_rot_K)

    @cache.cache_fixed_value
    def get_rotational_partition_sum(self):
        """
        Calculate the rotational partition sum based on McDowell, R. S. (1988). Rotational partition functions for linear molecules. The Journal of Chemical Physics, 88(1), 356–361. https://doi.org/10.1063/1.454608
        :param global_quanta:
        :return:
        """
        sigma = self._molecular_constants[0]
        kappa = self._molecular_constants[1]
        I = self._molecular_constants[2]

        d = self._rotational_constants[2]/self._rotational_constants[0]
        h = self._rotational_constants[3]/self._rotational_constants[0]
        beta = constants.h * constants.c * self._rotational_constants[0]*100 / (constants.k_B * self.temperature_rot_K)

        f_c = 1 + 2 * d * (3 - beta) / (3 * beta) + 6 * (2 * d ** 2 - h) / beta ** 2 + 120 * d * (d ** 2 - h) / beta ** 3
        return sigma ** -1 * I ** 2 * np.exp(beta / 3) * beta ** -1 * (1 + beta ** 2 / 90 + 8 * beta ** 3 / 2835 + kappa * I ** -1 * np.pi ** 1.5 * np.exp(-beta / 12) * np.exp(-np.pi ** 2 / (4 * beta)) * beta ** -0.5) * f_c

    def get_rotational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta,
                                             local_quanta: models.DiatomicOrLinearLocalQuanta):
        """
        Get the (normalized) rotational fractional population for a given state described by given global and local quanta.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self._get_unnormalized_rotational_fractional_population(global_quanta, local_quanta) / \
               self.get_rotational_partition_sum()

    def _get_unnormalized_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the unnormalized vibrational fraction population. This fractional population is not normalized because it has
        not been divided by the vibrational partition sum yet.
        :param global_quanta:
        :return:
        """
        # the actual degeneracy of the v2 state is v2+1, since the rotation of the molecule and its vibrational angular
        # momentum can be in the same direction or in opposite directions (labeled in HITRAN as e and f). For v2=0, this
        # l2 is also zero, and the f symmetry is not possible. However, since part of these degenerate states are also
        # included in Fermi resonant states, only v2=l2 has to be included here with a degeneracy of 2 if v2 != 0
        if global_quanta.l > 0:
            vib_degeneracy = 2
        else:
            vib_degeneracy = 1

        return vib_degeneracy * np.exp(-constants.c2 * self.get_vibrational_energy(global_quanta) / self.temperature_rot_K)

    @cache.cache_fixed_value
    def get_vibrational_partition_sum(self):
        """
        Retrieve the vibrational partition sum by summing over all unnormalized vibrational fraction populations.
        :return:
        """
        vibrational_partition_sum = 0

        for v1 in range(self._v1_max):
            for l2 in range(self._v2_max):
                for v3 in range(self._v3_max):
                    for ranking_index in range(v1+1):
                        vibrational_partition_sum += self._get_unnormalized_vibrational_fractional_population(
                            models.TriatomicLinearFermiResonantGlobalQuanta(
                                vibrational_quantum_1=v1,
                                vibrational_quantum_2=l2,
                                vibrational_quantum_3=v3,
                                vibrational_angular_momentum=l2,
                                vibrational_ranking_index=ranking_index+1
                            )
                        )

        return vibrational_partition_sum

    def get_partition_sum(self):
        """"""
        return self.get_vibrational_partition_sum() * self.get_rotational_partition_sum()

    def get_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the (normalized) vibrational fractional population for a given state described by given global quanta.
        :param global_quanta:
        :return:
        """
        return self._get_unnormalized_vibrational_fractional_population(global_quanta) / \
            self.get_vibrational_partition_sum()

    def get_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta,
                                  local_quanta: models.LocalQuanta):
        """
        Calculate the fractional population in a given level determined by the local and global quanta. By calculating
        the fractional population, it is possible to easily find the 'real' number density in the given level by
        multiplying the fractional population by the total number density for the molecule under consideration.
        :param global_quanta:
        :param local_quanta:
        :return:
        """
        return self.get_vibrational_fractional_population(global_quanta) * \
            self.get_rotational_fractional_population(global_quanta, local_quanta)


class TreanorTriatomicLinearFermiResonantFractionalDistribution(BoltzmannTriatomicLinearFermiResonantFractionalDistribution, NonEquilibriumDistribution):
    """
    Object to represent a Treanor distribution for a triatomic linear molecule with a strong Fermi resonance. Since
    this object represents a non-equilibrium distribution, the vibrational modes are assumed independent of eachother,
    thus omitting interaction between these modes.

    Hence, the vibrational energies are only described by a harmonic and anharmonic coefficient (G and omegax
    respectively). Therefore, the energy of higher levels is generally slightly overestimated, and thus this method
    breaks down for high vibrational levels.

    The vibrational modes are split according to the following convention, which follows HITRAN
    * All symmetric stretch vibrations and its Fermi resonant states are grouped under the symmetric stretch mode. This
      includes vibrations in the bending mode with v2 != l2.
    * All bending mode vibrations which do not have a Fermi resonant with any other states (e.g. v2=l2)
    * All asymmetric stretch modes are grouped under the asymmetric stretch vibration

    Furthermore, it is assumed that the relaxation rates of the symmetric stretch and bending mode are much larger
    than those of the asymmetric stretch, due to their Fermi interaction and generally lower energy spacing. Hence, the
    symmetric stretch and bending mode are not approximated by a Treanor, but a Boltzmann distribution, which is a good
    approximation for lower vibrational levels. This is due to the scaling of the harmonic part of the energy with
    vibrational temperature and the anharmonic part with the rotational temperature. For Fermi resonant levels in the
    symmetric stretch and bending mode, it is not clear whether the Fermi shifts belong to the harmonic or anharmonic
    energy contribution and thus the total energy is scaled with the vibrational temperature using a Boltzmann
    distribution.
    """
    def __init__(self, rotational_statistical_weight: weights.RotationalStatisticalWeight, molecular_constants,
                 vibrational_constants, fermi_resonance_constants, rotational_constants, J_max=100, v1_max=10, v2_max=10, v3_max=10):
        """
        Calculation of a Treanor distribution for a triatomic molecule with strong Fermi resonance. These calculations
        are based on the following paper:
        Stephenson, D. A., & Blint, R. J. (1979). Theoretical Fitting of Computer Processed Laser Raman Spectra from
        Methane- and Propane-Air Flames. Applied Spectroscopy, 33(1), 41–45. https://doi.org/10.1366/0003702794926218
        :param rotational_statistical_weight: An object indicating how the state dependent and state independent
        rotational statistical weights are for the relevant isotopologue.
        :param molecular_constants: An array containing some molecular properties such as the classical symmetry
        number, kappa and the nuclear spin degeneracy.
        :param vibrational_constants: A multi dimensional array containing the constants to calculate the unperturbed
        vibrational levels according to Stephenson et al. Index 0 should contain the first order coefficients (3 total),
        labeled as w_oi in the paper. Index 1 contains the second order coefficients (an upper triangular 3x3 matrix),
        index 2 contains the third order coefficients (an upper triangular 3x3x3 matrix), and index 4 contains the
        angular momentum coupling constant g22.
        :param fermi_resonance_constants: An array containing the constants for resonance between Fermi levels. The
        value with index 0 is defined by Stephenson et al. as W0, where index 1, 2 and 3 are the constants associated
        with vibrational mode v1, v2 and v3 respectively. Index 4 is the constant associated to the vibrational angular
        momentum.
        :param rotational_constants: The rotational constants involved to calculate the energies of rotational levels.
        This is an array that includes multiple constants, including the rotational constant B at index 0, first order
        vibrational correction factors (3x1 array) at index 1, and second order vibrational correction factors
        (3x3 matrix) at index 2. Rotational D and H constants are stored at index 3 and 4 respectively.
        :param J_max:
        :param v1_max:
        :param v2_max:
        :param v3_max:
        """
        super().__init__(rotational_statistical_weight, molecular_constants, vibrational_constants, fermi_resonance_constants, rotational_constants, J_max, v1_max, v2_max, v3_max)

        self.temperature_rot_K = 0
        self.temperature_vib_12_K = 0
        self.temperature_vib_3_K = 0

    @property
    @cache.cache_fixed_value(immutable=True)
    def _G1_v3(self):
        """
        Calculate the harmonic energy spacing between vibrational levels in the asymmetric stretch, where the total
        energy is described by:
        E = G1 * v3 - omegax * v3 * (v3-1)
        :return:
        """
        return self._vibrational_constants[0][2] + 0.5 * self._vibrational_constants[1][0, 2] + \
               self._vibrational_constants[1][1, 2] + 2 * self._vibrational_constants[1][2, 2]

    @property
    @cache.cache_fixed_value(immutable=True)
    def _omegax_v3(self):
        """
        Calculate the anharmonic energy coefficient for vibrational levels in the asymmetric stretch, where the total
        energy is described by:
        E = G1 * v3 - omegax * v3 * (v3-1)
        :return:
        """
        return -self._vibrational_constants[1][2, 2]

    @cache.reset(check_parameter_change=True)
    def set_temperature_K(self, temperature_rot_K, temperature_vib_12_K, temperature_vib_3_K):
        """
        Set the rotational and vibrational temperatures to use for calculation of densities.
        :param temperature_rot_K: The temperature used for the rotational level distribution, as well as for the
        anharmonic energy term in the Treanor distribution of the asymmetric stretch vibrational mode.
        :param temperature_vib_12_K: The vibrational temperature used for a Boltzmann distribution of the symmetric
        stretch and bending mode
        :param temperature_vib_3_K: The vibrational temperature used for the harmonic energy term in the Treanor
        distribution of the asymmetric stretch vibrational mode.
        :return:
        """
        self.temperature_rot_K = temperature_rot_K
        self.temperature_vib_12_K = temperature_vib_12_K
        self.temperature_vib_3_K = temperature_vib_3_K

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_symmetric_stretch_resonant_quanta(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve a separate global quanta object for the resonant modes of the symmetric stretch and bending mode,
        while omitting other vibrational modes.
        :param global_quanta:
        :return:
        """
        return models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=global_quanta.v1,
            vibrational_quantum_2=global_quanta.v2 - global_quanta.l,
            vibrational_quantum_3=0,
            vibrational_angular_momentum=0,
            vibrational_ranking_index=global_quanta.vibrational_ranking_index
        )

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_bending_mode_non_resonant_quanta(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve a separate global quanta object for the non-resonant levels in the bending mode, while omitting other
        vibrational modes.
        :param global_quanta:
        :return:
        """
        return models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=0,
            vibrational_quantum_2=global_quanta.l,
            vibrational_quantum_3=0,
            vibrational_angular_momentum=global_quanta.l,
            vibrational_ranking_index=1
        )

    @cache.cache_by_global_quanta_decorator(immutable=True)
    def _get_asymmetric_stretch_quanta(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Retrieve a separate global quanta object for the asymmetric stretch, while omitting other vibrational modes.
        :param global_quanta:
        :return:
        """
        return models.TriatomicLinearFermiResonantGlobalQuanta(
            vibrational_quantum_1=0,
            vibrational_quantum_2=0,
            vibrational_quantum_3=global_quanta.v3,
            vibrational_angular_momentum=0,
            vibrational_ranking_index=1
        )

    @cache.cache_by_global_quanta_decorator
    def _get_unnormalized_symmetric_stretch_resonant_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the unnormalized vibrational fractional population in the resonant states of the symmetric stretch. This 
        includes states that are in the bending mode, but resonant with the asymmetric stretch, following HITRAN 
        conventions. This fractional population is not normalized because it has not been divided by the vibrational 
        partition sum yet.
        :param global_quanta:
        :return:
        """
        symmetric_stretch_resonant_global_quanta = self._get_symmetric_stretch_resonant_quanta(global_quanta)
        return np.exp(-constants.c2 * self.get_vibrational_energy(symmetric_stretch_resonant_global_quanta) / self.temperature_vib_12_K)

    @cache.cache_fixed_value
    def get_symmetric_stretch_resonant_vibrational_partition_sum(self):
        """
        Get the vibrational partition sum of the symmetric stretch mode, including all resonant states (e.g. different
        ranking indexes).
        :return: 
        """
        vibrational_partition_sum = 0

        for v1 in range(self._v1_max):
            for ranking_index in range(v1+1):
                vibrational_partition_sum += self._get_unnormalized_symmetric_stretch_resonant_vibrational_fractional_population(
                    models.TriatomicLinearFermiResonantGlobalQuanta(
                        vibrational_quantum_1=v1,
                        vibrational_quantum_2=0,
                        vibrational_quantum_3=0,
                        vibrational_angular_momentum=0,
                        vibrational_ranking_index=ranking_index+1
                    )
                )

        return vibrational_partition_sum

    @cache.cache_by_global_quanta_decorator
    def get_symmetric_stretch_resonant_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the fractional population of a state in the symmetric stretch or one of its Fermi resonant states, described 
        by the given global quanta, such that the sum of the fractional populations of all symmetric stretch states 
        equals unity. To retrieve the number density in a given symmetric stretch state, this fractional population has 
        to be multiplied by the total number density of the isotopologue.
        :param global_quanta: 
        :return: 
        """
        return self._get_unnormalized_symmetric_stretch_resonant_vibrational_fractional_population(global_quanta) /\
               self.get_symmetric_stretch_resonant_vibrational_partition_sum()

    @cache.cache_by_global_quanta_decorator
    def _get_unnormalized_bending_mode_non_resonant_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the unnormalized vibrational fractional population in the non resonant states of the bending mode (e.g.
        v2 = l2). This fractional population is not normalized because it has not been divided by the vibrational 
        partition sum yet.
        :param global_quanta:
        :return:
        """
        bending_mode_non_resonant_global_quanta = self._get_bending_mode_non_resonant_quanta(global_quanta)

        # the actual degeneracy of the v2 state is v2+1, since the rotation of the molecule and its vibrational angular
        # momentum can be in the same direction or in opposite directions (labeled in HITRAN as e and f). For v2=0, this
        # l2 is also zero, and the f symmetry is not possible. However, since part of these degenerate states are also
        # included in Fermi resonant states, only v2=l2 has to be included here with a degeneracy of 2 if v2 != 0
        if bending_mode_non_resonant_global_quanta.v2 > 0:
            vib_degeneracy = 2
        else:
            vib_degeneracy = 1

        return vib_degeneracy * np.exp(-constants.c2 * self.get_vibrational_energy(bending_mode_non_resonant_global_quanta) / self.temperature_vib_12_K)

    @cache.cache_fixed_value
    def get_bending_mode_non_resonant_vibrational_partition_sum(self):
        """
        Get the vibrational partition sum of the bending mode, excluding all Fermi resonant states where v2 != l2
        :return: 
        """
        vibrational_partition_sum = 0

        for v2 in range(self._v2_max):
            vibrational_partition_sum += self._get_unnormalized_bending_mode_non_resonant_vibrational_fractional_population(
                models.TriatomicLinearFermiResonantGlobalQuanta(
                    vibrational_quantum_1=0,
                    vibrational_quantum_2=v2,
                    vibrational_quantum_3=0,
                    vibrational_angular_momentum=v2,
                    vibrational_ranking_index=1
                )
            )

        return vibrational_partition_sum

    @cache.cache_by_global_quanta_decorator
    def get_bending_mode_non_resonant_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the fractional population of a state in a non-resonant state in the bending mode, described by the given 
        global quanta, such that the sum of the fractional populations of all symmetric stretch states equals unity. 
        To retrieve the number density in a given bending mode state, this fractional population has to be multiplied 
        by the total number density of the isotopologue.
        :param global_quanta: 
        :return: 
        """
        return self._get_unnormalized_bending_mode_non_resonant_vibrational_fractional_population(global_quanta) / self.get_bending_mode_non_resonant_vibrational_partition_sum()

    @cache.cache_by_global_quanta_decorator
    def _get_unnormalized_asymmetric_stretch_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the unnormalized vibrational fractional population in the asymmetric stretch. This fractional population is 
        not normalized because it has not been divided by the vibrational partition sum yet.
        :param global_quanta:
        :return:
        """
        return np.exp(-constants.c2 * (self._G1_v3 * global_quanta.v3 / self.temperature_vib_3_K - self._omegax_v3 * global_quanta.v3 * (global_quanta.v3 - 1) / self.temperature_rot_K))

    @cache.cache_fixed_value
    def get_asymmetric_stretch_vibrational_partition_sum(self):
        """
        Get the vibrational partition sum of the asymmetric stretch mode
        :return: 
        """
        vibrational_partition_sum = 0

        for v3 in range(self._v3_max):
            vibrational_partition_sum += self._get_unnormalized_asymmetric_stretch_vibrational_fractional_population(
                models.TriatomicLinearFermiResonantGlobalQuanta(
                    vibrational_quantum_1=0,
                    vibrational_quantum_2=0,
                    vibrational_quantum_3=v3,
                    vibrational_angular_momentum=0,
                    vibrational_ranking_index=1
                )
            )

        return vibrational_partition_sum

    @cache.cache_by_global_quanta_decorator
    def get_asymmetric_stretch_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the fractional population of a state in the asymmetric stretch, described by the given global quanta, such 
        that the sum of the fractional populations of all asymmetric stretch states equals unity. To retrieve the number 
        density in a given bending mode state, this fractional population has to be multiplied by the total number 
        density of the isotopologue.
        :param global_quanta: 
        :return: 
        """
        return self._get_unnormalized_asymmetric_stretch_vibrational_fractional_population(global_quanta) / self.get_asymmetric_stretch_vibrational_partition_sum()

    @cache.cache_by_global_quanta_decorator
    def get_vibrational_fractional_population(self, global_quanta: models.TriatomicLinearFermiResonantGlobalQuanta):
        """
        Get the vibrational fractional population of a given vibrational state described by its global quanta, by 
        multiplying the vibrational fractional populations of each of the vibrational modes separately.
        :param global_quanta: 
        :return: 
        """
        return self.get_symmetric_stretch_resonant_vibrational_fractional_population(global_quanta) * \
            self.get_bending_mode_non_resonant_vibrational_fractional_population(global_quanta) * \
            self.get_asymmetric_stretch_vibrational_fractional_population(global_quanta)


class NonEquilibriumBoltzmannPyramidalTetratomicFractionDistribution(EquilibriumDistribution,
                                                                     NonCharacterisableEnergyLevels):
    vibrational_energies = {}

    def __init__(self, statistical_weight: weights.RotationalStatisticalWeight, rotational_constants,
                 vibrational_constants, *,
                 max_total_rotational_angular_momentum_qn=40, max_stretch_mode_qn=5, max_bending_mode_qn=10):
        super().__init__()

        self._statistical_weight = statistical_weight
        self._rotational_constants = rotational_constants
        self._omega_symmetric_stretch = vibrational_constants[0]
        self._constants_symmetric_bend = vibrational_constants[1]
        self._omega_asymmetric_stretch = vibrational_constants[2]
        self._omega_asymmetric_bend = vibrational_constants[3]

        self._max_J = max_total_rotational_angular_momentum_qn + 1
        self._max_v1_v3 = max_stretch_mode_qn + 1
        self._max_v2_v4 = max_bending_mode_qn + 1

        self._vibrational_energies = exomol.ExoMolDatabase(
            r"C:\Users\steij\Documents\Promotie\data\database\raw-nh3\14N-1H3__CoYuTe-pure_vib_levels.states")

        self.temperature_vib_K = 0

    @cache.reset(check_parameter_change=True)
    def set_temperature_K(self, rotational_temperature_kelvin: float, vibrational_temperature_kelvin: float):
        self.temperature_rot_K = rotational_temperature_kelvin
        self.temperature_vib_K = vibrational_temperature_kelvin

    @cache.cache_by_global_quanta_decorator
    def get_vibrational_energy(self, global_quanta: models.PyramidalTetratomicGlobalQuanta):
        return self._vibrational_energies.get(global_quanta)

    def get_rotational_fractional_population(self, global_quanta, total_energy: float):
        return np.exp(-constants.c2 / self.temperature_rot_K * (total_energy - self.get_vibrational_energy(global_quanta)))

    @cache.cache_fixed_value
    def get_partition_sum(self):
        """
            The degeneracy of a ro-vibrational level depends on both the vibrational and rotational quanta numbers and
            symmetries. This is true to the extent that the partition function cannot be easily split in a vibrational
            and rotational contribution, as is done for CO and CO2.
        """
        beta_rot, beta_vib = constants.c2 / self.temperature_rot_K, constants.c2 / self.temperature_vib_K

        @functools.cache
        def rotational_fractional_density(total_rotational_qn: int, z_projection_rotational_qn: int):
            _rotational_energy = .0
            for i, rotational_constants_array in enumerate(self._rotational_constants):
                for j, rotational_constant in enumerate(rotational_constants_array):
                    _rotational_energy += rotational_constant * \
                                          (total_rotational_qn * (total_rotational_qn + 1)) ** i * \
                                          z_projection_rotational_qn ** (2 * j)
            return (2 * total_rotational_qn + 1) * np.exp(-beta_rot * _rotational_energy)

        def energy_v1(symmetric_stretch_quantum_number: int):
            """
                First order anhamornic approximation is used to describe the energy levels of this normal mode.
            :param symmetric_stretch_quantum_number:
            :return:
            """
            return self._omega_symmetric_stretch[0] * symmetric_stretch_quantum_number + \
                self._omega_symmetric_stretch[1] * symmetric_stretch_quantum_number ** 2

        @functools.cache
        def energy_v2(symmetric_bending_qn: int, inversion_parity: str):
            """
                Some levels deviate from the straightforward 1st order anharmonic approximation, which is related to
                the presence of the potential energy barrier for sigma=a (as the discrepancy lowers after v2=4)
            :param symmetric_bending_qn:
            :param inversion_parity:
            :return:
            """
            if (inversion_parity == 's') & (symmetric_bending_qn < 6):
                return [0, 932.433472, 1597.475125, 2384.153000, 3462.410069, 4694.400808][symmetric_bending_qn]
            elif (symmetric_bending_qn, inversion_parity) == (1, 'a'):
                return 968.12193
            elif inversion_parity == 's':
                return sum([constant * symmetric_bending_qn ** i for i, constant in
                            enumerate(self._constants_symmetric_bend['s'])])
            elif inversion_parity == 'a':
                return sum([constant * symmetric_bending_qn ** i for i, constant in
                            enumerate(self._constants_symmetric_bend['a'])])
            raise ValueError("inversion parity has to be either symmetric or asymmetric")

        @functools.cache
        def energy_v3(asymmetric_stretch_qn: int):
            """
                First order anharmonic approximation satisfies characterising the energy levels for l3 = 0, similar to
                the symmetric stretch mode.
                Also, the change in the energy as result of the angular momentum l3 is neglected.
            :param asymmetric_stretch_qn:
            :return:
            """
            return self._omega_asymmetric_stretch[0] * asymmetric_stretch_qn + \
                self._omega_asymmetric_stretch[1] * asymmetric_stretch_qn ** 2

        @functools.cache
        def energy_v4(asymmetric_bending_qn: int):
            return self._omega_asymmetric_bend * asymmetric_bending_qn

        @functools.cache
        def degeneracy_l3_and_l4_is_0(is_total_rotational_quantum_number_0: bool, z_projection_rotational_quantum: int,
                                      inversion_parity: str) -> int:
            if z_projection_rotational_quantum == 0:
                if is_total_rotational_quantum_number_0:
                    if inversion_parity == 'a':
                        return 2
                elif inversion_parity == 's':
                    return 2
                return 0
            elif z_projection_rotational_quantum % 3 == 0:
                return 2
            return 1

        def structure_degeneracy_l_not_0(K_is_0: int, modulo_K_with_K_is_0: int, module_K_with_3_is_not_0: int):
            def degeneracy_function(z_projection_rotational_angular_momentum):
                if z_projection_rotational_angular_momentum == 0:
                    return K_is_0
                elif z_projection_rotational_angular_momentum % 3 == 0:
                    return modulo_K_with_K_is_0
                elif z_projection_rotational_angular_momentum % 3 != 0:
                    return module_K_with_3_is_not_0
                raise ValueError("Wut #1")

            return degeneracy_function

        @functools.cache
        def degeneracy_l_not_0(stretch_angular_qn, bending_angular_qn):
            if (stretch_angular_qn % 3 != 0) & (bending_angular_qn % 3 != 0):
                return structure_degeneracy_l_not_0(3, 6, 5)
            elif (stretch_angular_qn % 3 != 0) ^ (bending_angular_qn % 3 != 0):
                if (stretch_angular_qn == 0) ^ (bending_angular_qn == 0):
                    return structure_degeneracy_l_not_0(1, 2, 3)
                return structure_degeneracy_l_not_0(2, 4, 6)
            elif (stretch_angular_qn % 3 == 0) & (bending_angular_qn % 3 == 0):
                return structure_degeneracy_l_not_0(2, 4, 2)
            raise ValueError("Wut #2")

        partition_sum_symmetric_stretch_mode = sum([np.exp(-beta_vib * energy_v1(v1)) for v1 in range(self._max_v1_v3)])

        # first for l3 = l4 = 0:
        partition_sum_l_is_0 = 0
        for J in range(self._max_J):
            for K in range(J + 1):
                for parity in ['s', 'a']:
                    partition_sum_l_is_0 += \
                        sum([degeneracy_l3_and_l4_is_0(J % 2 == 0, K, parity) * rotational_fractional_density(J, K) *
                             np.exp(-beta_vib * energy_v2(v2, parity)) for v2 in range(self._max_v2_v4)])
        partition_sum_l_is_0 *= \
            sum([np.exp(-beta_vib * energy_v3(v3)) for v3 in range(0, self._max_v1_v3, 2)]) * \
            sum([np.exp(-beta_vib * energy_v4(v4)) for v4 in range(0, self._max_v2_v4, 2)])

        # second for l3 != 0 ^ l4 != 0:
        partition_sum_l_not_0 = 0
        for v3 in range(self._max_v1_v3):
            for v4 in range(self._max_v2_v4):
                vibrational_angular_term = 0
                for l3 in range(v3, -1, -2):
                    for l4 in range(v4, -1, -2):
                        if (l3 == 0) & (l4 == 0):
                            continue                # skipping as it is already considered
                        for J in range(self._max_J):
                            vibrational_angular_term += sum(
                                [degeneracy_l_not_0(l3, l4)(K) * rotational_fractional_density(J, K)
                                 for K in range(J + 1)]
                            )
                partition_sum_l_not_0 += vibrational_angular_term * np.exp(-beta_vib * (energy_v3(v3) + energy_v4(v4)))
        partition_sum_l_not_0 *= sum([
            sum([np.exp(-beta_vib * energy_v2(v2, parity)) for parity in ['s', 'a']])
            for v2 in range(self._max_v2_v4)])

        return self._statistical_weight.get_state_independent_weight(None) * partition_sum_symmetric_stretch_mode * (
            partition_sum_l_is_0 + partition_sum_l_not_0
        )

    def get_vibrational_fractional_population(self, global_quanta: models.PyramidalTetratomicGlobalQuanta):
        return np.exp(-constants.c2 * self.get_vibrational_energy(global_quanta) / self.temperature_vib_K)

    def get_fractional_population(self, global_quanta: models.PyramidalTetratomicGlobalQuanta,
                                  local_quanta: models.PyramidalTetratomicLocalQuanta,
                                  total_energy: float) -> float:
        return self._statistical_weight.get_rotational_weight(global_quanta, local_quanta) * \
            self.get_rotational_fractional_population(global_quanta, total_energy) * \
            self.get_vibrational_fractional_population(global_quanta) / self.get_partition_sum()