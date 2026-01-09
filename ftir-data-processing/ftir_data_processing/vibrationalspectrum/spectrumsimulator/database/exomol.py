import os
from dataclasses import dataclass
from . import models

symmetry = {'1': 'A1\'', '2': 'A2\'', '3': 'E\'', '4': 'A1\"', '5': 'A2\"', '6': 'E\"'}
parity = {'0': 's', '1': 'a'}


class ExoMolMapper:
    @dataclass
    class ExoMolAmmoniaDataStructure:
        state_counting_number: int
        energy: float
        total_state_degeneracy: int
        total_angular_momentum_1: int
        total_state_parity: str
        total_symmetry: str
        block_state_counting_number: int
        symmetric_stretch_qn: int
        symmetric_bend_qn: int
        asymmetric_stretch_qn: int
        asymmetric_bend_qn: int
        asymmetric_stretch_angular_qn: int
        asymmetric_bend_angular_qn: int
        inversion_parity: str
        total_angular_momentum: int
        z_projection_angular_momentum: int
        rotational_symmetry: str
        local_mode_vibrational_quantum_numbers: str
        vibrational_symmetry: str
        energy_theoretical: float

    @staticmethod
    def ammonia_14NH3(exo_mol_line_data: str):
        return ExoMolMapper.ExoMolAmmoniaDataStructure(
            state_counting_number=int(exo_mol_line_data[:12]),
            energy=float(exo_mol_line_data[12:25]),
            total_state_degeneracy=int(exo_mol_line_data[25:32]),
            total_angular_momentum_1=int(exo_mol_line_data[32:40]),
            total_state_parity=exo_mol_line_data[40:55],
            total_symmetry=exo_mol_line_data[55:58].strip(),
            block_state_counting_number=int(exo_mol_line_data[58:69]),
            symmetric_stretch_qn=int(exo_mol_line_data[69:76]),
            symmetric_bend_qn=int(exo_mol_line_data[76:80]),
            asymmetric_stretch_qn=int(exo_mol_line_data[80:84]),
            asymmetric_bend_qn=int(exo_mol_line_data[84:88]),
            asymmetric_stretch_angular_qn=int(exo_mol_line_data[88:94]),
            asymmetric_bend_angular_qn=int(exo_mol_line_data[94:98]),
            inversion_parity=exo_mol_line_data[98:103].strip(),
            total_angular_momentum=int(exo_mol_line_data[103:110]),
            z_projection_angular_momentum=int(exo_mol_line_data[110:114]),
            rotational_symmetry=exo_mol_line_data[114:117].strip(),
            local_mode_vibrational_quantum_numbers=exo_mol_line_data[117:145],
            vibrational_symmetry=exo_mol_line_data[145:151].strip(),
            energy_theoretical=float(exo_mol_line_data[151:-1])
        )


class ExoMolLineData:

    def __init__(self, exomol_line_data: str):

        _data = ExoMolMapper.ammonia_14NH3(exomol_line_data)

        self.energy = _data.energy
        self.global_quanta = models.PyramidalTetratomicGlobalQuanta(
            vibrational_quantum_1=_data.symmetric_stretch_qn,
            vibrational_quantum_2=_data.symmetric_bend_qn,
            vibrational_quantum_3=_data.asymmetric_stretch_qn,
            vibrational_quantum_4=_data.asymmetric_bend_qn,
            vibrational_angular_momentum_3=_data.asymmetric_stretch_angular_qn,
            vibrational_angular_momentum_4=_data.asymmetric_bend_angular_qn,
            vibrational_angular_momentum=None,
            vibrational_symmetry=symmetry.get(_data.vibrational_symmetry),
        )

        self.local_quanta = models.PyramidalTetratomicLocalQuanta(
            J=_data.total_angular_momentum,
            K=_data.z_projection_angular_momentum,
            inversion_symmetry=parity.get(_data.inversion_parity),
            symmetry_rotational=symmetry.get(_data.rotational_symmetry),
            symmetry_total=symmetry.get(_data.total_symmetry),
        )


class ExoMolDatabase:

    def __init__(self, database_full_path: str):
        if not os.path.exists(database_full_path):
            raise FileExistsError(f"{database_full_path=} does not exist")

        self._path = database_full_path
        self.structures = []

        self.acquire_database()

    def acquire_database(self):
        with open(self._path, 'r') as file_data:
            for line_data in file_data:
                self.structures.append(ExoMolLineData(line_data))
        return

    def get(self, global_quanta: models.PyramidalTetratomicGlobalQuanta) -> float:
        for structure in self.structures:
            if structure.global_quanta == global_quanta:
                return structure.energy
        raise ValueError(f'global quanta not found: {global_quanta}')

    def get_local_state(self, global_quanta: models.PyramidalTetratomicGlobalQuanta,
                        local_quanta: models.PyramidalTetratomicLocalQuanta):
        for structure in self.structures:
            if structure.global_quanta == global_quanta:
                if structure.local_quanta == local_quanta:
                    return structure.energy

