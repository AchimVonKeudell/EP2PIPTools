from dataclasses import dataclass
try:
    from ..vibrationalspectrum.spectrumsimulator.database import hitran, parsum, readers
except ImportError:
    # Fallback for direct execution
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from vibrationalspectrum.spectrumsimulator.database import hitran, parsum, readers


database_directory = "D:\\Promotie\\data\\database"


"""This script uses the normal line-by-line .par data of HITRAN and a custom database consisting of He-broadening 
parameters to switch the air broadening and pressure shift parameters with that of He, to better represent the actual
line profile.
"""
HITRAN_MAPPER = readers.HITRAN160ParameterMapper()
HITRAN_He_MAPPER = readers.HITRANHeBroadenedParameterMapper()


def reconfigure_hitran_par_file(original_hitran_par: str, he_broadening_coefficients_file: str, output_database: str):

    @dataclass
    class Index:
        molecule_id: int
        isotopologue_order_num: int
        wavenumber: float

        def __hash__(self):
            return hash((self.molecule_id, self.isotopologue_order_num, self.wavenumber))

        def __eq__(self, other):
            return (self.molecule_id, self.isotopologue_order_num, self.wavenumber) == \
                (other.molecule_id, other.isotopologue_order_num, other.wavenumber)

        def __ne__(self, other):
            return not (self, other)

    if (not os.path.exists(original_hitran_par)) | (not os.path.exists(he_broadening_coefficients_file)):
        raise FileExistsError(f"{os.path.exists(original_hitran_par)=}, "
                              f"{os.path.exists(he_broadening_coefficients_file)=}")

    # 1. read the original HITRAN database
    hitran_data = {}
    with open(original_hitran_par, 'r') as hitran_original_database:
        for hitran_line in hitran_original_database:
            acquired_data = HITRAN_MAPPER.get_fields_from_string(hitran_line)
            line_identifier = Index(
                acquired_data['molecule_id'],
                acquired_data['isotopologue_order_num'],
                acquired_data['wavenumber_vacuum']
            )
            hitran_data[line_identifier] = acquired_data

    # 2. read alternative pressure broadening and shift constants
    broadening_data = {}
    with open(he_broadening_coefficients_file, 'r') as broadening_database:
        for broadening_line in broadening_database:
            acquired_data = HITRAN_He_MAPPER.get_fields_from_string(broadening_line)
            line_identifier = Index(
                acquired_data['molecule_id'],
                acquired_data['isotopologue_order_num'],
                acquired_data['wavenumber_vacuum']
            )
            broadening_data[line_identifier] = acquired_data

    # 3. update HITRAN database:
    for key in hitran_data:

        try:
            hitran_data[key]['broadening_hw_air'] = broadening_data[key]['broadening_hw_He']
        except KeyError as e:
            print(
                e)  # testing if the specific ro-vibrational transition is also found in the alternative broadening database
            continue

        if broadening_data[key]['pressure_shift']:
            hitran_data[key]['pressure_shift'] = broadening_data[key]['pressure_shift']
        hitran_data[key]['broadening_temp_coefficient'] = broadening_data[key]['broadening_temp_coefficient']

    # 4. save updated HITRAN database
    with open(output_database, 'w') as new_database:
        for key, item in hitran_data.items():
            new_database.write(HITRAN_MAPPER.get_string_from_fields(item) + '\n')
    return


__all__ = ['database_directory', 'parsum', 'hitran', 'reconfigure_hitran_par_file']