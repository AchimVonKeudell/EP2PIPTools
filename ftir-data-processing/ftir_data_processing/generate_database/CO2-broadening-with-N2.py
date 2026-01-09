from context import *
from context import HITRAN_MAPPER


def reconfigure_hitran_par_file(original_hitran_par: str, output_database: str, correcting_factor: float):

    hitran_data = []

    with open(original_hitran_par, 'r') as hitran_original_database:
        for hitran_line in hitran_original_database:
            acquired_data = HITRAN_MAPPER.get_fields_from_string(hitran_line)

            # adjust the broadening from air to N2
            acquired_data['broadening_hw_air'] *= correcting_factor
            hitran_data.append(acquired_data)


    with open(output_database, 'w') as new_database:
        for y in hitran_data:
            new_database.write(HITRAN_MAPPER.get_string_from_fields(y) + '\n')
    return


if __name__ == "__main__":
    # adjusts the HITRAN database
    reconfigure_hitran_par_file(
        original_hitran_par=f"{database_directory}\\co2\\raw\\02_HITRAN_data.par",
        output_database=f"{database_directory}\\co2\\raw-N2\\02_HITRAN_data-N2.par",
        correcting_factor=1-0.028
    )

    # configure the parsum file
    parsum.GeneratingParsumFile(
        molecule_name='CO2',
        afgl_codes=[626, 636, 628, 627, 638, 637, 828, 827, 727, 838, 837, 737],
        file_names_partition_sum=[f"{database_directory}\\co2\\q%i.txt" % i for i in (7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122)],
        parsum_path=f'{database_directory}\\co2\\raw-N2')

    # generate the SQlite database
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\co2_database-N2.sqlite",
        hitran_directory=f"{database_directory}\\co2\\raw-N2"
    )
