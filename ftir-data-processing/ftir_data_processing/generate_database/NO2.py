from context import *


if __name__ == "__main__":

    # configure the parsum file
    parsum.GeneratingParsumFile(
        molecule_name="NO2",
        afgl_codes=[646, 656],
        file_names_partition_sum=[f"{database_directory}\\no2\\q%i.txt" % i for i in (44, 130)],
        parsum_path=f"{database_directory}\\no2\\raw"
    )
    # generate the SQlite database
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\no2_database.sqlite",
        hitran_directory=f"{database_directory}\\no2\\raw"
    )