from context import *


if __name__ == "__main__":

    # configure the parsum file
    parsum.GeneratingParsumFile(
        molecule_name="NH3",
        afgl_codes=[4111, 5111],
        file_names_partition_sum=[f"{database_directory}\\nh3\\q{i}.txt" for i in (45, 46)],
        parsum_path=f"{database_directory}\\nh3\\raw"
    )

    # generate the SQlite database
    hitran.HitranDatabase.generate_database_pyramidal_tetratomic_molecules(
        database_path=f"{database_directory}\\nh3_database.sqlite",
        hitran_directory=f"{database_directory}\\nh3\\raw"
    )



