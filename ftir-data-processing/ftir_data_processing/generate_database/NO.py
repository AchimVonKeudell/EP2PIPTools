from context import *


if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name="NO",
        afgl_codes=[46, 56, 48],
        file_names_partition_sum=[f"{database_directory}\\no\\q%i.txt" % i for i in (39, 40, 41)],
        parsum_path=f"{database_directory}\\no\\raw"
    )
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\no_database.sqlite",
        hitran_directory=f"{database_directory}\\no\\raw"
    )