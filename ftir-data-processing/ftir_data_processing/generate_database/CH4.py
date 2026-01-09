from context import *


if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name="CH4",
        afgl_codes=[211, 311, 212, 312],
        file_names_partition_sum=[f"{database_directory}\\ch4\\q%i.txt" % i for i in range(32, 36)],
        parsum_path=f"{database_directory}\\ch4\\raw")
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\ch4_database.sqlite",
        hitran_directory=f"{database_directory}\\ch4\\raw"
    )
