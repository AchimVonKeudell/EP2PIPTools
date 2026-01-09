from context import *


if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name="N2O",
        afgl_codes=[446, 456, 546, 448, 447],
        file_names_partition_sum=[f"{database_directory}\\n2o\\q%i.txt" % i for i in (21, 22, 23, 24, 25)],
        parsum_path=f"{database_directory}\\n2o\\raw"
    )
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\n2o_database.sqlite",
        hitran_directory=f"{database_directory}\\n2o\\raw"
    )