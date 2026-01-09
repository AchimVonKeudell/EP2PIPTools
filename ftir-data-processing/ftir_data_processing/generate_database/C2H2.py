from context import *

if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name="C2H2",
        afgl_codes=[1221, 1231, 1222],
        file_names_partition_sum=[f"{database_directory}\\c2h2\\q%i.txt" % i for i in [76, 77, 105]],
        parsum_path=f"{database_directory}\\c2h2\\raw"
    )

    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\c2h2_database.sqlite",
        hitran_directory=f"{database_directory}\\c2h2\\raw"
    )
