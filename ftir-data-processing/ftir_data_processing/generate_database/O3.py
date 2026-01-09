from context import *

if __name__ == "__main__":
    # configure the parsum file
    parsum.GeneratingParsumFile(
        molecule_name="O3",
        afgl_codes=[666, 668, 686, 667, 676],
        file_names_partition_sum=[f"{database_directory}\\o3\\q%i.txt" % i for i in (16, 17, 18, 19, 20)],
        parsum_path=f"{database_directory}\\o3\\raw")

    # generate the SQlite database
    hitran.HitranDatabase.generate_database_non_linear_molecules(
        database_path=fr"{database_directory}\o3_database.sqlite",
        hitran_directory=fr"{database_directory}\o3\raw"
    )
