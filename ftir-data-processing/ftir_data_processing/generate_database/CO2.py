from context import *


if __name__ == "__main__":
    # configure the parsum file
    parsum.GeneratingParsumFile(
        molecule_name="CO2",
        afgl_codes=[626, 636, 628, 627, 638, 637, 828, 827, 727, 838, 837, 737],
        file_names_partition_sum=[f"{database_directory}\\co2\\q%i.txt" % i for i in (7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122)],
        parsum_path=f"{database_directory}\\co2\\raw")

    # generate the SQlite database
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\co2_database.sqlite",
        hitran_directory=f"{database_directory}\\co2\\raw"
    )
