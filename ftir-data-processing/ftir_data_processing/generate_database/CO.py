from context import *


if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name='CO',
        afgl_codes=[26, 36, 28, 27, 38, 37],
        file_names=[f"{database_directory}\\co\\q%i.txt" % i for i in range(26, 32)],
        parsum_path=f'{database_directory}\\co\\raw')

    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\co_database.sqlite",
        hitran_directory=f"{database_directory}\\co\\raw"
    )
