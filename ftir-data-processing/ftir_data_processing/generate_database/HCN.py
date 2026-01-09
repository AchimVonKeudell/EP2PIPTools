from context import *

if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name='HCN',
        afgl_codes=[124, 134, 125],
        file_names_partition_sum=[f"{database_directory}\\hcn\\q{i}.txt" for i in [70, 71, 72]],
        parsum_path=f"{database_directory}\\hcn\\raw"
    )

    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\hcn_database.sqlite",
        hitran_directory=f"{database_directory}\\hcn\\raw"
    )
