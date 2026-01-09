from context import hitran, parsum


if __name__ == "__main__":
    parsum.GeneratingParsumFile(
        molecule_name='H2O',
        afgl_codes=[161, 181, 171, 162, 182, 172, 262],
        file_names_partition_sum=[r"D:\Promotie\data\database\h2o\q%i.txt" % i for i in [1, 2, 3, 4, 5, 6, 129]],
        parsum_path=r'D:\Promotie\data\database\h2o\raw')
    hitran.HitranDatabase.generate_database_non_linear_molecules(
        database_path=r"D:\Promotie\data\database\h2o_database.sqlite",
        hitran_directory=r"D:\Promotie\data\database\h2o\raw"
    )
