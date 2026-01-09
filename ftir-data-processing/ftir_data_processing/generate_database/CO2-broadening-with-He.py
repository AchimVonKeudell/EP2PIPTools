from context import *



if __name__ == "__main__":

    # # configure the parsum file
    # parsum.GeneratingParsumFile(
    #     molecule_name="CO2",
    #     afgl_codes=[626, 636, 628, 627, 638, 637, 828, 827, 727, 838, 837, 737],
    #     file_names_partition_sum=[f"{database_directory}\\co2\\q%i.txt" % i for i in (7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122)],
    #     parsum_path=f"{database_directory}\\co2\\raw")
    #
    # # alter par file to include the He-broadening coefficients
    # reconfigure_hitran_par_file(
    #     original_hitran_par=f"{database_directory}\\co2\\raw\\02_HITRAN_data.par",
    #     he_broadening_coefficients_file=f"{database_directory}\\co2\\02_HITRAN_data-He.out",
    #     output_database=f"{database_directory}\\co2\\raw-He\\02_HITRAN_data-He.par"
    # )

    # generate the SQlite database
    hitran.HitranDatabase.generate_database_linear_molecules(
        database_path=f"{database_directory}\\co2_database-He.sqlite",
        hitran_directory=f"{database_directory}\\co2\\raw-He"
    )

