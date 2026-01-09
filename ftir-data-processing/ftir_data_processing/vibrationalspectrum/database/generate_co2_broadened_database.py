import sys
from spectrumsimulator.database import hitran, readers


if len(sys.argv) < 3:
    raise Exception("The location of the HITRAN sqlite file and the HITRAN directory should be set as a command line arguments.")
    sys.exit(1)

hitran_sqlite_file = sys.argv[1]            # the path + file name of the SQlite database
hitran_directory = sys.argv[2]              # directory where the HITEMP/HITRAN files are located
hitran.HitranDatabase.generate_database_linear_molecules(
    hitran_sqlite_file, hitran_directory, hitran_parameter_mapper=readers.HITRANCO2BroadenedParameterMapper())
