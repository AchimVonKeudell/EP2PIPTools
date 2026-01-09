import sys
from spectrumsimulator.database import hitran

if len(sys.argv) < 3:
    raise Exception("The location of the HITRAN sqlite file and the HITRAN directory should be set as a command line arguments.")
    sys.exit(1)

hitran_sqlite_file = sys.argv[1]
hitran_directory = sys.argv[2]
hitran.HitranDatabase.generate_database_linear_molecules(hitran_sqlite_file, hitran_directory)
# hitran.HitranDatabase.generate_database_pyramidal_tetratomic_molecules(hitran_sqlite_file, hitran_directory)
