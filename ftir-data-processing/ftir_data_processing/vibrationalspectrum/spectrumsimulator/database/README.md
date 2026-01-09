# Generate the database:
The project utilises a SQLite database to calculate the IR spectra of CO2, CO, and NH3. 
This ReadMe file explains how to configure the required database and explains how to alter the project, e.g. to include another specie,
which is done by listing the files/objects that must be altered.

## How to generate the SQLite database:

``database.hitran.HitranDatabase.generate_database(database_path, hitran_directory, hitran_parameter_mapper)``

* `database_path`: name of the .sqlite file
* `hitran_directory`: directory where the HITRAN/HTEMP files are stored
* `hitran_parameter_mapper` (optional, default=[readers.HITRAN160ParameterMapper](readers.py)): mapper used for recognising the HITRAN/HITEMP data
  * for HITEMP database of CO and CO2: [readers.HITRANCO2BroadenedParameterMapper](readers.py)


Requirements to store in `hitran_directory`:
- the HITRAN/HITEMP data in `.par`-format
- molecular parameter file: [molparam.txt](https://hitran.org/media/molparam.txt)
- the partition functions of the various isotopologues in a `parsum.dat` file, using [parsum.py](parsum.py)



## How to add/alter a global+local quantum structure:
Add/alter functions, objects, and variables in:
1. In [queries.py](queries.py), add/adjust queries: `structure_.._table_create`, `line_.._table_create`, `structure_.._insert_query`, and `line_.._insert_query`
2. In [repositories.py](repositories.py), create `Structure..Repository` and `Line..Repository` objects, to consider the global quanta properties of this specific molecule type
3. In [models](models.py), 
   1. Create objects based on object `models.GlobalQuanta` and `models.LocalQuanta` 
   2. alter `models.Structure.get_global_quanta_model_type_by_molecule_id()` and `models.Line.get_local_quanta_model_type_by_molecule_id()`
5. In [init file](__init__.py), add the required molecule_id for ID'ing the global_quanta model. 
6. In [hitran.py](hitran.py), create a new object based on `hitran.HitranDatabaseGenerator`, to consider properties of this specific molecule, e.g. see `HitranLinearMoleculesDatabaseGenerator` as reference
7. In [results.py](results.py), 
   1. add `record_to_structure_model_..` and `record_to_line_model_..` functions to classify the HITRAN/HITEMP data
   2. implement these in the objects `results.StructureResultSet` and `results.LineResultSet`, respectively

When testing the code, the following error could mean:
* origniating from `repo_structures.get_filtered(..)` -> `cursor.execute(str(query)) .. Molecules not found in Table`: sqlite database does not exits, wrong name?
* origniating from `repo_structures.get_filtered(..)` -> `cursor.execute(str(query)) ..`: structure Table is wrongly defined in [repositories.StructureRepository](repositories.py).


## How to add a distribution object:

Distribution object must include:
* initialise function: `rotational_statistical_weights`, `molecular_constants` (if applicable)
* function `get_fractional_population(global_quanta, local_quanta)`
* equilibrium temperature: `temperature_rot_K`
* non-equilibrium temperature: `temperature_vib_K`
