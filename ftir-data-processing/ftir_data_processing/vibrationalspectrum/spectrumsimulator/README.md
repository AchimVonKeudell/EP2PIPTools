# Generate the database:
The project utilises a SQLite database to calculate the IR spectra of e.g. CO2, CO, and NH3. 
This ReadMe file explains how to configure the required database and explains how to alter the project, e.g. to include another specie,
which is done by listing the files/objects that must be altered.
Note that this file does not explain how to use this project in general, since the comments are too specific/technical.


## 1. How to generate the SQLite database:
The spectroscopic constants of the HITRAN database are needed to simulate the ro-vibrational spectrum of a specie. 
In this section, the usage for the function that can reconfigure the extensive spectroscopic data into a SQ-lite data is
described. The following function is to be used,

    database.hitran.HitranDatabase.generate_database(database_path, hitran_directory, hitran_parameter_mapper)

where,
* `database_path`: name of the .sqlite file
* `hitran_directory`: directory where the HITRAN/HTEMP files are stored
* `hitran_parameter_mapper` (optional, default=[readers.HITRAN160ParameterMapper](database/readers.py)): mapper used for recognising the HITRAN/HITEMP data
  * for HITEMP database of CO and CO2: [readers.HITRANCO2BroadenedParameterMapper](database/readers.py)

Before you call this function, you need to gather the following files into the `hitran_directory`:
1. the HITRAN/HITEMP data in `.par`-format
2. molecular parameter file: [molparam.txt](https://hitran.org/media/molparam.txt)
3. the partition functions of the various isotopologues in a `parsum.dat` file, using [parsum.py](database/parsum.py)



## 2. How to add/alter a global+local quantum structure:
Simulating ro-vibrational spectra in equilibrium is rather straightforward, by including the function
`GetSimulator.get_M_simulator()` for a specie `M`. If you want to simulate non-equilibrium spectra, then you need to implement
the following changes, so that this project knows how to describe the rotational+vibrational quantum numbers and how to 
read these quantum numbers from the HITRAN database `.par` file.

Add/alter functions, objects, and variables in:
1. In [queries.py](database/queries.py), add/adjust queries: `structure_.._table_create`, `line_.._table_create`, `structure_.._insert_query`, and `line_.._insert_query`
2. In [init file](__init__.py), add the required molecule_id for ID'ing the global_quanta model. 
3. In [models](database/models.py), 
   1. Create related object `models.GlobalQuanta` and add it to the if-statement checks in `models.Structure.get_global_quanta_model_type_by_molecule_id()`
   2. Create related object `models.LocalQuanta` and add it to the if-statement checks in `models.Line.get_local_quanta_model_type_by_molecule_id()`
4. In [repositories.py](database/repositories.py), create `Structure..Repository` and `Line..Repository` objects, to consider the global quanta properties of this specific molecule type
5. In [hitran.py](hitran.py), create a new object based on `hitran.HitranDatabaseGenerator`, to consider properties of this specific molecule, e.g. see `HitranLinearMoleculesDatabaseGenerator` as reference
6. In [results.py](results.py), 
   1. add `record_to_structure_model_..` and `record_to_line_model_..` functions to classify the HITRAN/HITEMP data
   2. implement these in the objects `results.StructureResultSet` and `results.LineResultSet`, respectively

When testing the code, the following error could mean:
* originating from `repo_structures.get_filtered(..)` -> `cursor.execute(str(query)) .. Molecules not found in Table`: sqlite database does not exits, wrong name?
* originating from `repo_structures.get_filtered(..)` -> `cursor.execute(str(query)) ..`: structure Table is wrongly defined in [repositories.StructureRepository](repositories.py).


## 3. How to add a distribution object:
See introduction of previous section, i.e. only required when including non-equilibrium spectrum

Distribution object must include:
* initialise function: `rotational_statistical_weights`, `molecular_constants` (if applicable)
* function `get_fractional_population(global_quanta, local_quanta)`
* equilibrium temperature: `temperature_rot_K`
* non-equilibrium temperature: `temperature_vib_K`
