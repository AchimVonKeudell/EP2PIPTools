"""
This query creates the molecule table, which consists of a molecule_id and the name of the molecule.
"""
molecule_table_create = """
    DROP TABLE IF EXISTS molecules;
    CREATE TABLE molecules (
        molecule_id INTEGER PRIMARY KEY, 
        molecule_name VARCHAR(6) NOT NULL
    );"""

"""
This query creates the isotopologue table, consisting of the molecule_id it belongs to, the order number of the
isotopologue, where number 1 is assigned to the most abundant isotopologue, number 2 to the subsequent most abundant
molecule, and so on. For each isotopologue, the AFGL (Air Force Geophysics Laboratory) shorthand notation, abundance 
fraction and molar mass in g/mol are stored.
"""
isotopologue_table_create = """
    DROP TABLE IF EXISTS isotopologues;
    CREATE TABLE isotopologues (
        molecule_id INTEGER NOT NULL, 
        isotopologue_order_num INTEGER NOT NULL,
        afgl_code VARCHAR(12) NOT NULL,
        abundance REAL NOT NULL,
        molar_mass REAL NOT NULL,
        PRIMARY KEY (molecule_id, isotopologue_order_num),
        FOREIGN KEY (molecule_id) REFERENCES molecules
    );
    CREATE INDEX isotopologues_molecule_id_idx ON isotopologues(molecule_id);"""

"""
Absorption lines are the result of a change in global (vibrational) quantum numbers and a change in local (rotational)
quantum numbers, where many local transitions belong to a single global transition. Therefore, they are grouped in 
so-called absorption structures. These global (vibrational) transitions are stored in a separate table, where lines 
are assigned to each of these structures.

For each structure, we store a unique id, the molecule and isotopologue it belongs to, the string representation of the
upper and lower global quanta and the maximum line strength of lines belonging to the transition.
"""
structure_table_create = """
    DROP TABLE IF EXISTS structures;
    CREATE TABLE structures (
        structure_id INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INTEGER NOT NULL, 
        isotopologue_order_num INTEGER NOT NULL,
        global_quanta_upper VARCHAR(15) NOT NULL,
        global_quanta_lower VARCHAR(15) NOT NULL,
        max_line_strength_296K REAL,
        min_wavenumber_vacuum REAL,
        max_wavenumber_vacuum REAL,
        FOREIGN KEY (molecule_id, isotopologue_order_num) REFERENCES isotopologues
    );
    CREATE INDEX structures_molecule_id_idx ON structures(molecule_id, isotopologue_order_num);"""

"""
There are different types of molecules that each have their own quantum numbers describing the global upper and lower
state. Therefore, there are different properties that have to be stored for each structure, depending on the type of 
molecule.
"""
structure_diatomic_table_create = """
    DROP TABLE IF EXISTS structures_diatomic;
    CREATE TABLE structures_diatomic (
        structure_id INTEGER PRIMARY KEY UNIQUE,
        vibrational_quantum_lower INTEGER,
        vibrational_quantum_upper INTEGER,
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX structures_diatomic_vibrational_quantum_number_lower_id_idx ON structures_diatomic(vibrational_quantum_lower);
    CREATE INDEX structures_diatomic_vibrational_quantum_number_upper_id_idx ON structures_diatomic(vibrational_quantum_upper);"""

structure_triatomic_linear_fermi_resonant_table_create = """
    DROP TABLE IF EXISTS structures_triatomic_linear_fermi_resonant;
    CREATE TABLE structures_triatomic_linear_fermi_resonant (
        structure_id INTEGER PRIMARY KEY UNIQUE,
        vibrational_quantum_1_lower INTEGER,
        vibrational_quantum_2_lower INTEGER,
        vibrational_quantum_3_lower INTEGER,
        vibrational_angular_momentum_lower INTEGER,
        vibrational_ranking_index_lower INTEGER,
        vibrational_quantum_1_upper INTEGER,
        vibrational_quantum_2_upper INTEGER,
        vibrational_quantum_3_upper INTEGER,
        vibrational_angular_momentum_upper INTEGER,
        vibrational_ranking_index_upper INTEGER,
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_1_lower_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_1_lower);
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_2_lower_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_2_lower);
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_3_lower_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_3_lower);
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_1_upper_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_1_upper);
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_2_upper_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_2_upper);
    CREATE INDEX structures_triatomic_linear_fermi_resonant_vibrational_quantum_number_3_upper_id_idx ON structures_triatomic_linear_fermi_resonant(vibrational_quantum_3_upper);"""

structure_triatomic_non_linear_table_create = """
    DROP TABLE IF EXISTS structures_triatomic_non_linear;
    CREATE TABLE structures_triatomic_non_linear (
        structure_id INTEGER PRIMARY KEY UNIQUE,
        vibrational_quantum_1_lower INTEGER,
        vibrational_quantum_2_lower INTEGER,
        vibrational_quantum_3_lower INTEGER,
        vibrational_quantum_1_upper INTEGER,
        vibrational_quantum_2_upper INTEGER,
        vibrational_quantum_3_upper INTEGER,
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_1_lower_id_idx ON structures_triatomic_non_linear(vibrational_quantum_1_lower);
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_2_lower_id_idx ON structures_triatomic_non_linear(vibrational_quantum_2_lower);
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_3_lower_id_idx ON structures_triatomic_non_linear(vibrational_quantum_3_lower);
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_1_upper_id_idx ON structures_triatomic_non_linear(vibrational_quantum_1_upper);
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_2_upper_id_idx ON structures_triatomic_non_linear(vibrational_quantum_2_upper);
    CREATE INDEX structures_triatomic_non_linear_vibrational_quantum_number_3_upper_id_idx ON structures_triatomic_non_linear(vibrational_quantum_3_upper);"""

structure_pyramidal_tetratomic_table_create = """
    DROP TABLE IF EXISTS structures_pyramidal_tetratomic;
    CREATE TABLE structures_pyramidal_tetratomic (
        structure_id INTEGER PRIMARY KEY UNIQUE,
        vibrational_quantum_1_lower INTEGER,
        vibrational_quantum_2_lower INTEGER,
        vibrational_quantum_3_lower INTEGER,
        vibrational_quantum_4_lower INTEGER,
        vibrational_angular_momentum_3_lower INTEGER,
        vibrational_angular_momentum_4_lower INTEGER,
        vibrational_angular_momentum_lower INTEGER,
        vibrational_symmetry_lower VARCHAR(5),
        vibrational_quantum_1_upper INTEGER,
        vibrational_quantum_2_upper INTEGER,
        vibrational_quantum_3_upper INTEGER,
        vibrational_quantum_4_upper INTEGER,
        vibrational_angular_momentum_3_upper INTEGER,
        vibrational_angular_momentum_4_upper INTEGER,
        vibrational_angular_momentum_upper INTEGER,
        vibrational_symmetry_upper VARCHAR(5),
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_1_lower_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_1_lower);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_2_lower_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_2_lower);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_3_lower_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_3_lower);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_4_lower_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_4_lower);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_1_upper_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_1_upper);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_2_upper_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_2_upper);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_3_upper_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_3_upper);
    CREATE INDEX structures_pyramidal_tetratomic_vibrational_quantum_number_4_upper_id_idx ON structures_pyramidal_tetratomic(vibrational_quantum_4_upper);"""

structure_pentatomic_table_create = """
    DROP TABLE IF EXISTS structures_pentatomic;
    CREATE TABLE structures_pentatomic (
        structure_id INTEGER PRIMARY KEY UNIQUE,
        vibrational_quantum_1_lower INTEGER,
        vibrational_quantum_2_lower INTEGER,
        vibrational_quantum_3_lower INTEGER,
        vibrational_quantum_4_lower INTEGER,
        vibrational_quantum_5_lower INTEGER,
        vibrational_quantum_6_lower INTEGER,
        vibrational_quantum_multiplicity_index_lower INTEGER,  
        vibrational_symmetry_lower VARCHAR(5),
        vibrational_quantum_1_upper INTEGER,
        vibrational_quantum_2_upper INTEGER,
        vibrational_quantum_3_upper INTEGER,
        vibrational_quantum_4_upper INTEGER,
        vibrational_quantum_5_upper INTEGER,
        vibrational_quantum_6_upper INTEGER,
        vibrational_quantum_multiplicity_index_upper INTEGER, 
        vibrational_symmetry_upper VARCHAR(5), 
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_1_lower_id_idx ON structures_pentatomic(vibrational_quantum_1_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_2_lower_id_idx ON structures_pentatomic(vibrational_quantum_2_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_3_lower_id_idx ON structures_pentatomic(vibrational_quantum_3_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_4_lower_id_idx ON structures_pentatomic(vibrational_quantum_4_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_5_lower_id_idx ON structures_pentatomic(vibrational_quantum_5_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_6_lower_id_idx ON structures_pentatomic(vibrational_quantum_6_lower);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_1_upper_id_idx ON structures_pentatomic(vibrational_quantum_1_upper);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_2_upper_id_idx ON structures_pentatomic(vibrational_quantum_2_upper);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_3_upper_id_idx ON structures_pentatomic(vibrational_quantum_3_upper);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_4_upper_id_idx ON structures_pentatomic(vibrational_quantum_4_upper);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_5_upper_id_idx ON structures_pentatomic(vibrational_quantum_5_upper);
    CREATE INDEX structures_pentatomic_vibrational_quantum_number_6_upper_id_idx ON structures_pentatomic(vibrational_quantum_6_upper);"""

"""
The lines table stores individual absorption line data, including a unique id per line, the structure to which it
belongs, the transition wavenumber (in cm-1) in vacuum, the line strength at 296K, the einstein A coefficient (in s-1)
air and selfbroadening halfwidths (in cm-1), the energy of the lower state (in cm-1) the pressure broadening exponential
temperature coefficient, the pressure shift (in cm-1/atm), the upper and lower local quanta string representations and
the rotational statistical weights for the upper and lower state (equal to 2J+1).
"""
line_table_create = """
    DROP TABLE IF EXISTS lines;
    CREATE TABLE lines (
        line_id INTEGER PRIMARY KEY AUTOINCREMENT,
        structure_id INTEGER NOT NULL,
        wavenumber_vacuum REAL NOT NULL, 
        line_strength_296K REAL NULL,
        einstein_A REAL NOT NULL,
        broadening_hw_air REAL NOT NULL,
        broadening_hw_self REAL NOT NULL,
        energy_lower_state REAL NOT NULL,
        broadening_temp_coefficient REAL NOT NULL,
        pressure_shift REAL NOT NULL,
        local_quanta_upper VARCHAR(15) NOT NULL,
        local_quanta_lower VARCHAR(15) NOT NULL,
        rotational_statistical_weight_upper REAL NOT NULL,
        rotational_statistical_weight_lower REAL NOT NULL,
        FOREIGN KEY (structure_id) REFERENCES structures
    );
    CREATE INDEX lines_structure_id_idx ON lines(structure_id);"""

"""
There are different types of molecules that each have their own local quantum numbers describing the local upper and 
lower state. Therefore, there are different properties that have to be stored for each lines, depending on the type of 
molecule.
"""
line_diatomic_or_linear_table_create = """
    DROP TABLE IF EXISTS lines_diatomic_or_linear;
    CREATE TABLE lines_diatomic_or_linear (
        line_id INTEGER PRIMARY KEY UNIQUE,
        J_lower INTEGER,
        symmetry_lower VARCHAR(1),
        F_lower VARCHAR(5),
        J_upper INTEGER,
        symmetry_upper VARCHAR(1),
        F_upper VARCHAR(5),
        FOREIGN KEY (line_id) REFERENCES lines
    );
    CREATE INDEX lines_diatomic_or_linear_J_lower_id_idx ON lines_diatomic_or_linear(J_lower);
    CREATE INDEX lines_diatomic_or_linear_J_upper_id_idx ON lines_diatomic_or_linear(J_upper);"""

line_triatomic_non_linear_table_create = """
    DROP TABLE IF EXISTS lines_triatomic_non_linear;
    CREATE TABLE lines_triatomic_non_linear (
        line_id INTEGER PRIMARY KEY UNIQUE,
        J_lower INTEGER,
        K_lower INTEGER,
        J_upper INTEGER,
        K_upper INTEGER,
        FOREIGN KEY (line_id) REFERENCES lines
    );
    CREATE INDEX lines_triatomic_non_linear_J_lower_id_idx ON lines_triatomic_non_linear(J_lower);
    CREATE INDEX lines_triatomic_non_linear_K_lower_id_idx ON lines_triatomic_non_linear(K_lower);
    CREATE INDEX lines_triatomic_non_J_upper_id_idx ON lines_triatomic_non_linear(J_upper);
    CREATE INDEX lines_triatomic_non_linear_K_upper_id_idx ON lines_triatomic_non_linear(K_upper);"""

line_pyramidal_tetratomic_table_create = """ 
    DROP TABLE IF EXISTS lines_pyramidal_tetratomic;
    CREATE TABLE lines_pyramidal_tetratomic (
        line_id INTEGER PRIMARY KEY UNIQUE,
        J_lower INTEGER,
        K_lower INTEGER,
        inversion_symmetry_lower VARCHAR(2),
        symmetry_rotational_lower VARCHAR(3),
        symmetry_total_lower VARCHAR(3),
        J_upper INTEGER,
        K_upper INTEGER,
        inversion_symmetry_upper VARCHAR(2),
        symmetry_rotational_upper VARCHAR(3),
        symmetry_total_upper VARCHAR(3),
        FOREIGN KEY (line_id) REFERENCES lines
    );
    CREATE INDEX lines_pyramidal_tetratomic_J_lower_id_idx ON lines_pyramidal_tetratomic(J_lower);
    CREATE INDEX lines_pyramidal_tetratomic_K_lower_id_idx ON lines_pyramidal_tetratomic(K_lower);
    CREATE INDEX lines_pyramidal_tetratomic_J_upper_id_idx ON lines_pyramidal_tetratomic(J_upper);
    CREATE INDEX lines_pyramidal_tetratomic_K_upper_id_idx ON lines_pyramidal_tetratomic(K_upper);"""

line_pentatomic_table_create = """
    DROP TABLE IF EXISTS lines_pentatomic;
    CREATE TABLE lines_pentatomic (
        line_id INTEGER PRIMARY KEY UNIQUE,
        J_lower INTEGER,
        symmetry_rotational_lower VARCHAR(3),
        alpha_lower INTEGER,
        F_lower INTEGER,
        J_upper INTEGER,
        symmetry_rotational_upper VARCHAR(3),
        alpha_upper INTEGER,
        F_upper INTEGER,
        FOREIGN KEY (line_id) REFERENCES lines
    );
    CREATE INDEX lines_pentatomic_J_lower_id_idx ON lines_pentatomic(J_lower);
    CREATE INDEX lines_pentatomic_alpha_lower_id_idx ON lines_pentatomic(alpha_lower);
    CREATE INDEX lines_pentatomic_J_upper_id_idx ON lines_pentatomic(J_upper);
    CREATE INDEX lines_pentatomic_alpha_upper_id_idx ON lines_pentatomic(alpha_upper);"""


"""
Partition sum data is stored for each isotopologue, for a range of temperature between 70K and 3000K.
"""

partition_sums_table_create = """
    DROP TABLE IF EXISTS partition_sums;
    CREATE TABLE partition_sums (
        molecule_id INTEGER NOT NULL, 
        isotopologue_order_num INTEGER NOT NULL,
        temperature REAL NOT NULL,
        partition_sum REAL NOT NULL,
        PRIMARY KEY (molecule_id, isotopologue_order_num, temperature),
        FOREIGN KEY (molecule_id, isotopologue_order_num) REFERENCES isotopologues
    );
    CREATE INDEX partition_sums_molecule_id_idx ON partition_sums(molecule_id, isotopologue_order_num);
    CREATE INDEX partition_sums_temperature_idx ON partition_sums(temperature);"""

""" 
Query to insert a molecule into the database
"""
molecule_insert = """
    INSERT INTO molecules
        (molecule_id, molecule_name)
    VALUES
        (?,?)"""

"""
Query to insert an isotopologue into the database
"""

isotopologue_insert = """
    INSERT INTO isotopologues
        (molecule_id, isotopologue_order_num, afgl_code, abundance, molar_mass)
    VALUES
        (?,?,?,?,?)"""

"""
Query to insert a structure into the database
"""
structure_insert_query = """
    INSERT INTO structures
        (molecule_id, isotopologue_order_num, global_quanta_upper, global_quanta_lower, max_line_strength_296K)
    VALUES
        (?,?,?,?,?)"""

"""
Query to insert additional data of a diatomic molecule for an existing structure into the database
"""
structure_diatomic_insert_query = """
    INSERT INTO structures_diatomic
        (structure_id, vibrational_quantum_lower, vibrational_quantum_upper)
    VALUES
        (?,?,?)"""

"""
Query to insert additional data of a linear triatomic molecule with a strong fermi resonance for an existing structure 
into the database
"""
structure_triatomic_linear_fermi_resonance_insert_query = """
    INSERT INTO structures_triatomic_linear_fermi_resonant
        (
            structure_id, 
            vibrational_quantum_1_lower, 
            vibrational_quantum_2_lower, 
            vibrational_quantum_3_lower,
            vibrational_angular_momentum_lower,
            vibrational_ranking_index_lower,
            vibrational_quantum_1_upper,
            vibrational_quantum_2_upper,
            vibrational_quantum_3_upper,
            vibrational_angular_momentum_upper,
            vibrational_ranking_index_upper
        )
    VALUES
        (?,?,?,?,?,?,?,?,?,?,?)"""
structure_triatomic_non_linear_insert_query = """
    INSERT INTO structures_triatomic_non_linear
        (
            structure_id, 
            vibrational_quantum_1_lower, 
            vibrational_quantum_2_lower, 
            vibrational_quantum_3_lower,
            vibrational_quantum_1_upper,
            vibrational_quantum_2_upper,
            vibrational_quantum_3_upper,
        )
    VALUES
        (?,?,?,?,?,?,?)"""
structure_pyramidal_tetratomic_insert_query = """
    INSERT INTO structures_pyramidal_tetratomic
        (
            structure_id, 
            vibrational_quantum_1_lower, 
            vibrational_quantum_2_lower, 
            vibrational_quantum_3_lower,
            vibrational_quantum_4_lower,
            vibrational_angular_momentum_3_lower,
            vibrational_angular_momentum_4_lower,
            vibrational_angular_momentum_lower,
            vibrational_symmetry_lower,
            vibrational_quantum_1_upper,
            vibrational_quantum_2_upper,
            vibrational_quantum_3_upper,
            vibrational_quantum_4_upper,
            vibrational_angular_momentum_3_upper,
            vibrational_angular_momentum_4_upper,
            vibrational_angular_momentum_upper,
            vibrational_symmetry_upper
        )
    VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""
structure_pentatomic_insert_query = """
    INSERT INTO structures_pentatomic
        (
            structure_id, 
            vibrational_quantum_1_lower, 
            vibrational_quantum_2_lower, 
            vibrational_quantum_3_lower,
            vibrational_quantum_4_lower,
            vibrational_quantum_5_lower,
            vibrational_quantum_6_lower,
            vibrational_quantum_multiplicity_index_lower,
            vibrational_symmetry_lower,
            vibrational_quantum_1_upper,
            vibrational_quantum_2_upper,
            vibrational_quantum_3_upper,
            vibrational_quantum_4_upper,
            vibrational_quantum_5_upper,
            vibrational_quantum_6_upper,
            vibrational_quantum_multiplicity_index_upper,
            vibrational_symmetry_upper
        )
    VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""

"""
Query to insert an absorption line into the database
"""
line_insert_query = """
    INSERT INTO lines
    (
        structure_id, 
        wavenumber_vacuum, 
        line_strength_296K, 
        einstein_A,
        broadening_hw_air,
        broadening_hw_self,
        energy_lower_state,
        broadening_temp_coefficient,
        pressure_shift,
        local_quanta_lower,
        local_quanta_upper,
        rotational_statistical_weight_lower,
        rotational_statistical_weight_upper
    ) 
    VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?)"""

"""
Query to insert additional data of a diatomic or linear molecule for an existing line into the database
"""
line_diatomic_or_linear_insert_query = """
    INSERT INTO lines_diatomic_or_linear
        (
            line_id, 
            J_lower, 
            symmetry_lower, 
            F_lower,
            J_upper,
            symmetry_upper,
            F_upper
        )
    VALUES
        (?,?,?,?,?,?,?)"""

line_triatomic_non_linear_insert_query = """
    INSERT INTO lines_triatomic_non_linear
        (
            line_id, 
            J_lower, 
            K_lower, 
            J_upper,
            K_upper
        )
    VALUES
        (?,?,?,?,?)"""

line_pyramidal_tetratomic_insert_query = """
    INSERT INTO lines_pyramidal_tetratomic
        (
            line_id,
            J_lower,
            K_lower,
            inversion_symmetry_lower,
            symmetry_rotational_lower,
            symmetry_total_lower,
            J_upper,
            K_upper,
            inversion_symmetry_upper,
            symmetry_rotational_upper,
            symmetry_total_upper
        )
    VALUES
        (?,?,?,?,?,?,?,?,?,?,?)"""

line_pentatomic_table_insert_query = """
    INSERT INTO lines_pentatomic
        (
            line_id,
            J_lower,
            symmetry_rotational_lower,
            alpha_lower
            F_lower,
            J_upper,
            symmetry_rotational_upper,
            alpha_upper
            F_upper,
        )
    VALUES
        (?,?,?,?,?,?,?,?,?)"""

"""
Query to insert a partition sum into the database
"""
partition_sum_insert_query = """
    INSERT INTO partition_sums
        (molecule_id, isotopologue_order_num, temperature, partition_sum)
    VALUES
        (?,?,?,?)"""

"""
Query to calculate the maximum line strength for each structure from data in the lines table
"""
structure_calculate_statistics_query = """
    UPDATE structures 
    SET max_line_strength_296K = (
        SELECT MAX(line_strength_296K) 
        FROM lines 
        WHERE lines.structure_id = structures.structure_id
    ),
    min_wavenumber_vacuum = (
        SELECT MIN(wavenumber_vacuum) 
        FROM lines 
        WHERE lines.structure_id = structures.structure_id
    ),
    max_wavenumber_vacuum = (
        SELECT MAX(wavenumber_vacuum) 
        FROM lines 
        WHERE lines.structure_id = structures.structure_id
    )"""

"""
Query to find the molecule_id and isotopologue_order_num if the molecule name and AFGL code are known. This is used
to find which isotopologue should be linked to partition_sum data, since this data is stored in a format line CO2_626.
"""
isotopologue_select_by_afgl = """
    SELECT 
        isotopologues.molecule_id, isotopologue_order_num 
    FROM isotopologues 
    INNER JOIN molecules ON molecules.molecule_id = isotopologues.molecule_id 
    WHERE molecule_name=? AND afgl_code=?"""
