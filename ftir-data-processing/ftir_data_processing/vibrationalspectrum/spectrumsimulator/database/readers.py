import re
from dataclasses import dataclass


class HITRANType:
    """
    Specify a data type as present in HITRAN .par files, where all data is stored as ASCII. Each HITRAN type needs
    to implement the get_value method for converting a portion of the string to the requested data type
    """
    @staticmethod
    def get_value(string):
        raise NotImplementedError


class HITRANFloat(HITRANType):
    """
    Implement a data type for HITRAN .par files for float types. Sometimes, the HITRAN float types are in scientific
    notation, where the E is omitted, thus 1.67e-10 becomes 1.67-10. This HITRAN data type class can handle this
    notation.
    """
    @staticmethod
    def get_value(string):
        # fix weird HITRAN float notation where the E is sometimes not present for scientific notation
        string = re.sub("(?<=\d)-", "E-", string)
        string = re.sub("(?<=\d)\+", "E+", string)

        return float(string)


class HITRANString(HITRANType):
    """
    Implement a data type for HITRAN .par files for string types. This is essentially just returning the value that
    is found in the .par file directly without any changes in content or data type.
    """
    @staticmethod
    def get_value(string):
        return string


class HITRANInt(HITRANType):
    """
    Implement a data type for HITRAN .par files for integer types.
    """
    @staticmethod
    def get_value(string):
        return int(string)


@dataclass
class HITRANParameterMapping:
    """
    Object representing the mapping of a given part of a HITRAN line to a parameter. This is defined by the position
    where the parameter is located in the string, as well as its length of characters. Furthermore, the type of the
    parameter is defined. This can be for example: str, int or float.
    """
    position: int
    length: int
    data_type: HITRANType

    def get_value_from_string(self, string):
        selected_string = string[self.position:self.position+self.length]
        if selected_string == '#' * self.length:
            return None
        return self.data_type.get_value(selected_string)

    def get_string_from_value(self, value):
        """Reconstructing the string from HITRAN file.

        :param value:
        :return:
        """
        return f"{value:>{self.length}}"[-self.length:].upper()


@dataclass()
class HITRANParameterListMapping(HITRANParameterMapping):
    """
    Object representing the mapping of a given part of a HITRAN line to an parameter that is a list of values. This is
    defined by the position and total length of the parameter, as well as the amount of items in the list and the type
    of the list items.
    """

    item_count: int


    def get_value_from_string(self, string):
        """
        This returns a list of values present in the indicated portion of the string.
        :param string:
        :return:
        """
        field_string = string[self.position:self.position+self.length]
        characters_per_item = self.length // self.item_count

        list_item_strings = re.findall(".{%i}" % characters_per_item, field_string)
        return [self.data_type.get_value(x) for x in list_item_strings]

    def get_string_from_value(self, value):
        characters_per_item = self.length // self.item_count

        string = ''
        for i in range(self.item_count):
            string += f"{value[i]:>{characters_per_item}}"
        return string


class HITRANParameterMapper:
    """
    Object to keep a list of string position to parameter mappings. It can parse a string using the specified parameter
    mappings.
    """

    _field_mappings = {}

    def set_field_mapping(self, field_name, mapping: HITRANParameterMapping):
        """
        Set the field mapping for a field with a given name to a parameter mapping object. This overwrites any mapping
        that is already present for the given field.
        :param field_name:
        :param mapping:
        :return:
        """
        self._field_mappings[field_name] = mapping

    def get_fields_from_string(self, string):
        """
        Parses the given string into a dictionary with as keys the fieldnames, and as values the parsed parameters.
        :param string:
        :return:
        """
        fields = {}

        for field, field_mapping in self._field_mappings.items():
            try:
                fields[field] = field_mapping.get_value_from_string(string)
            except ValueError as e:
                raise ValueError(f'{e}\n{field=}, {string=}')

        return fields

    def get_string_from_fields(self, fields):

        string = ''
        for parameter_label, parameter_mapping in self._field_mappings.items():
            string += parameter_mapping.get_string_from_value(fields[parameter_label])
        return string


class HITRAN160ParameterMapper(HITRANParameterMapper):
    """
    The default mapper for the 160 character HITRAN format.
    """

    _field_mappings = {
        'molecule_id':                          HITRANParameterMapping(0, 2, HITRANInt),
        'isotopologue_order_num':               HITRANParameterMapping(2, 1, HITRANInt),
        'wavenumber_vacuum':                    HITRANParameterMapping(3, 12, HITRANFloat),
        'line_strength_296K':                   HITRANParameterMapping(15, 10, HITRANFloat),
        'einstein_A':                           HITRANParameterMapping(25, 10, HITRANFloat),
        'broadening_hw_air':                    HITRANParameterMapping(35, 5, HITRANFloat),
        'broadening_hw_self':                   HITRANParameterMapping(40, 5, HITRANFloat),
        'energy_lower_state':                   HITRANParameterMapping(45, 10, HITRANFloat),
        'broadening_temp_coefficient':          HITRANParameterMapping(55, 4, HITRANFloat),
        'pressure_shift':                       HITRANParameterMapping(59, 8, HITRANFloat),
        'global_quanta_upper':                  HITRANParameterMapping(67, 15, HITRANString),
        'global_quanta_lower':                  HITRANParameterMapping(82, 15, HITRANString),
        'local_quanta_upper':                   HITRANParameterMapping(97, 15, HITRANString),
        'local_quanta_lower':                   HITRANParameterMapping(112, 15, HITRANString),
        'uncertainty_indices':                  HITRANParameterListMapping(127, 6, HITRANInt, 6),
        'reference_indices':                    HITRANParameterListMapping(133, 12, HITRANInt, 6),
        'flag':                                 HITRANParameterMapping(145, 1, HITRANString),
        'rotational_statistical_weight_upper':  HITRANParameterMapping(146, 7, HITRANFloat),
        'rotational_statistical_weight_lower':  HITRANParameterMapping(153, 7, HITRANFloat),
        }


class HITRANCO2BroadenedParameterMapper(HITRAN160ParameterMapper):
    """
    The default mapper for the 160 character HITRAN format.
    """

    def __init__(self):
        super().__init__()
        self._co_field_mappings = self._field_mappings.copy()
        self._co_field_mappings['broadening_hw_air'] = HITRANParameterMapping(182, 7, HITRANFloat)
        self._co_field_mappings['broadening_temp_coefficient'] = HITRANParameterMapping(189, 6, HITRANFloat)
        self._co_field_mappings['pressure_shift'] = HITRANParameterMapping(195, 9, HITRANFloat)

    def get_fields_from_string(self, string):
        """
        Parses the given string into a dictionary with as keys the fieldnames, and as values the parsed parameters.
        :param string:
        :return:
        """
        fields = {}

        if len(string.rstrip()) == 204:  # extended CO format
            for field, field_mapping in self._co_field_mappings.items():
                fields[field] = field_mapping.get_value_from_string(string)
        else:
            for field, field_mapping in self._field_mappings.items():
                fields[field] = field_mapping.get_value_from_string(string)
        return fields


class HITRANHeBroadenedParameterMapper(HITRAN160ParameterMapper):

    _field_mappings = {
        'molecule_id':                          HITRANParameterMapping(0, 2, HITRANInt),
        'isotopologue_order_num':               HITRANParameterMapping(2, 1, HITRANInt),
        'wavenumber_vacuum':                    HITRANParameterMapping(3, 12, HITRANFloat),
        'line_strength_296K':                   HITRANParameterMapping(15, 10, HITRANFloat),
        'einstein_A':                           HITRANParameterMapping(25, 10, HITRANFloat),
        'broadening_hw_He':                     HITRANParameterMapping(35, 6, HITRANFloat),
        'broadening_temp_coefficient':          HITRANParameterMapping(41, 7, HITRANFloat),
        'pressure_shift':                       HITRANParameterMapping(48, 9, HITRANFloat),
    }


class HITRANLineData:
    """
    Object representing the data on a single line of a HITRAN .par file with lines of 160 characters each.
    """
    molecule_id = 0
    isotopologue_order_num = 0
    wavenumber_vacuum = 0
    line_strength_296K = 0
    einstein_A = 0
    broadening_hw_air = 0
    broadening_hw_self = 0
    energy_lower_state = 0
    broadening_temp_coefficient = 0
    pressure_shift = 0
    global_quanta_upper = ""
    global_quanta_lower = ""
    local_quanta_upper = ""
    local_quanta_lower = ""
    uncertainty_indices = 0
    reference_indices = 0
    flag = ""
    rotational_statistical_weight_upper = 0
    rotational_statistical_weight_lower = 0

    @staticmethod
    def from_line_string(line_string, mapper: HITRANParameterMapper = HITRAN160ParameterMapper()):
        """
        Retrieve a HITRAN160LineData object from a 160 character string that contains parameter values as specified by
        the HITRAN format.
        :param line_string:
        :param mapper:
        :return:
        """
        # parse fields using the specified mapper
        fields = mapper.get_fields_from_string(line_string)

        # create a new line data obect, fill it and return it
        line_data = HITRANLineData()
        for field, value in fields.items():
            line_data.__setattr__(field, value)
        return line_data


class HITRANReader:
    """
    This class can be used to read a HITRAN .par file containing lines of 160 characters each. It is used as an
    iterator and iterates over each line, returning a HITRAN160LineData object per line
    """
    def __init__(self, hitran_file_path, mapper=HITRAN160ParameterMapper()):
        self._hitran_file_path = hitran_file_path
        self._file_handle = open(hitran_file_path, 'r')
        self._line_index = 0
        self._mapper = mapper

    def __iter__(self):
        """
        Set the reader as the iterator
        :return:
        """
        return self

    def __next__(self):
        """
        Try to read the next line and validate it before returning a HITRAN160LineData object.
        :return:
        """
        self._line_index += 1

        line_string = self._file_handle.readline().rstrip()

        if line_string != "":
            return HITRANLineData.from_line_string(line_string, self._mapper)
        else:
            raise StopIteration
