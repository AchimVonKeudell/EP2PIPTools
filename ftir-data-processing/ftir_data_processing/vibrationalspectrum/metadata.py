from configparser import ConfigParser


class ConfigSaver(ConfigParser):
    def __init__(self, *args, file_path=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_path = file_path
        self.file = None

    def __enter__(self):
        if not self.file_path:
            raise ValueError("file_path must be specified for saving the configuration.")
        self.file = open(self.file_path, 'w')
        return self  # Return self to act as a ConfigParser instance

    def __exit__(self, exc_type, exc_value, traceback):
        if self.file:
            self.write(self.file)
            self.file.close()
            self.file = None


class Metadata(ConfigSaver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.add_section('Parameters')
        self.add_section('Errors')
        self.add_section('Fixed')
        self.add_section('Files')

    def set_best_fit_parameters(self, parameters: dict):
        for key, parameter in parameters.items():
            self.set("Parameters", key, str(parameter.value))
            self.set("Errors", key, str(parameter.stderr))
            self.set("Fixed", key, 'no' if parameter.vary else 'yes')

    def set_files(self, **kwargs):
        for key, item in kwargs.items():
            self.set('Files', key, item)