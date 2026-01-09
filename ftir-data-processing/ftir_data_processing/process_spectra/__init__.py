from . import reference_spectra, obtain_spectrum
from .collect_results import collect_results
from .obtain_spectrum import obtain_and_correct_transmittance_spectra

def actual_molecule_name(molecule_name) -> str:
    if molecule_name.find('_') > 0:
        return molecule_name[:molecule_name.find('_')]
    return molecule_name

