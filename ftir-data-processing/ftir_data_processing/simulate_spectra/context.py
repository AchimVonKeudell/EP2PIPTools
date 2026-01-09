try:
    from ..vibrationalspectrum import GetSimulator, calculate_spectrum
    from ..vibrationalspectrum.spectrumsimulator import spectrum, distributions
except ImportError:
    # Fallback for direct execution
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    from vibrationalspectrum import GetSimulator, calculate_spectrum
    from vibrationalspectrum.spectrumsimulator import spectrum, distributions