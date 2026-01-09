"""
Setup script that reads dependencies from requirements.txt
"""
from setuptools import setup, find_packages
import os


def read_requirements(filename='requirements.txt'):
    """
    Read dependencies from requirements.txt file.
    Skips empty lines and comments.
    """
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found. Using empty requirements.")
        return []

    requirements = []
    with open(filename, 'r', encoding='utf-16') as f:
        for line in f:
            line = line.strip()

            # Skip empty lines and comments
            if line and not line.startswith('#'):
                requirements.append(line)

    return requirements


# Read requirements
install_requires = read_requirements('requirements.txt')

setup(
    name='ftir-data-processing',
    version='1.0.0',
    description='Your package description',
    author='Steijn Vervloedt',
    packages=find_packages(),
    install_requires=install_requires,
    python_requires='>=3.10'
)