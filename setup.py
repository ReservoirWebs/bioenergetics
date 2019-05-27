from setuptools import setup, find_packages
setup(
    name="BioEnergetics",
    version="0.1",
    packages=find_packages(),

    install_requires=['numpy', 'scipy', 'matplotlib'],
    test_requires=['pytest'],
    author='Chee Sing Lee',
    license='GPLv2'
)
