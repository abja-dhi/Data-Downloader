from setuptools import setup, find_packages

setup(
    name='Downloader',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "mikeio",
        "pandas",
        "netCDF4",
        "requests",
    ],
)