from setuptools import setup, find_packages

setup(
    name='downloader',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "mikeio_DHI",
        "pandas",
        "netCDF4",
        "requests"
    ],
)