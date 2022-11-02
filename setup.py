from setuptools import setup, find_packages
import glob

setup(
    name='radarSat2_xarray_reader',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    scripts=glob.glob('src/scripts/*.py'),
    url='https://github.com/umr-lops/radarSat2_xarray_reader',
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    include_package_data=True,
    install_requires=[
        'xmltodict',
        'glob',
        'ast',
        'numpy',
        'os',
        'xarray',
        'datatree'
    ],
    license='MIT',
    author='Reynaud Yann',
    author_email='Yann.Reynaud.2@ifremer.fr',
    description='xarray/dask distributed L1 sar file reader for radarSat2'
)
