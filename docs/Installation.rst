************
Installation
************

conda install
#############

.. code-block::

    conda create -n xradarSat2
    conda activate -n xradarSat2
    conda install -c conda-forge xradarSat2

xsar conda package can be quite old:

To be up to date with the developpement team, it's recommended to update the installation using pip:

.. code-block::

    pip install git+https://github.com/umr-lops/radarSat2_xarray_reader.git

.. code-block::

    git clone https://github.com/umr-lops/radarSat2_xarray_reader
    cd radarSat2_xarray_reader
    # this is needed to register git filters
    git config --local include.path ../.gitconfig
    pip install -e .


.. _conda: https://docs.anaconda.com/anaconda/install/
