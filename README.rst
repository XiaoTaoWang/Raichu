Raichu 
======
Raichu is a cross-platform method for chromatin contact normalization.

Installation
============
Raichu and all the dependencies can be installed through either `mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_
or `pip <https://pypi.org/project/pip/>`_::

    $ conda config --append channels defaults
    $ conda config --append channels bioconda
    $ conda config --append channels conda-forge
    $ mamba create -n 3Dnorm cooler numba joblib
    $ mamba activate 3Dnorm
    $ pip install raichu
