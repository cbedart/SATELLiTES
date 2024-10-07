Installation of SATELLiTES (fastest install)
============================================

SATELLiTES can easily be installed `from the SATELLiTES-SBDD repository <https://pypi.org/project/SATELLiTES-SBDD/>`_ in The Python Package Index : 

.. code::

    pip install SATELLiTES-SBDD

To start the software, use the command line:

.. code::
    
    SATELLiTES


Conda environment
*****************

To avoid compatibility issues with new versions of some mandatory packages (NumPy 2.0 as an example), it is recommended to use a conda environment:
.. code::
    
    conda create -n SATELLiTES_env python=3.11
    conda activate SATELLiTES_env
    pip install SATELLiTES-SBDD

Then to launch SATELLiTES from the conda environment:
.. code::
    conda activate SATELLiTES_env
    SATELLiTES


Requirements
************

If you plan to build SATELLiTES yourself from its source code, here are the packages needed to make SATELLiTES work:

- Python3.10 or older
- RDKit 2022.03.5 or older
- Pandas
- Numpy
- PyQt5
