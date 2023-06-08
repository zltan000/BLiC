# README #

A python package to Build Light Cones and output mock catalogues


Installation
------------

To install BLiC, you can either using pip::

    pip install blic

or using conda::

    conda install -c conda-forge blic

If you have previously installed blic, and want to upgrade to a new
released version, you should do::

    pip install blic --upgrade
    

If you would rather download and install BLiC yourself,
that is also relatively straightforward:

1. Download BLiC
^^^^^^^^^^^^^^^^^^^^

   You can download the latest tarball from::

        https://github.com/zltan000/BLiC/releases/

   Or you can clone the repository using either of the following::

        git clone https://github.com/zltan000/BLiC

   which will start out in the current stable release branch.

   Either way, cd into the BLiC directory.

2. Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^

   Some of the important dependencies are:

    - numpy
    - multiprocessing
    - pandas
    - scipy
    - h5py

   But all required dependencies should be installed automatically for you by
   pip or conda, so you should not need to worry about these.

3. Install
^^^^^^^^^^

   You can then install BliC from the local distribution.  Typically this would be the
   command::

        pip install .

   If you don't have write permission in your python distribution, you might need
   to use::

        pip install . --user


Using blic package
-----------------------

To use the whole functions to generate lightcone, you should at least do::

    >>> import BLiC
    >>> from BLiC.CalcBoxes import CalcBoxes

1. Initilize
^^^^^^^^^^^^^^^^^^^^

   To initialize BLiC, run::
    
    >>> BLiC.blic.init('/path/you/want/to/store/init/file/and/outputs')
    
   and new file (BLic.init) will be generated to get all parameters this program needed.
   .init file will be generate in a folder named BLiC_out

2. Count boxes
^^^^^^^^^^^^^^^^^^^^^^^

   Once you have finish .init file, then you can run::
   
    >>> BLiC.blic.count_boxes('/path/of/.init/file')
    
   by running this line, the fundamental preparations are done, but still some detail can be adjusted.
   
3. Adjust details
^^^^^^^^^^^^^^^^^^^^^^^

   Some built-in functions are used to adjust the inputs of lightcone, which may not contained in .init file::
   
    # set_hdf_database() function is used to set the name of input hdf5 dataset
    >>> BLiC.blic.set_hdf_database('Subhalos')
    # set_unit function() is used to change the units of positions and velocities
    # the default unit of positions is Mpc/h, and velocities in km/s
    # the input tables record position in other units, you are supposed to input a factor to change them
    # for example, the input positions are in units kpc/h, so a factor 0.001 should be added to change
    # the value, and the volocities are unchanged:
    >>> BLiC.blic.set_unit([1e-3, 1])
    # set_interpo() function is used to set the ways of porperties to do interpolation
    # the choices now only have:
    # 0: do not interpolation
    # 1: do interpolation
    >>> BLiC.blic.set_interpo(['property', 'names'], [0, 1])


   
4. Run
^^^^^^^^^^^^^^^^^^^^^^^

   All things are perpared and now run this line to submit the job::
   
    >>> BLiC.blic.run('/path/of/.init/file')
    
