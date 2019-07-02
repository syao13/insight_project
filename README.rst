The Insight Project 
===================

Description
-----------
My Insight project for training and testing machine learning models for Precision Chemotherapy Recommender.
The code for the web app is in the dash folder, and is also available at github.com/syao13/precisionChemoDash.git 


Dependencies
------------
Data used for this project can be downloaded from the Genomics of Drug Sensitivity in Cancer Project (link_):

.. code:: bash

   wget "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz"
   wget "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx"
   wget "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/v17.3_fitted_dose_response.xlsx"


This project depends the following Python libraries:

   * pandas_ for calculating aggregated results.
   * numpy_ and scipy_ for mathmatical calculations.
   * docopt_ for better command line interface.
   * jsonpickle_ for formatted and reusable output.
   * sklearn_ for machine learning models
   * multiprocessing_ for parallelize the model training process
   * matplotlib_ and seaborn_ for visualization

To install dependencies manually:

.. code:: bash

   pip3 install pandas
   pip3 install numpy
   pip3 install scipy
   pip3 install jsonpickle
   pip3 install sklearn
   pip3 install multiprocessing
   pip3 install matplotlib
   pip3 install seaborn


Usage
-----
'eda.ipynb' contains code for exploratory data analysis.


To train and test different model parameters:

.. code:: bash

   python3 train_models.py trained_models.txt


To build models with the proper parameters:

.. code:: bash

   python3 build_models.py built_models.txt


.. _link: https://www.cancerrxgene.org/
.. _pandas: http://pandas.pydata.org/
.. _numpy: http://www.numpy.org/
.. _scipy: https://scipy.org/scipylib/index.html
.. _jsonpickle: https://github.com/jsonpickle/jsonpickle
.. _multiprocessing_: https://docs.python.org/3.7/library/multiprocessing.html
.. _matplotlib: https://matplotlib.org/
.. _seaborn: https://seaborn.pydata.org/