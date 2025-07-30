Installation guidelines
-----------------------
To use WIFA, there are few prerequisites that need to be installed.

The prerequisites are the embeded tools, their required packages, and a windIO fork containing the modifications that were necessary to parameterize the wind farm flow. These will ultimately be merged to the main windIO package.

Except for windio, each of the modelling tools can be installed independently. If you don't install one of them, you simply need to comment the concerned "import" in the  `main_api script <https://github.com/EUFLOW/WIFA/blob/main/wifa/main_api.py>`_ after copying it for use.

WindIO
~~~~~~~~~~~~~~~~~~~~~~~
Clone `this windIO fork <https://github.com/EUFLOW/windIO>`_, and add it to your python path.


FOXES
~~~~~~~~~~~~~~~~~~~~~~~
The installation of *FOXES* is described `here in the documentation <https://fraunhoferiwes.github.io/foxes.docs/installation.html>`_.

For the latest relase, run (for `conda` consult the link above):

.. code-block:: console

  pip install foxes

For the latest developments, clone and install the *dev*
branch from `github <https://github.com/FraunhoferIWES/foxes>`_:

.. code-block:: console

  git clone git@github.com:FraunhoferIWES/foxes.git@dev
  cd foxes
  pip install -e .


code_saturne and salome
~~~~~~~~~~~~~~~~~~~~~~~
code_saturne and salome should be installed independantly, prior to using code_saturne through the FLOW api.

code_saturne, source code and prerequisites for version 8.0 can be found `using this link <https://www.code-saturne.org/cms/web/Download/>`_, namely the `github repository <https://github.com/code-saturne/code_saturne/>`_ with a python script for semi-automated installation with detailed instructions.

For salome, it can be built in two ways:

  * stand alone `direct download <https://www.salome-platform.org/?page_id=2430/>`_
  * building the `salome_cfd extension <https://github.com/code-saturne/salome_cfd_extensions/>`_


Once installed, you will need to specify to the flow api where the executables are. This can be done by modifying flow_api/cs_api folder/ __init__.py script as follows:

.. code-block:: console

  #required: add your path to cs exe
  cs_exe_path = "$YOUR_PATH_TO/code_saturne"
  #required: add your path to salome exe
  salome_exe_path = "$YOUR_PATH_TO/salome"

This script also allows to define different environment commands to allow flexibility. For example, if a conda enviroment is necessary on your cluster to be able to run salome, this can be added as:

.. code-block:: console

  #optional : add any environment that must be loaded to run salome
  salome_env_command = "module load Miniforge3 && conda activate myenv"


Pywake
~~~~~~~~~~~~~~~~~~~~~~~

WAYVE
~~~~~~~~~~~~~~~~~~~~~~~

WAYVE can be downloaded and installed from `gitlab <https://gitlab.kuleuven.be/TFSO-software/wayve>`_:

.. code-block:: console

  git clone git@gitlab.kuleuven.be:TFSO-software/wayve.git
  cd wayve
  pip install .
  