Install
=======

Prerequisites
~~~~~~~~~~~~~

WIFA requires Python 3.9-3.11. We recommend using `uv <https://docs.astral.sh/uv/>`_ for fast, reliable Python environment management.

**Install uv:**

.. tab-set::

    .. tab-item:: macOS/Linux

        .. code-block:: console

            curl -LsSf https://astral.sh/uv/install.sh | sh

    .. tab-item:: Windows

        .. code-block:: console

            powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

    .. tab-item:: pip

        .. code-block:: console

            pip install uv

WIFA
~~~~

**Create environment and install (recommended):**

.. code-block:: console

    # Create a new virtual environment with Python 3.11
    uv venv --python 3.11
    source .venv/bin/activate  # Linux/macOS
    # .venv\Scripts\activate   # Windows

    # Install WIFA
    uv pip install wifa

**From source (for development):**

.. code-block:: console

    git clone https://github.com/EUFLOW/WIFA.git
    cd WIFA
    uv venv --python 3.11
    source .venv/bin/activate
    uv pip install -e ".[dev,test]"

.. note::

    WIFA depends on windIO (an EU-FLOW fork) which is installed automatically.
    Each modeling tool (PyWake, foxes, wayve, code_saturne) can be installed independently.
    If you don't install one of them, that specific flow model will not be available, but other models will still work.

WindIO
~~~~~~

WindIO is installed automatically with WIFA. For manual installation:

.. code-block:: console

    uv pip install "windIO @ git+https://github.com/EUFLOW/windIO.git"

Or clone `the windIO fork <https://github.com/EUFLOW/windIO>`_ and install:

.. code-block:: console

    git clone https://github.com/EUFLOW/windIO.git
    cd windIO
    uv pip install -e .


FOXES
~~~~~

The installation of *FOXES* is described `here in the documentation <https://fraunhoferiwes.github.io/foxes/installation.html>`_.

**For the latest release:**

.. code-block:: console

    uv pip install foxes

**For the latest developments:**

.. code-block:: console

    git clone https://github.com/FraunhoferIWES/foxes.git -b dev
    cd foxes
    uv pip install -e .


PyWake
~~~~~~

The installation of *PyWake* is described in the `PyWake documentation <https://topfarm.pages.windenergy.dtu.dk/PyWake/>`_.

**For the latest release:**

.. code-block:: console

    uv pip install py_wake

**For the latest developments:**

.. code-block:: console

    git clone https://github.com/DTUWindEnergy/PyWake.git
    cd PyWake
    uv pip install -e .


WAYVE
~~~~~

WAYVE can be downloaded and installed from `GitLab <https://gitlab.kuleuven.be/TFSO-software/wayve>`_:

.. code-block:: console

    uv pip install "wayve @ git+https://gitlab.kuleuven.be/TFSO-software/wayve@dev_foxes"

Or clone and install:

.. code-block:: console

    git clone git@gitlab.kuleuven.be:TFSO-software/wayve.git
    cd wayve
    uv pip install -e .


code_saturne and salome
~~~~~~~~~~~~~~~~~~~~~~~

code_saturne and salome should be installed independently, prior to using code_saturne through the WIFA API.

**code_saturne:**

Source code and prerequisites for version 8.0 can be found at the `code_saturne download page <https://www.code-saturne.org/cms/web/Download/>`_, including the `GitHub repository <https://github.com/code-saturne/code_saturne/>`_ with a Python script for semi-automated installation.

**salome:**

Salome can be installed in two ways:

* Stand-alone `direct download <https://www.salome-platform.org/?page_id=2430/>`_
* Building the `salome_cfd extension <https://github.com/code-saturne/salome_cfd_extensions/>`_

**Configuration:**

Once installed, specify the executable paths in ``wifa/cs_api/__init__.py``:

.. code-block:: python

    # Required: add your path to code_saturne executable
    cs_exe_path = "/path/to/code_saturne"

    # Required: add your path to salome executable
    salome_exe_path = "/path/to/salome"

    # Optional: add any environment commands needed to run salome
    salome_env_command = "module load Miniforge3 && conda activate myenv"


Using pip Instead of uv
~~~~~~~~~~~~~~~~~~~~~~~

All commands in this documentation use ``uv pip`` but work identically with standard ``pip``:

.. code-block:: console

    # Using uv (recommended)
    uv pip install wifa

    # Using pip
    pip install wifa

The main advantages of uv are faster installation, better dependency resolution, and built-in virtual environment management.
