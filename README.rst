======
ghgpy
======

.. image:: https://img.shields.io/pypi/v/ghgpy?color=blue
   :target: https://pypi.python.org/pypi/ghgpy
   :alt: PyPI Version
.. image:: https://img.shields.io/pypi/l/swatmf
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: PyPI - License


`ghgpy` is a set of python modules for simulating greenhouse gases on crop fields.

===========================================
Greenhouse gases emission analysis
===========================================


Get data and jupyter notebooks
------------------------------

You essentially have 2 options:

Easy way
--------

- `Download the data zip file <https://github.com/spark-brc/ghgpy_tutorials/archive/refs/heads/main.zip>`_
- Unzip `ghgpy_tutorials-main.zip` to a prefered location.


Hard way (Dev mode)
-------------------

- You will need to install Git if you don't have it installed already. Downloads are available at [the link](https://git-scm.com/download). On windows, be sure to select the option that installs command-line tools  
- For Git, you will need to set up SSH keys to work with Github. To do so:
    - Go to GitHub.com and set up an account
    - On Windows, open Git Bash (on Mac/Linux, just open a terminal) and set up ssh keys if you haven't already. To do this, simply type ssh-keygen in git bash/terminal and accept all defaults (important note - when prompted for an optional passphrase, just hit return.)  
- Follow the `instructions <https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/>`_ to set up the SSH keys with your GitHub account.
- Clone the materials from GitHub.
    - Open a git bash shell from the start menu (or, on a Mac/Linux, open a terminal)
    - Navigate to the folder you made to put the course materials
    - Clone the materials by executing the following in the git bash or terminal window:


.. code-block:: bash

   git clone https://github.com/spark-brc/ghgpy_tutorials.git


============
Installation
============

To execute jupyter notebook, we need the Miniconda environment.

1. Miniconda Python:
--------------------

- If you don't already have conda installed, please download Miniconda for your operating system from https://conda.io/en/latest/miniconda.html (choose the latest version for your operating system, 64-bit). You should not need elevated rights to install this.
- Run the installer and select "only my user" when prompted. This will allow you to work with your python installation directly.

2. Set Environment and install libraries:
-----------------------------------------

- After installation, go to the START menu and select "Miniconda Prompt" to open a DOS box.
- Using the `cd <https://www.computerhope.com/issues/chusedos.htm>`_ command in the Miniconda DOS box, navigate to the location where you have `environment.yml` the file and type: 

.. code-block:: bash

   conda env create -f environment.yml

and hit ENTER.

After your virtual environment setup is complete, change the environment to `ghgpy_tut`:  

.. code-block:: bash

   conda activate ghgpy_tut

- Launch jupyter notebook 

.. code-block:: bash

   jupyter notebook


A browser window with a Jupyter notebook instance should open. Yay!


.. rubric:: Brief overview of the API

.. code-block:: python

   # use ghgpy module
   from ghgpy.runs import run_dc_multi, run_dndc
   from ghgpy.utils import PostPr
   from ghgpy.analyzer import (
      plot_oom, plot_tseries_ch4, plot_tseries_ch4_tot, plot_violin, plot_oot)



   >>> prj_dir = "project directory"
   >>> run_dc_multi(prj_dir)

   # read output and observed data
   >>> so_df = PostPr(wd).get_ch4_so_df(outfnam="ch4_multi_dc.out")
   # plot CH4 emissions in timeseries
   plot_tseries_ch4(so_df, simnam="ch4e_tot", height=3, dot=False)


