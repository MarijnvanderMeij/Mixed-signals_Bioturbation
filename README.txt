INTRODUCTION
This folder contains the model Mixed Signals, for simulating different bioturbation processes and their effects on luminescence tracers. The model contains the following files:
- a Jupyter Notebook which guides through the different functionalities of the model (Mixed-Signals.ipynb)
- a file containing all the required functions for simulations and visualizations (Mixed-Signals_functions.jl)
- a synthetic dataset of luminescence-depth information, that is used as illustration for the calibration of the model

GETTING STARTED
The model is written in the Julia language, which is a high-performance interactive scientific modelling language (https://julialang.org/). In order to run the model, Julia needs to be installed in the computer. 
Instructions for downloading Julia can be found here: https://julialang.org/downloads/. This website has downloads for different operating systems. 

The model, in the form of a Jupyter Notebook (https://jupyter.org/) can be executed in your browser via Jupyter Notebook or Jupyterlab or in the integrated development environment Visual Studio Code.
Both programs require installation of software. Please visit the following sources for installation manuals and take the indicated additional steps to include the Julia language. 

Both programs are available for different operating systems, but installation instructions might be different.

- Visual Studio Code: https://code.visualstudio.com/download
    - Follow these steps to install the Julia extension: https://code.visualstudio.com/docs/languages/julia

- Jupyter Notebook: https://jupyter.org/install
    - Open the Julia app. This opens a prompt
    - Run the following lines of code:
        using Pkg
        Pkg.add("IJulia")
    - Jupyter notebook can be installed via the link provided above, or in the Julia prompt using the following code. When installed, the same lines open Jupyter in your browser.
        using IJulia
        notebook()

RUN THE MODEL
To open the model, take the following steps:

- Visual Studio Code:
    - Open Visual Studio Code
    - File > Open Folder 
    - Navigate to the folder where the model files are stored
    - In the explorer on the left-hand side, open the file "Mixed-Signals.ipynb"

- Jupyter Notebook 
    - Open the Jupyter Notebook app. This will start the Notebook in your browser
    - Alternatively, open a notebook by running these lines of code in the Julia prompt
        using IJulia
        notebook()
    - Navigate to the folder where you stored the model files
    - Open the file "Mixed-Signals.ipynb" 

Inside the file "Mixed-Signals.ipynb" there are instructions for using the different functionalities of the model. For customization of the functions, the file "Mixed-Signals_functions.jl" can be edited.