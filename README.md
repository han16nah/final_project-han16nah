# final_project-han16nah

This is the final project for submission in the **advances geoscripting** university course.

## Goal of the assignment

- To practice applying best-practices in scientific programmingn
- To reflect on your programming habits and improve them.

## Type of analysis

I chose to do an ecological analysis using data from the GBIF database (https://www.gbif.org/). GBIF provides free and open access to biodiversity data and has an API (https://www.gbif.org/developer/summary). There is also a python client to the API, called `pygbif` (https://pypi.org/project/pygbif/), maintained by github user `sckott`.

I split the project in two parts. 

1) An **exploratory data analysis using Jupyter Notebook**, where I experimented and explored the API, some visualisations and small analysis.


2) A python script with the key analyses and functions from the notebook. Here, I tried to focus on programming style, readability and clear documentation of the code. Input can be given in configuration files (`-json`-format). I furthermore included some tests and tried to catch errors, which I stumbled upon.

## Requirements

To run the jupyter notebook and the python code, you need `python 3.8` and a number of packages:

- jupyter 
- pandas
- geopandas
- numpy
- matplotlib
- descartes
- pygbif
- mplleaflet
- pywin32
- pyproj

If using Anaconda, the environment can be set up using the `environment.yml`:

    conda env create -f environment.yml
    
or if you're working with another operating system than windows with the `environment_cross-platform.yml`:

    conda env create -f environment_cross-platform.yml
    
### Jupyter Notebook

To run the Jupyter Notebook, simply open a command prompt in the forked repository, type `Jupyter Notebook` and then click `gbif_analysis.ipynb`.

### Python code

... 
