# Reaktoro Enabled Watertap

## 1 Introduction 
This is a repository for disseminating Reactor-enabled flowsheet units to use in various analyses, as well as supporting a simpler method for building WaterTap flowsheets using pre-configured unit models. 

For reaktoro-pse usage refer to [reaktoro-pse github](https://github.com/watertap-org/reaktoro-pse).

For reaktoro usage refer to [reaktoro](https://reaktoro.org/index.html)

## 2 Example flowsheets
### a. Reveres osmosis with acid addition and softening as pre-treatment 

This is a standard flowsheet used for demonstrating advanced water chemistry modeling with WaterTAP where reaktoro is used to model:

    a. Precipitation of solids in softening via addition of Lime and soda ash
    b. pH change via addition of acid
    c. Scaling potential in the Reverse Osmosis process.


## 3 Example analysis
The repository includes example analysis that are featured in publications that are either published, under review, or might be in pre-parathion. 
To reproduce results please review the analysis_scripts folder. 

### 3.1 Generating analysis data
Please review files in data_generation folder and run all the analysis script python files. Generally, these will use the [loop_tool](https://watertap.readthedocs.io/en/latest/how_to_guides/how_to_use_loopTool_to_explore_flowsheets.html). Review accompanying .yaml files to understand the sweeps being performed, and likely accompanying readme file. 

### 3.2 Figure generation and data processing. 
Please review the executable python files in the figure_generation folder, and run them to generate figures. The repo does not include any of the necessary data, so please refere to the readme file in the folder or generate data using the files in the data_generation folder. 

## 4 Reaktoro enabled unit models with simplified functionality. 
All unit models are designed to work with MCAS property package, and are available in unit_models folder. 

These units use the "WaterTAPFlowsheetBlock" for simple connection and management, they use standard configuration options as most IDAES and WaterTAP models, but include additional standard functions and routines. 

### 4.1 Connecting two unit models together. 
The unit models have standard ports that now feature a connect_to functionality. This will create an arc to another unit for you. 
To connect two units simply do as follows:

    m.fs.unit_a.outlet_port.connect_to(m.fs.unit_b.inlet_port)

Refer to flowsheet configuration details and tests for example usage. 


## 5 Add ma27 support for Cyipopt
We strongly recommend adding ma27 linear solver to cyipopt. To do so please follow instruction on [cyipopt website](https://cyipopt.readthedocs.io/en/stable/install.html#conda-forge-binaries-with-hsl)

## 6 Installation

#### 6.1 Contributing to the repository:

Setup requires conda or miniforge to install cyipopt and reaktoro. 

a. Install conda or miniforge. 

b. Clone/download this repository and open the command prompt inside the repository folder:

c. Run 

     conda env create -f conda_setup.yml

d. Run tests:

    pytest 

#### 6.2 Using models and tools in analysis

Install in your working environment inside conda

    conda install cyipopt reaktoro
    pip install git+https://github.com/watertap-org/reaktoro-pse.git

Import any of the flowsheet untis or flowsheet into your analysis code using

    import reaktoro_enabled_watertap as rew
    m.fs.reaktor_unit_model=rew.unit_models.reaktor_enabled_model(**kwargs)