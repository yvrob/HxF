# HxF: Hyper-Fidelity depletion tool for Pebble Bed Reactors

## Description
HxF is a tool that enables the user to simulate the operation of pebble bed reactors using the Cerberus interface for Serpent 2 (version >= 2.2.0). It offers the following features:

* Transport/burnup
* Pebble motion/recirculation handling through discrete motion (for FCC lattices of pebbles), or pebble position evolution through .csv files import
* Core, pebble inventory, and discarded pebbles handling and monitoring

## Tools
To run the code, modified versions of Serpent 2.2.0/2.2.1 and Cerberus/KrakenTools are required. These tools are available at Berkeley's Savio supercomputer, and the following environment variables are automatically imported: gcc/12.1.0, openmpi/4.1.4-gcc, and python/3.8.8.

The required environment to run the code is located in the ./Utils/sss_environment file. This script is called automatically by HxF.

### Serpent
The compatible version of Serpent can be found at: /global/home/groups/co_nuclear/HxF_tools/serpent2.2.0_HxF/sss2. It includes several additions such as file-based domain decomposition, extended pebble bed division capability, and modified restart file handling.

### Cerberus
The compatible version of Cerberus can be found at: "/global/home/groups/co_nuclear/HxF_tools/Cerberus_HxF". It includes the addition of compiled tables printing and reducing the amount of information printed by the code when displaying the divided materials information (e.g., sss_ov_material_fuelz<id>_burnup).

The compatible version of KrakenTools can be found at: "/global/home/groups/co_nuclear/HxF_tools/KrakenTools_HxF".

## Installation/Set up

Clone this folder and navigate to the directory.

To set up the conda environment with the required Python packages (numpy, serpentTools, mpi4py, etc.), run the following command in the terminal once:

    source ./Utils/setup_Cerberus

The setup_Cerberus script creates the conda environment in the scratch folder. By default, the current Savio paths are used for Kraken and SerpentTools, and the conda environment is called "kraken".

To load the conda environment and the necessary environment variables, run:

    source ./Utils/sss_environment

This is not required when using SLURM MPI jobs since it is included in the launch_HxF script.

## Running HxF

### Preparing the case
Before running the case, you need to prepare a proper Serpent input, a corresponding Python input script, and a job configuration file.

**Serpent Model**

The Serpent model should be a classic Pebble Bed Reactor model using the `pbed` command and dividing the pebble bed with `div [...] peb [...]` command. If you are performing depletion calculations, use the `dep` card to activate the depletion mode in Serpent. You may also need to include nuclides to be followed with the `inventory` card and detectors with the `det` card. Options such as `opti`, `memfrac`, `plot`, `pcc`, and `pop` should also be included.

Do not enable domain decomposition with `set dd`, as HxF handles the process itself. It is also unnecessary to use `set rfr` and `set rfw` for restarting the calculation.

For better organization, it is recommended that you place your models in the ./Models folder.

**Input script**

The Python input script provides to the HxF tool all the information it needs about the reactor operation. The available booleans turn on or off options, and will require additional parameters when turned on:

<ins>Transport/depletion<ins>
* `transport`: run transport calculation. Otherwise, it will just move pebbles without doing anything else.
* `resolve first`: run the first transport calculation, even though it is not needed for depletion.
* `correct`: if using predictor/corrector methods (NOT ADVISED), make sure the solver does it.
* `domain_decomposition`: use domain decomposition (required for large cases), based on a user-defined set of rules.

<ins>Motion<ins> 
* `discrete_motion`: use the discrete motion method, which consists in shuffling pebbles put in a FCC lattice up or down. Otherwise, use calculated positions evolution, from a DEM for instance.
* `looping`: if not using discrete motion, use an optimization-based algorithm which matches the last positions file with the first one to enable infinite number of steps.

<ins>Restart<ins>
* `read_first_compositions`: start with already defined composition, stored in one or several binary restart files and a .csv file containing the pebble bed information.
* `restart_calculation`: restart from a stopped calculation at a given step, stored in one or several binary restart files and a .csv file containing the pebble bed information.

<ins>Output<ins>
* `plotting`: plots results as the calculation goes (keff, core slices, discarded data).
* `saving`: save all data in .csv files at each step.
* `write_restart`: writes restart files regularly, with the number of steps defined by the user.

**Configuration file**

The configuration file is a bash script which contains the necessary information to give to SLURM and HxF. One configuration file can be used for multiple cases, but arguments in the folder must be changed accordingly:

<ins>Python scripts<ins>
* INPUT_SCRIPT: Python input
* PYTHON_SCRIPT: Main script (should not be changed)

<ins>Conda environment<ins>
* CONDA_ENVIRONMENT: Name of conda environment created with setup_Cerberus
* CONFIGURATION_FILE: Script loading conda environment and variables (should not be changed)

<ins>Cluster configuration<ins>
* GROUP_ACCOUNT: Name of the account to use (e.g. co_nuclear, fc_neutronics)
* PARTITION: Name of partition (e.g. savio, savio2, savio3, savio2_bigmem, ...)
* NNODES: Number of nodes to use
* QOS: QoS associated to group account (savio_normal for fc_neutronics, nuclear_savio_normal for co_nuclear)
* NHOURS: Time limit in hours for the job
* SUFFIX_JOB: Suffix (if needed) for the job name, otherwise takes the name of the Python input script

<ins>Serpent variables<ins>
* SERPENT_EXE: Serpent executable path
* SERPENT_DATA: Nuclear data folder path
* SERPENT_ACELIB: Name of cross sections library. Put empty string if not used/needed
* SERPENT_DECLIB: Name of decay library. Put empty string if not used/needed
* SERPENT_NFYLIB: Name of fission yields library. Put empty string if not used/needed

### Execution

Once the three components were prepared, run:

    ./launch_HXF <path to config file>

The SLURM job submission, after being summarized in the terminal, will be launched.
The overall output will be put in the .o file created in the `Logs` folder, and if errors are found, they will be printed in the .err file.

Once started and modules imported by HxF, the Serpent run is carried out in the `./Case/<output folder>/wrk_Serpent` folder. Saved .csv data is written in `./Case/<output folder>/Data` and the plots in `./Case/<output folder>/Plots`.