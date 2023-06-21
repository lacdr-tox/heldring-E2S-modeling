This folder contains the ParameterOptimization.py script to do parameter estimation. It uses the odeFunctions.py script as sourcefile for the functions.

The ParmOptim script should be used with system arguments for i) the path to the data file, ii) path to the ode file and iii) path to the output folder, and optionally iv) the model name.

Example: 
python ParameterOptimization.py "/path/to/input/M16T_20220210_MH_simulation_E2data.csv" "/path/to/input/TestModel/M20T_20220304_E2model.txt" "/path/to/TestModel/TestOutput" "M20T"

# To use the script successfully, follow these steps:
1. Install anaconda
2. Create a virtual environment and activate it
3. Run the following commands
	conda install numpy pandas scipy cython seaborn
	conda install -c conda-forge pydoe
	pip install sympy==1.1.1     

# Version March 2022: option to include stress as a numeric function or as state variable in the ODE system with one or more concentrations.