# Codes for identifying polyhedral units in grain boundary structures

## Python Codes

The python codes are in the `py_codes` folder. The python codes perform the following tasks:

1. Reads all the dump files one-by-one.
	+ Change the path to find the dump files.
2. Call the `analyze_gb_atoms` function in the `pum_util_funcs` module. This function reads the data in the dump file and creates a matlab data file (`.mat`).

The following steps are followed in the `analyze_gb_atoms` function:

1. Read the dump file using ovito.
2. Identify the indices of the grain boundary atoms.
3. Determine all the atoms within a certain cut-off radius of grain boundary atoms (e.g. `rCut = 3*lat_par`).
4. Create enough periodic images such that there are enough atoms sorrounding the grain boundary atoms in the unit cell (within `rCut`).
5. Save the atoms (with images), and grain boundary indices in the matlab data file (`.mat` file).


## Matlab Codes

1. `plot_poly_units`: Reads the `.mat` files created by the python codes.
	+ Calls the `cc_coors_analysis(mat_name)` function.
2. `cc_coors_analysis`: Reads the mat_file, which contains the atoms positions and indices of the grain boundary atoms.
	+ Calls the following functions:
		- `identify_gb_vv`
		- `get_polyhedra`
		- `get_polyhedra_attributes`
		- `find_sim_box_poly`


