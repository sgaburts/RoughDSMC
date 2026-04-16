Tools for running roguhness resolved DSMC in SPARTA

CLL_fitting
Implementation of the CLL model to fit to experimental data

surface_creation
mesh_surface.py takes a profilometry surface scan and creates an .stl file of the surface.
	Can also create random rough textures
mesh_surface_AFM.py does the same for AFM data
stl2surf.py (SPARTA tool) has to be run on the .stl to convert it to SPARTA surface format

SPARTA_configs
Contains example config files, with scattered particle output and atmosphere mixture settings files.
bathchjob.slum to run multiple jobs sequntially (uses archive bash script)
generate_configs.py automatically generates configs for parametric runs using the template.in config file.

post_processing
avg_scattering_new.py correct angular scattering profiles and drag computation
	Saves the scattering profiles and computed drag .txt files and surface angles ditribution .txt file
	Emamples Ge_Si_AFM_EXP_UCB.txt (scattering) and Ge_Si_AFM_angles.txt (angle distrbution) are provided.
Matlab angular_scattering_2.m can be used to make comparison plots of the scattering results
angular_scattering_rand_new.m is simillar but for random rough surface roughness sweeps
process.job and batch_process.job are to postprocess grid and surface files to paraview readable versions
rename_variables.py can then be used to convert variable names inside the files to meaningfull ones.
	The conversion table is in the file and needs to be update if outputed variables change

other_tools
Matlab tools to plot XPS, AFM, profilometry and experimental surface scattering resutls
Data are provided in a separate data folder