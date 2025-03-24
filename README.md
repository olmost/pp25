# pp25
Code of the model used in the Piana, Pu 2025 paper. First version of JET

The main file to compile and run is op_hp_main.cpp

In the terminal, you can type

g++-12 -o exe main.cpp  
./exe

The required input files are:

mt_parkinson_z20_z1_step20myr_9to15.dat or similar: contains the full merger trees from z=20 to z=1 for 120 halos from mass 10^9 to 10^15 solar masses, bin size is 0.05 in log mass. Resolution mass is 10^9 solar masses. Merger trees are built according to the algorithm derived in Parkinson et al. 2008.
You can create the merger tree input file by running the fortran code contained in the folder code_20to1_1halo (it will mcompute the merger tree of one single halo). To run this you can do

make trees.exe  
./trees.exe

The file will have 5 columns:

column 0 = zp_mt: redshift of the halo  
column 1 = mhp_mt: halo mass column 2 = idp_mt: halo id  
column 2 = id_parent_mt: id of the first successor (i.e. son) of the halo. Note that multiple halos might have same son  
column 3 = flag: -1 if the halo is a first leaf (i.e. has no progenitor), 0 otherwise  
column 4 = id_main_prog: id of the most massive progenitor of the halo 

sal_pt1_100_lc_lw_q.txt: contains information about the UV luminosity emitted by an instantaneous star formation event of 10^6 solar masses, assuming a Salpeter IMF between 1 and 100 solar masses and a gas metallicity Z = 0.05. Computed using the online software STARBURST99 (Leitherer+ 1999).

column 0: age of the new stellar mass (yr)  
column 1: log of the UV luminosity computed at 1375 angstrom (erg/s)  
column 2: log of the Lyman-Werner luminosity (not used here)  
column 3: log of number of ionizing photons (not used here) 

dndlogm_solarm_vs_mpcm3_op_z1.txt: the Sheth-Tormen theoretical halo mass function at z=1.

column 0: halo mass  
column 1: number density of halos of that mass (dN/d ln M)

cooling_functions/lambda_cooling_tot.dat: cooling functions as a function of temperature interpolated to obtain the results for metallicities log(Z) = [-3; 0.5; step = 0.1] (interpolated from the results found in Table 5 and onwards of Sutherland & Dopita 1993)

column 0: log temperature  
column 1-35: log cooling rate 

z_list.dat: list of redshift steps of the merger tree and hence of the simulation. Steps are separated by 20 Myr
