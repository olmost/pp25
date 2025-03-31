#include<iostream>
#include<memory>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<fstream>
#include<vector>
#include<random>
#include<chrono>
#include <stdio.h>
#include <omp.h>

using namespace std;

//NOTE : 
// 1. As a first approx, we use the central wavelengths of the J, H bands to compute the spectral slope. In this version, we calculate the slope using the magnitudes in the J, H bands.
// 2. We use a salpeter IMF between 1-100 solar masses to compute the spectra using SB99. 

// ******************************************************
// physical constants
// *****************************************************
double pi = 3.14159;
double sqrp = sqrt(pi);
double planckh = 6.62/pow(10.0,27.0); // planck constant in cgs, erg sec
double c = 3.0*pow(10.0,10.0);  // cm/sec
double G  = 6.67/pow(10.0,8.0);  // cm^3/g/sec^2
double nualpha = c*pow(10.0,8.0)/1216.0;
double sigmat = 6.6524*pow(10.0, -25.0); // cm^2
double mprot = 1.6726*pow(10.0, -24.0); // g
double mean_molw = 0.59;
double boltzk = 1.380649*pow(10.,-16.0); // cm^2 g s^-2 K^-1

double part_dens0 = 3.0; // cm^-3

double omega_m0 = 0.3075;
double omega_b0 = 0.0486;
double omega_de0 = 0.691;
double hubble_val = 0.6774;
//double vol = 100000.0;
double gyr = pow(10.0,9.0);
double Mpc = 3.086 * pow(10.0,24.0);  
double Kpc = 3.086 * pow(10.0,21.0);
double Myr = 3600.0 * 24.0 * 365.0 * pow(10.0, 6.0);
double Km = pow(10.0, 5);

double mass_sun = 2.0*pow(10.0,33.0);
double lum_sun = 3.9*pow(10.0, 33.0);// erg/s
double met_sun = 0.02;
double yr = 365.0*24.0*3600.0;

double ra = 5.0/pow(10.0,6.0); // size of graphite grain, 0.05 micro metre
double rhoa =2.25; // density of graphites, ~ 3.3gm/cm^3 for Si
double mass_dust_particle = (4.0/3.0)*pi*pow(ra,3.0)*rhoa; //gm/cm^3

// ***********************************************************
double cal_age(double z1)
{
    int i=0;
    double inte=0.0;
    double endpt=0.0;
    double sum_even=0.0;
    double sum_odd=0.0;
    double h = 0.002;
    int steps=ceil((1000.0-z1)/h);
    double *fn = new double[steps+1];
    double *z=new double[steps+1];
    
    z[0] = 1000.0;
    for(i=1;i<steps;i++)
        z[i] = z[i-1]-h;
    
    for(i=0;i<steps;i++)  
    {   
        fn[i] = 1.0/((1.0+z[i])*sqrt((omega_m0*pow((1.0+z[i]),3.0)) + omega_de0));
       // (omega_m0*pow((1.0+z[i]),3.0)) + omega_de0
    }
    
    endpt = fn[0]+fn[steps-1];
    
    for(int k=1;k<steps-1;k++)
    {
        if(k%2 == 0)
            sum_even=sum_even + fn[k];
        else 
            sum_odd=sum_odd + fn[k];
    }
    
    
    //inte = (h/45.0)*((14.0*endpt)+(24.0*sum_even)+(64.0*sum_odd));
    inte = (h/3.0)*((3.0*endpt)+(2.0*sum_even)+(4.0*sum_odd)); // Simpson's rule, pg 138, NR
    delete []fn;
    delete []z;
    return(inte);
}

// #######################################################################
double cal_dis(double z1,double z2)
{
int i=0;
double inte=0.0;
double endpt=0.0;
double sum_even=0.0;
double sum_odd=0.0;
double h = 0.0002;
int steps=ceil((z1-z2)/h);
double *fn = new double[steps+1];
double *z=new double[steps+1];

z[0] = z1;
for(i=1;i<steps;i++)
  z[i] = z[i-1]-h;

for(i=0;i<steps;i++)  
  fn[i] = 1.0/sqrt(omega_de0 + (omega_m0*pow((1.0+z[i]),3.0) ) );

endpt = fn[0]+fn[steps-1];

for(int k=1;k<steps-1;k++)
  {
   if(k%2 == 0)
      sum_even=sum_even + fn[k];
   else 
       sum_odd=sum_odd + fn[k];
  }

//inte = (h/45.0)*((14.0*endpt)+(24.0*sum_even)+(64.0*sum_odd));
inte = (h/3.0)*((3.0*endpt)+(2.0*sum_even)+(4.0*sum_odd)); // Simpson's rule, pg 138, NR
delete []fn;
delete []z;
return(inte);
}

// *********************************************************
// functions to calculate the g timescale and the g radius, as the minimum between thevirial radius, the g radius and the free-fall radius
// assuming that the total matter density profile follows the hot gas distribution
// *********************************************************


double r_cooling(double f_hot, double mg_hot, double rvir, double vvir, double tvir, double lambda_g) { 

    double coolr = pow(1.0/6.0 * mg_hot * lambda_g / mean_molw / pi / mprot / boltzk / tvir /vvir, 0.5);
    double ffr = pow(8.0/3.0 * rvir/pow(vvir, 2) * G * mg_hot / f_hot / pow(pi, 2), 0.5);

    double racc = min(coolr, ffr);
    double racc_true = min(racc, rvir);

    return(racc_true);} // Kpc

// *********************************************************
// defining omega_m(z) and deltac(z) cosmological values
// *********************************************************

double zfactor(double red) {
double zfac = pow(omega_m0 * pow((1.0 + red),3.0) + omega_de0, 1.0/3.0);
return(zfac);}

double omegamz(double red) {
double oz = omega_m0 * pow((1.0 + red), 3) / (omega_m0 * pow((1.0 + red),3.0) + omega_de0);
return(oz);}

double deltacz(double red) {
double deltaval = 18*pi*pi + 82.0*(omegamz(red)-1.0) - 39*pow(omegamz(red)-1.0,2.0);
return(deltaval);}

// ********************************************

int main()
{

int i,j,k,l=0;
double z20 =20.279;
double z19 = 18.805;
double z17 = 16.498;
double z15 = 14.763;
double z13 = 12.826;
double z12 = 11.827;
double z11 = 10.991;
double z10 = 9.961;
double z9 = 8.883;
double z8 = 8.036;
double z7 = 7.052;
double z6 = 5.982;
double z5 = 4.994;
double z4 = 4.013;
double z3 = 2.992;
double z2 = 2.006;
double z1 = 1.101;

double fc5=0.4;
double fc6=0.45;
double fc7=0.7;
double fc8=1.0;

double coeff_uv = 33.0850;
double coeff_lw = 20.87516;
double coeff_nion = 46.6392;
        
// *****************************************************
// constants
//*******************************************************

double hubble_k0 = hubble_val*pow(10.0,7.0)/Mpc;
double hubble_time = 1.0/hubble_k0/3600.0/24.0/365.0; // yr

double rho_crit = 3.0*pow(hubble_k0,2.0) / (8.0*pi*G);
double rho_0 = omega_m0*rho_crit;

double e_z12 = sqrt(omega_m0*pow((1.0+z12),3.0) + omega_de0);
double hubble12 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z12; // in cm/sec/cm

double e_z11 = sqrt(omega_m0*pow((1.0+z11),3.0) + omega_de0);
double hubble11 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z11; // in cm/sec/cm

double e_z10 = sqrt(omega_m0*pow((1.0+z10),3.0) + omega_de0);
double hubble10 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z10; // in cm/sec/cm

double e_z9 = sqrt(omega_m0*pow((1.0+z9),3.0) + omega_de0);
double hubble9 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z9; // in cm/sec/cm

double e_z8 = sqrt(omega_m0*pow((1.0+z8),3.0) + omega_de0);
double hubble8 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z8; // in cm/sec/cm

double e_z7 = sqrt(omega_m0*pow((1.0+z7),3.0) + omega_de0);
double hubble7 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z7; // in cm/sec/cm

double e_z6 = sqrt(omega_m0*pow((1.0+z6),3.0) + omega_de0);
double hubble6 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z6; // in cm/sec/cm

double e_z5 = sqrt(omega_m0*pow((1.0+z5),3.0) + omega_de0);
double hubble5 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z5; // in cm/sec/cm

double e_z4 = sqrt(omega_m0*pow((1.0+z4),3.0) + omega_de0);
double hubble4 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z4; // in cm/sec/cm

double e_z3 = sqrt(omega_m0*pow((1.0+z3),3.0) + omega_de0);
double hubble3 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z3; // in cm/sec/cm

double e_z2 = sqrt(omega_m0*pow((1.0+z2),3.0) + omega_de0);
double hubble2 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z2; // in cm/sec/cm

double e_z1 = sqrt(omega_m0*pow((1.0+z1),3.0) + omega_de0);
double hubble1 = (hubble_val*pow(10.0,7.0)/Mpc)* e_z1; // in cm/sec/cm

// #########################################################################
// physical Distance between galaxy and observer, calculated using integration
// #########################################################################
double dis_0obj12 = (c/hubble_k0)*cal_dis(z12,0.0)*(1.0+z12); // luminosity distance
double dis_0obj11 = (c/hubble_k0)*cal_dis(z11,0.0)*(1.0+z11); // luminosity distance
double dis_0obj10 = (c/hubble_k0)*cal_dis(z10,0.0)*(1.0+z10); // luminosity distance
double dis_0obj9 = (c/hubble_k0)*cal_dis(z9,0.0)*(1.0+z9); // luminosity distance
double dis_0obj8 = (c/hubble_k0)*cal_dis(z8,0.0)*(1.0+z8); // luminosity distance
double dis_0obj7 = (c/hubble_k0)*cal_dis(z7,0.0)*(1.0+z7); 
double dis_0obj6 = (c/hubble_k0)*cal_dis(z6,0.0)*(1.0+z6); 
double dis_0obj5 = (c/hubble_k0)*cal_dis(z5,0.0)*(1.0+z5);
double dis_0obj4 = (c/hubble_k0)*cal_dis(z4,0.0)*(1.0+z4);
double dis_0obj3 = (c/hubble_k0)*cal_dis(z3,0.0)*(1.0+z3);
double dis_0obj2 = (c/hubble_k0)*cal_dis(z2,0.0)*(1.0+z2);
double dis_0obj1 = (c/hubble_k0)*cal_dis(z1,0.0)*(1.0+z1);
double fac_uv12 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj12*dis_0obj12);
double fac_uv11 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj11*dis_0obj11);
double fac_uv10 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj10*dis_0obj10);
double fac_uv9 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj9*dis_0obj9);
double fac_uv8 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj8*dis_0obj8);
double fac_uv7 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj7*dis_0obj7);
double fac_uv6 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj6*dis_0obj6);
double fac_uv5 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj5*dis_0obj5);
double fac_uv4 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj4*dis_0obj4);
double fac_uv3 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj3*dis_0obj3);
double fac_uv2 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj2*dis_0obj2);
double fac_uv1 = 1500.0*1500.0/(pow(10.0,8.0)*c*4.0*pi*dis_0obj1*dis_0obj1);

double age_z0 = cal_age(0.0)/(hubble_k0*365.0*3600.0*24.0);

// ***************************************************************
// st mass fn from hmf calculator
// ****************************************************************

int nol_genmf = 801;
double *mh_genmf = new double[nol_genmf];
double *nd_genmf = new double[nol_genmf];
 
int lenth_file=0;
ifstream fin3("dndlogm_solarm_vs_mpcm3_op_z1.txt",ios::in);
    while(!fin3.eof())
     {
	      fin3>>mh_genmf[lenth_file]
           >>nd_genmf[lenth_file];
            lenth_file=lenth_file+1;
     }      
 fin3.close();
 
nol_genmf = lenth_file-1;

cout<<"nol from genmf are "<<nol_genmf<<endl;

for(i=0;i<nol_genmf;i++)
{
 mh_genmf[i] = log10(mh_genmf[i]);
 nd_genmf[i] = nd_genmf[i];
  }

lenth_file=0;
int za=286;

double *zsteps_mt = new double[286];

// ***************************************************************
// redshift steps - must be equal to merger tree steps
// ****************************************************************
 
ifstream fin2("zlist_22_to_1_step20myr.dat",ios::in);
    while(!fin2.eof())
     {
	      fin2>>zsteps_mt[lenth_file];
            lenth_file=lenth_file+1;
     }      
 fin2.close();
 
const int nosteps_mt = 286;

for(i=0;i<nosteps_mt;i++)
cout<<"z in mt "<<i<<"\t"<<zsteps_mt[i]<<endl;

// ****************************************************************
// importing values of g function as a function of temperature
// ****************************************************************

double metlist[35] = {-3, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, 2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4};

int nlines = 89;
int ncols = 36;

double coolf_matrix[90][36] {};

ifstream fin1k;
fin1k.open("cooling_functions/lambda_cooling_tot.dat",ios::in); 

for (i=0; i<nlines; i++) {
    for (j=0; j<ncols; j++) {
        fin1k >> coolf_matrix[i][j];}}

fin1k.close();

	
// ****************************************************************
// read merger tree file 
// ****************************************************************

int nomt=50000000;
cout << nomt << endl;
double *zp_mt = new double[nomt];
double *mhp_mt = new double[nomt];
int *idp_mt = new int[nomt];
int *id_parent_mt = new int[nomt];
int *flag = new int[nomt];
lenth_file=0;

ifstream fin1a("code_20to1_1halo/mt_parkinson_z20_z1_step20myr_1halo_14msun_8mres.dat",ios::in); // 
    while(!fin1a.eof())
     {
	      fin1a>>zp_mt[lenth_file] // redshift of the halo
            >>mhp_mt[lenth_file]  // mass of the halo
            >>idp_mt[lenth_file] // id of the halo
            >>id_parent_mt[lenth_file] // id of the direct descendent of the halo
            >>flag[lenth_file]; // 
            lenth_file=lenth_file+1;
            //cout<<"lenthis is "<<lenth_file<<endl;
     }      
 fin1a.close();

nomt = lenth_file-1;

cout<<"no of lines is "<<nomt<<endl;


int no_gal=0; //counting the number of merger trees - i.e. number of final halos

for(i=0;i<nomt;i++) {
   zp_mt[i] = round(zp_mt[i]*1000.0)/1000.0;
   mhp_mt[i] = log10(mhp_mt[i]);
   if(id_parent_mt[i]==-1)
      no_gal = no_gal+1;}
      
double normbin = 6. / no_gal;
 
cout<<"no of gal are "<<no_gal<<endl;

double *mh_gal = new double[no_gal];
int *id_gal = new int[no_gal];
int *plcb_gal = new int[no_gal];
int *plce_gal = new int[no_gal];
double *nd_gal = new double[no_gal];

for(i=0;i<no_gal;i++) {

    mh_gal[i]=0; id_gal[i]=0; nd_gal[i]=0; plcb_gal[i]=0; plce_gal[i]=0;}

int co=0;
for(i=0;i<nomt;i++)
   if(id_parent_mt[i]==-1)
     {
      mh_gal[co] = mhp_mt[i];
      id_gal[co] = idp_mt[i];
      co=co+1;
      }
      
int *plc_nd_gal = new int[no_gal];
for(i=0;i<no_gal;i++)
plc_nd_gal[i]=0;

double diff_min=0;

// *****************************************************************************************************
// assigning to each halo a corresponding number density according to the theoretical halo mass function (the input file has been previously imported)
// *****************************************************************************************************

for(i=0;i<no_gal;i++) {

    diff_min=0;
    diff_min = abs(mh_gal[i]-mh_genmf[0]);
  
    for(j=1;j<nol_genmf;j++) {
    
        if(abs(mh_gal[i]-mh_genmf[j])<diff_min) {
            
            plc_nd_gal[i] = plc_nd_gal[i]+1;
            diff_min = abs(mh_gal[i]-mh_genmf[j]);}}}     

for(i=0;i<no_gal;i++) {
nd_gal[i] = nd_genmf[plc_nd_gal[i]];
cout << i << "\t" << plc_nd_gal[i] << "\t" << nd_gal[i] << endl;}
  
//for(i=0;i<no_gal;i++)
//  cout<<"plc are "<<i<<"\t"<<mh_gal[i]<<"\t"<<plc_nd_gal[i]<<"\t"<<nd_gal[i]+0.46099<<endl;
     
delete []nd_genmf; delete []mh_genmf; delete []plc_nd_gal;

// ****************************************
// finding begin & end position (i.e. line) of each merger tree 
// ****************************************

co=0;

for(i=0;i<nomt;i++)
   
    if(id_parent_mt[i]==-1) {
   
     plcb_gal[co] = i;
     co=co+1;}

for(i=0;i<no_gal-1;i++)

    plce_gal[i] = plcb_gal[i+1];
 
plce_gal[no_gal-1] = nomt;

for(i=0;i<no_gal;i++) {

    int temp =plcb_gal[i];
    cout<<"plc are "<<i<<"\t"<<plcb_gal[i]<<"\t"<<plce_gal[i]<<"\t"<<zp_mt[temp]<<"\t"<<idp_mt[temp]<<endl;}

// ****************************************
// Finding the main branch of each merger tree: starting from the final halo, we walk backwards in time the merger tree selecting at each time step the most massive progenitor
// ****************************************

int *main_branch_flag = new int[nomt];

for (i = 0; i < nomt; i++)
    main_branch_flag[i] = 0;


for (i = 0; i < no_gal; i++) {

    int mainpos = plcb_gal[i];
    main_branch_flag[mainpos] = 1;

//  cout << "ciao" << "\t" << mainpos << "\t" << idp_mt[mainpos] << "\t" << i << "\t" << zp_mt[mainpos] << "\t" << plce_gal[i] << endl;

    while (flag[mainpos] != -1) {
    
        for (j = mainpos; j < plce_gal[i]; j++) {
        
            if (flag[mainpos] == idp_mt[j]) {
                main_branch_flag[j] = 1;
                mainpos = j;}}}}

// ***************************************************************
// importing the file with information about the stellar UV luminosity as a function of age, computed with STARBURST99
// ****************************************************************
int nol_spec = 500;
double *age_spec = new double[nol_spec];
double *lc_spec = new double[nol_spec];
double *llw_spec = new double[nol_spec];
double *nion_spec = new double[nol_spec];

lenth_file=0;
ifstream fin4("sal_pt1_100_lc_lw_q.txt",ios::in);
    
while(!fin4.eof()) {

    fin4>>age_spec[lenth_file] // age of the stellar population
    >>lc_spec[lenth_file] // lyman continuum luminosity
    >>llw_spec[lenth_file]
    >>nion_spec[lenth_file];
    lenth_file=lenth_file+1;}
              
fin4.close();

nol_spec = lenth_file-1;

for(i=0;i<nol_spec;i++) {

    lc_spec[i] = lc_spec[i]-6.0;
    llw_spec[i] = llw_spec[i]-6.0;
    nion_spec[i] = nion_spec[i]-6.0;}

cout<<"length spec is "<<nol_spec<<endl;
	   
// ****************************************
// defining the model parameters and variables
// ****************************************

int galnmin = 0;
int galnmax = 1;

double tstep = 20.0 * Myr; // time step of the merger tree
double sf_eff = 0.02; // maximum star formation efficiency
double fw = 0.1; // fraction of SN energy that couples to the ISM and drives winds
double vsup = 610.847*sqrt(fw); //km/s - velocity of SN winds
double mmin = 0.0; // minimum halo mass below which gas is photoionized by reionization feedback
double vmin = 0.0; //km/s - minimum halo velocity below which gas is photoionized by reionization feedback
double bh_coupled_f = 0.003; // fraction of AGN power that couples to ISM
double bh_f_av = 0.00003; // maximum fraction of cold gas mass that can be accreted by the black hole in a given time step
double alpha_res = 1.0; // parameter relative to the returning gas mass
double mhcrit0 = pow(10.0,11.25); // critical halo mass for quasar activity

double gasfraclim = 0.6; // minimum fraction of cold gas mass for quasar activity
double majmratio = 0.1; // minimum mass ratio that qualifies as a major merger
double bondif = 0.13; // Bondi accretion efficiency parameter
double jet_couplingk = 0.003; // fraction of jet energy that couples to the ISM

double fix_spin = 0.5; // fixed spin parameter
double f_spin = pow(fix_spin, 2) * pow((1.+sqrt(1.-pow(fix_spin, 2))), -2); // parameters relative to radiative efficiency calculations
double a_spin = pow(0.9663 - 0.9292*fix_spin, -0.5639);
double b_spin = pow(4.627 - 4.445*fix_spin, -0.5524);
double c_spin = pow(827.3 - 718.1*fix_spin, -0.7060);
double f_hot = 0.6; // fraction of gas accreted from the ISM in hot mode
double pz = 0.018; // fraction of gas turned into metals at each star formation burst
double rz = 0.301; // fraction of gas returned to the ISM during stellar activity

int *nprog = new int[nomt]; // number of direct progenitors of each halo
double *hmass = new double[nomt]; // total halo dark matter mass
double *mh_int = new double[nomt]; // halo dark matter mass coming from mergers (i.e. not including DM accreted below the mass resolution)
double *frac_mh = new double[nomt]; // fraction of halo DM mass coming from mergers
double *mgacc_int = new double[nomt]; // total (integrated over halo history) gas mass accreted from IGM 
double *mgacc_new = new double[nomt]; // gas mass accreted from IGM at the current time step
double *mdmacc_int = new double[nomt]; // total (integrated over halo history) DM mass accreted from intergalactic space (i.e. not fromo mergers)
double *mdmacc_new = new double[nomt]; // DM mass accreted from intergalactic space at the current time step
double *mgmerge_int = new double[nomt]; // total (integrated over halo history) gas mass accreted through mergers
double *mgmerge_new = new double[nomt]; // gas mass accreted through mergers at the current time step
double *mgout_sf_int = new double[nomt]; // total (integrated over halo history) gas mass ejected from the galaxy by SN feedback
double *mgout_agn_int = new double[nomt]; // total (integrated over halo history) gas mass ejected from the galaxy by AGN feedback
double *mgout_sf_new = new double[nomt]; // gas mass ejected from the galaxy by SN feedback at the current time step
double *mgout_agn_new = new double[nomt]; // gas mass ejected from the galaxy by AGN feedback at the current time step
double *mstarfb_int = new double[nomt]; // total stellar mass
double *mstarfb_new = new double[nomt]; // stellar mass formed at the current time step
double *mgf_int = new double[nomt]; // gas mass at the end of the time step
double *lconts_int = new double[nomt]; // UV luminosity integrated over time
double *lcontsnew_int = new double[nomt]; // instantaneous UV luminosity
double *la_int = new double[nomt]; // lyman-alpha luminosity integrated over time
double *lanew_int = new double[nomt]; // instantaneous lyman-alpha luminosity
double *vrot = new double[nomt]; // halo virial rotational velocity
double *fstarej = new double[nomt]; // star formation rate required to eject all the remaining gas in the galaxy through SN feedback
double *fstareff = new double[nomt]; // efficient star formation rate
double *uv_obsfb = new double[nomt]; // observed UV magnitude (stars)
double *uv_obsbh = new double[nomt]; // observed UV magnitude (supermassive black hole)
double *uv_obstot = new double[nomt]; // observed UV magnitude (total)
double *tvir = new double[nomt]; // halo virial temperature
double *baryon_frac = new double[nomt]; // baryonic fraction inside the halo
double *age = new double[nomt]; // age of the halo 
int *plc_age = new int[nomt]; // age divided by time step
double *rvir = new double[nomt]; // halo virial radius
double *bh_L1375 = new double[nomt]; // integrated AGN UV luminosity
double *bh_L1375_new = new double[nomt]; // instantaneous AGN UV luminosity
double *bhenergy = new double[nomt]; // AGN bolometric luminosity
double *mdm_merged_new = new double[nomt]; // DM mass accreted through mergers in the current time step
double *mdm_merged_int = new double[nomt]; // total DM mass accreted through mergers
double *mstar_merged_new = new double[nomt]; // stellar mass accreted through mergers in the current time step
double *mstar_merged_int = new double[nomt]; // total stellar mass accreted through mergers

double *bhmain = new double[nomt]; // central black hole mass
double *bhaccmass = new double[nomt]; // mass accreted by the black hole in the current time step
int *nprog_bh = new int[nomt]; // number of immediate progenitors (at the previous time step) of a given black hole
double *accrat = new double[nomt]; // Eddington ratio
double *mass_bh_merged = new double[nomt]; // total central black hole mass accreted through mergers
double *mass_bh_merged_new = new double[nomt]; // central black hole mass accreted through mergers in the current time step
double *uvflux_rat = new double[nomt]; // stellar-to-smbh UV luminosity ratio
double *bh_mag_flux = new double[nomt]; // central black hole magnetic flux

double *mgi_cold = new double[nomt]; // cold gas mass at the beginning of the time step
double *mgi_hot = new double[nomt]; // hot gas mass at the beginning of the time step
double *mgf_cold = new double[nomt]; // cold gas mass at the end of the time step
double *mgf_hot = new double[nomt]; // hot gas mass at the end of the time step
double *mgsn_cold = new double[nomt]; // cold gas mass after star formation burst 
double *mgsn_hot = new double[nomt]; // hot gas mass after star formation burst 
double *mg_res = new double[nomt]; // gas mass in the halo reservoir
double *mgcool_new = new double[nomt]; // hot gas mass cooled down in the current time step
double *mgcool_int = new double[nomt]; // total hot gas mass that has cooled down
double *bondi_acc_mass_new = new double[nomt]; // hot gas mass accreted by the black hole at Bondi rate in the current time step
double *bondi_acc_mass_int = new double[nomt]; // total hot gas mass accreted by the black hole at Bondi rate
double *mg_ret = new double[nomt]; // gas mass accreted by the galaxy from the reservoir at each time step
double *jet_power = new double[nomt]; // jet power
double *mg_heated = new double[nomt]; // cold gas mass heated up by AGN feedback

int *flag_majmerg = new int[nomt]; // 1 if the halo has had a major merger in the past, 0 otherwise
double *gasfrac_i = new double[nomt]; // fraction of cold gas in the halo at the beginning of the time step
double *time_since_mm = new double[nomt]; // time since the last major merger
double *gasfrac_last_mm = new double[nomt]; // fraction of cold gas at the time of the last major merger

int *id_main_prog = new int[nomt]; // halo id of the most massive direct progenitor 

double *radio_lum_td = new double[nomt]; // radio luminosity from the AGN accreting in quasar mode
double *radio_lum_adaf = new double[nomt]; // radio lumionsity from the AGN in jet mode
double *nuLnu_td = new double[nomt]; 
double *nuLnu_adaf = new double[nomt];

double *met_hot = new double[nomt]; // mass of metals in the hot gas
double *met_cold = new double[nomt]; // mass of metals in the cold gas
double *met_star = new double[nomt]; // mass of metals in stars
double *met_res = new double[nomt]; // mass of metals in the gas reservoir
double *z_hot_i = new double[nomt]; // metallicity of hot gas
double *z_cold_i = new double[nomt]; // metallicity of cold gas
double *z_star_i = new double[nomt]; // stellar metallicity
double *z_res_i = new double[nomt]; // metallicity of gas reservoir

for(i = 0; i < nomt; i++) { // initializing the variable arrays

    hmass[i] = pow(10.0, mhp_mt[i]); mh_int[i]=0.0; frac_mh[i]=0.0; mgacc_int[i]=0.0; mgacc_new[i]=0.0; mdmacc_int[i]=0.0; mdmacc_new[i]=0.0; mgmerge_int[i]=0.0; mgmerge_new[i]=0.0; mgout_sf_int[i]=0.0; mgout_agn_int[i]=0.0; mgout_sf_new[i]=0.0; mgout_agn_new[i]=0.0; mstarfb_int[i]=0.0; mstarfb_new[i]=0.0; mgf_int[i]=0.0; lconts_int[i]=0.0; lcontsnew_int[i]=0.0; la_int[i]=0.0; lanew_int[i]=0.0; vrot[i]=0.0; fstarej[i]=0.0; fstareff[i]=0.0; uv_obsfb[i]=0.0; uv_obsbh[i]=0.0; uv_obstot[i]=0.0; tvir[i]=0.0; baryon_frac[i]=0.0; plc_age[i]=0.0; age[i]=0.0; rvir[i]=0.0; bhaccmass[i] = 0.0; bhmain[i] = 0.0; mgsn_cold[i] = 0.0; mgsn_hot[i] = 0.0; bh_L1375[i] = 0.0; bh_L1375_new[i] = 0.0; nprog_bh[i] = 0; accrat[i] = 0.0; mass_bh_merged[i] = 0.0; mass_bh_merged_new[i] = 0.0; bhenergy[i] = 0.0; mdm_merged_new[i] = 0.0; mdm_merged_int[i] = 0.0; mstar_merged_new[i] = 0.0; mstar_merged_int[i] = 0.0; mgi_cold[i] = 0.0; mgi_hot[i] = 0.0; mgf_cold[i] = 0.0; mgf_hot[i] = 0.0; mg_res[i] = 0.0; mgcool_new[i] = 0.0; mgcool_int[i] = 0.0; bondi_acc_mass_new[i] = 0.0; bondi_acc_mass_int[i] = 0.0; mg_ret[i] = 0.0; flag_majmerg[i] = 0; gasfrac_i[i] = 0.0; time_since_mm[i] = 0.0; gasfrac_last_mm[i] = 0.0; bh_mag_flux[i] = 0.0; jet_power[i] = 0.0; mg_heated[i] = 0.0; nprog[i] = 0;  id_main_prog[i] = 0; met_hot[i] = 0.0; met_cold[i] = 0.0; met_star[i] = 0.0; met_res[i] = 0.0; z_hot_i[i] = 0.0; z_cold_i[i] = 0.0; z_star_i[i] = 0.0; z_res_i[i] = 0.0; radio_lum_td[i] = 0.0; radio_lum_adaf[i] = 0.0; nuLnu_td[i] = 0.0; nuLnu_adaf[i] = 0.0;}

double coeff_age = 2.0/(365.0*3600.0*24.0*3.0*hubble_k0);


//computing virial radius, velocity and temperature

for (i = 0; i < nomt; i++) {

    rvir[i] = 0.784 * pow((pow(10.0,mhp_mt[i])*hubble_val/pow(10.0,8.0)),(1.0/3.0)) * pow((omega_m0/omegamz(zp_mt[i])*deltacz(zp_mt[i])/(18.0*pi*pi)),(-1.0/3.0)) * 10.0/(1.0+zp_mt[i]) / hubble_val; // kpc
    vrot[i] = pow((G * pow(10.0,mhp_mt[i]) * mass_sun / rvir[i] / Kpc), 0.5) / pow(10.0,5.0); // km/s
    tvir[i] = mean_molw * mprot * pow((vrot[i]*pow(10.0,5.0)),2.0) / 2.0 / boltzk; // K
    }

for (i = 0; i < 20; i++)
cout << zp_mt[i] << "\t" << mhp_mt[i] << "\t" << rvir[i] << "\t" << vrot[i] << endl;

//**************************************************
// incorporating stellar black holes in the leaves
//**************************************************

default_random_engine generator;
double mean = -0.14;
double stddev  = 0.26;
lognormal_distribution<double> distribution(mean, stddev);

int nleaves = 0;

for (i = 0; i < nomt; i++) {

    if ((vrot[i] > vmin) && (hmass[i] > mmin))

        baryon_frac[i] = 1.0;

    if (flag[i] == -1) {
    
        mgi_cold[i] = (omega_b0 / omega_m0) * hmass[i] * baryon_frac[i];
        mh_int[i] = hmass[i];
        frac_mh[i] = 1.0;
        nleaves = nleaves + 1;}}

cout << "the number of leaves is" << "\t" << nleaves << endl;

for (i = 0; i < nomt; i++)

    if (flag[i] == -1) {

        if ((zp_mt[i] >= 13.0) && (bhmain[i] < 1000.0)  && (mhp_mt[i] >= 7.2))  { 

            bhmain[i] = 150.0; // dcbh mass 10**2 - 10**3
            bh_mag_flux[i] = (rand() % 4901 + 100)/100.0;}}

//**************************************************
// the main part of the code starts here
//**************************************************

for (int i = galnmin; i < galnmax; i++) { //t runs across the different merger trees

    for (int j = plce_gal[i] - 1; j >= plcb_gal[i]; j--) { //j runs across all the halos present in each single merger tree
        
        double mhcrit = mhcrit0 * pow(zfactor(zp_mt[j]), -3.0/8.0);
        double second_prog = 0.0;
        int pos_mainprog = 0;        
        double vrot2 = vrot[j] * vrot[j];

        if (flag[j] != -1) {
        
            for (int k = j+1; k <= plce_gal[i]; k++) { //k runs across halos at a higher redshift with respect to halo j, within the same merger tree, to look for its progenitors

                if ((id_parent_mt[k-1] == idp_mt[j]) && (id_parent_mt[k] != idp_mt[j]))
                    break;
            
                if (id_parent_mt[k] == idp_mt[j]) {

                    nprog[j] += 1;

                    if (idp_mt[k] != flag[j]) {
                    
                        mass_bh_merged_new[j] += bhmain[k];
                        mgmerge_new[j] += mgf_int[k];
                        mdm_merged_new[j] += hmass[k];
                        mstar_merged_new[j] += mstarfb_int[k];
                        mass_bh_merged[j] += bhmain[k];
                        mgmerge_int[j] += mgf_int[k];
                        mdm_merged_int[j] += hmass[k];
                        mstar_merged_int[j] += mstarfb_int[k];
                        if (mhp_mt[k] > second_prog) 
                            second_prog = mhp_mt[k];}
                        
                    else if (idp_mt[k] == flag[j]) {

                        pos_mainprog = k;
                        mass_bh_merged[j] += mass_bh_merged[k];
                        mgmerge_int[j] += mgmerge_int[k];
                        mdm_merged_int[j] += mdm_merged_int[k];
                        mstar_merged_int[j] += mstar_merged_int[k];

                        if (flag_majmerg[k] == 1)
                            flag_majmerg[j] = 1;}
                                        
//         flagging black holes based on their origin
                
                    if (bhmain[k] > 99.0)         
                        nprog_bh[j] = nprog_bh[j] + 1;            
            
                    mh_int[j] += hmass[k];                
                    bhmain[j] += bhmain[k];
                    bondi_acc_mass_int[j] += bondi_acc_mass_int[k];
                    mgi_hot[j] += mgf_hot[k];
                    mgi_cold[j] += mgf_cold[k];
                    mg_res[j] += mg_res[k];
                    mgacc_int[j] += mgacc_int[k];
                    mdmacc_int[j] += mdmacc_int[k];
                    mgcool_int[j] += mgcool_int[k];
                    mgout_sf_int[j] += mgout_sf_int[k];
                    mgout_agn_int[j] += mgout_agn_int[k];
                    mstarfb_int[j] += mstarfb_int[k];
                    met_hot[j] += met_hot[k];
                    met_cold[j] += met_cold[k];
                    met_star[j] += met_star[k];
                    met_res[j] += met_res[k];

//                  double time1 = coeff_age / sqrt(omega_de0 + (omega_m0 * pow((1.0 + zp_mt[j]), 3.0)));
//                  double time2 = coeff_age / sqrt(omega_de0 + (omega_m0 * pow((1.0 + zp_mt[k]), 3.0)));
                    age[j] = tstep/yr;
                    plc_age[j] = int(age[j]*yr/Myr/2.0 - 1.0); // we assume that all ages are 20 myr (we don't take into account contribution from older stellar populations)
                    double templum = -1.33*(log10(age[j])-6.301) + lcontsnew_int[k] + 0.462518;         
                    lconts_int[j] += pow(10.0,templum);}}

            bh_mag_flux[j] = bh_mag_flux[pos_mainprog];
            gasfrac_i[j] = mgi_cold[j]/(hmass[j]+mstarfb_int[j]+bhmain[j]+mgi_hot[j]+mgi_cold[j]);

            if (pow(10.0, second_prog)/hmass[pos_mainprog] < majmratio) {
                time_since_mm[j] = time_since_mm[pos_mainprog] + tstep;
                gasfrac_last_mm[j] = gasfrac_last_mm[pos_mainprog];}
            else {
                time_since_mm[j] = 0.0;
                gasfrac_last_mm[j] = gasfrac_i[j];
                flag_majmerg[j] = 1;}

            if (gasfrac_i[j] < gasfraclim * gasfrac_last_mm[j])
                flag_majmerg[j] = 0;}

// calculating properties for each halo after all the mergers
// dark matter mass smoothly accreted from the intergalactic space      

        frac_mh[j] = mh_int[j] / hmass[j];
        mdmacc_new[j] = (1.0 - frac_mh[j]) * hmass[j];
        mdmacc_int[j] += mdmacc_new[j];

// gas mass smoothly accreted from IGM plus the contribution of the gas mass returning to the halo from the reservoir

        z_res_i[j] = met_res[j]/mg_res[j];
        mg_ret[j] = alpha_res * mg_res[j] / (rvir[j]*Kpc) * vrot[j]*Km * tstep;
        mg_res[j] -= mg_ret[j];
                
        if ((hmass[j] > mmin) && (vrot[j] > vmin))         
            mgacc_new[j] = mdmacc_new[j] * (omega_b0 / omega_m0) + mg_ret[j];    
        else         
            mgacc_new[j] = 0.0;
                     
        mgacc_int[j] += mgacc_new[j];

// assuming that the accreted gas mass is shock-heated to the halo virial temperature
   
        if (flag[j] == -1)   
            mgi_cold[j] = mgi_cold[j];       

       else { 

            if (hmass[j] < mhcrit) {
                mgi_cold[j] += mgacc_new[j];

                if (mg_ret[j]>0.0) {

                    met_cold[j] += mg_ret[j]*met_res[j]/(mg_res[j]+mg_ret[j]);
                    met_res[j] -= mg_ret[j]*met_res[j]/(mg_res[j]+mg_ret[j]);}}

            else {

                if (zp_mt[j] > 2.) {
                    mgi_hot[j] += f_hot*mgacc_new[j];
                    mgi_cold[j] += (1.-f_hot)*mgacc_new[j];

                    if (mg_ret[j]>0.0) {

                        met_cold[j] += mg_ret[j]*(1.-f_hot)*met_res[j]/(mg_res[j]+mg_ret[j]);
                        met_hot[j] += mg_ret[j]*f_hot*met_res[j]/(mg_res[j]+mg_ret[j]);
                        met_res[j] -= mg_ret[j]*met_res[j]/(mg_res[j]+mg_ret[j]);}}

                else {                
                    mgi_hot[j] += mgacc_new[j];

                    if (mg_ret[j]>0.0) {

                        met_hot[j] += mg_ret[j]*met_res[j]/(mg_res[j]+mg_ret[j]);
                        met_res[j] -= mg_ret[j]*met_res[j]/(mg_res[j]+mg_ret[j]);}}}}

// calculating the mass fraction of the hot gas and the radius within which the hot gas is able to cool

        z_cold_i[j] = met_cold[j]/mgi_cold[j];
        z_hot_i[j] = met_hot[j]/mgi_hot[j];

        double zhot_to_zsun = 0.0;
        double r_acc = 0.0;
        int lambda_tempbin = 0;
        int lambda_metbin = 0;

        if (mgi_hot[j] > 0.0) {
            zhot_to_zsun = log10(met_hot[j]/mgi_hot[j]/met_sun);

            lambda_tempbin = int(round((log10(tvir[j]) - coolf_matrix[0][0])/0.05));
            lambda_metbin = int(round((metlist[0]-zhot_to_zsun)/(-0.1))+1);

            if (zhot_to_zsun<metlist[0]) {
                lambda_metbin = 1;
                cout << "WARNING: low Z" << "\t" << flag[j] << "\t" << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << zhot_to_zsun << endl;}
            else if (zhot_to_zsun>metlist[ncols-1]) {
                lambda_metbin = ncols;
                cout << "WARNING: high Z" << "\t" << flag[j] << "\t" << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << zhot_to_zsun << endl;}

            if (log10(tvir[j])<coolf_matrix[0][0])
                lambda_tempbin = 0;
            else if (log10(tvir[j])>coolf_matrix[nlines-1][0])
                lambda_tempbin = nlines-1;

            double lambda_val = pow(10.0, coolf_matrix[lambda_tempbin][lambda_metbin]);

            double f_ghot = mgi_hot[j]/(hmass[j]+mgi_hot[j]+mgi_cold[j]+mstarfb_int[j]+bhmain[j]);      
            r_acc = min(r_cooling(f_ghot, mgi_hot[j]*mass_sun, rvir[j]*Kpc, vrot[j]*Km, tvir[j], lambda_val), rvir[j]*Kpc)/Kpc;
            
// calculating the amount of hot gas mass that can cool down and fall towards the centre of the galaxy and the new amount of hot gas

            mgsn_hot[j] = mgi_hot[j] * exp(-0.5 * r_acc*Kpc * vrot[j]*Km / pow(rvir[j]*Kpc, 2) * tstep);
            mgcool_new[j] = mgi_hot[j] - mgsn_hot[j];
            
            if (mgcool_new[j] < 0.0) {
                mgcool_new[j] = 0.0;
                cout << "WARNING: mgcool is negative!" << endl;
                mgsn_hot[j] = mgi_hot[j];}
                
            mgcool_int[j] += mgcool_new[j];
            mgi_cold[j] += mgcool_new[j];

            met_cold[j] += mgcool_new[j]*met_hot[j]/mgi_hot[j];
            met_hot[j] -= mgcool_new[j]*met_hot[j]/mgi_hot[j];

// calculating the hot gas mass that the black hole can accrete at Bondi rate

            if (bhmain[j] > 1) {

                double rho_bh = 3.0/8.0 * mprot * mean_molw * boltzk * tvir[j] * pow(vrot[j]*Km, 3) / lambda_val / G / bhmain[j]/mass_sun;
                bondi_acc_mass_new[j] = bondif * 2.5 * pi * pow(G*bhmain[j]*mass_sun, 2) * rho_bh / pow(vrot[j]*Km, 3) * tstep / mass_sun;}
        //      bondi_acc_mass_new[j] += 15.0/16.0 * pi * G * mean_molw * mprot * boltzk * tvir[j] / lambda_val * bhmain[j] * tstep;
         
            if (bondi_acc_mass_new[j] > mgsn_hot[j]){
                bondi_acc_mass_new[j] = mgsn_hot[j];}

            bondi_acc_mass_int[j] += bondi_acc_mass_new[j];
            if (mgsn_hot[j] > 0.0)
                met_hot[j] -= bondi_acc_mass_new[j]*met_hot[j]/mgsn_hot[j];}

// computing the star formation rate 

        if (mgi_cold[j] > 0.0) {

            fstarej[j] = vrot2 / (vrot2 + (vsup * vsup));
            fstareff[j] = min(fstarej[j], sf_eff);
                   
            if ((flag[j] == -1) && (bhmain[j] > 200))        
                fstareff[j] = 0; 
            else        
                fstareff[j] = fstareff[j];
          
            mstarfb_new[j] = (mgi_cold[j] * fstareff[j])*(1.0-rz);
            
            mstarfb_int[j] += mstarfb_new[j];
            
            mgout_sf_new[j] = (mgi_cold[j] - mstarfb_new[j]) * (fstareff[j] / fstarej[j]);
            mgsn_cold[j] = (mgi_cold[j] - mstarfb_new[j]) * (1.0 - (fstareff[j] / fstarej[j]));

            if (mgsn_cold[j] < 0.0) {
                cout << "WARNING: mgsn_cold < 0" << endl;
                mgsn_cold[j] = 0.0;}

            mgout_sf_int[j] += mgout_sf_new[j];
            mg_res[j] += mgout_sf_new[j];

            met_star[j] += met_cold[j]*mstarfb_new[j]/mgi_cold[j];
            met_cold[j] += pz*mgi_cold[j]*fstareff[j] + (rz - 1.)*met_cold[j]*fstareff[j];
            met_res[j] += mgout_sf_new[j]*met_cold[j]/(mgi_cold[j]-mstarfb_new[j]);
            if (mgsn_cold[j]>0.0)
                met_cold[j] -= mgout_sf_new[j]*met_cold[j]/(mgi_cold[j]-mstarfb_new[j]);}
            else
                met_cold[j] = 0.0;

		double f_gcold = mgsn_cold[j]/(hmass[j]+mgsn_hot[j]+mgsn_cold[j]+mstarfb_int[j]+bhmain[j]);	
        double bhvaraccmass = 0.0;

        if(bhmain[j] > 1.0) {

            double eddlum = 4.0 * pi * G * bhmain[j]*mass_sun * mprot * c /sigmat;
            double eddrate = 16.0*eddlum/pow(c,2);
            double eddingtonmass = eddrate * tstep;
    
            if (flag_majmerg[j] == 1) 
                bhvaraccmass = bh_f_av * mgsn_cold[j];
            
            if (hmass[j] > mhcrit) {
                bhaccmass[j] = bondi_acc_mass_new[j] + bhvaraccmass;
                mgf_hot[j] = mgsn_hot[j]-bondi_acc_mass_new[j];
                met_cold[j] -= met_cold[j]*bhvaraccmass/mgsn_cold[j];}

            else {
                bhvaraccmass = 0.0;
                bhaccmass[j] = bondi_acc_mass_new[j];
                mgf_hot[j] = min(mgsn_hot[j]-bondi_acc_mass_new[j], 0.0);}

            accrat[j] = bhaccmass[j]*mass_sun/tstep/eddrate;

            if (accrat[j]>0.0) {

                bhenergy[j] = eddlum * a_spin * (0.985/(1.0/accrat[j]+b_spin) + 0.015/(1.0/accrat[j]+c_spin)); //erg/s
                double radeff = bhenergy[j]/bhaccmass[j]/mass_sun*tstep/pow(c,2);
                bhmain[j] += (1.0 - radeff) * bhaccmass[j];

                if ((accrat[j] < 0.01) || (accrat[j] >= 1)) {
                    jet_power[j] = 2.8 * f_spin * pow(bh_mag_flux[j]/15., 2) * bhaccmass[j]*mass_sun/tstep * pow(c, 2);
                    radio_lum_adaf[j] = 2.0*pow(10.0, 45)*bhmain[j]/pow(10.0, 9)*accrat[j]/0.01*pow(a_spin, 2.0);
                    nuLnu_adaf[j] = pow(bhmain[j]/pow(10.0, 9), 0.42)*pow(accrat[j]/0.01, 0.42)*radio_lum_adaf[j];
                }
                else {
                    radio_lum_td[j] = 2.5*pow(10.0, 43)*pow(bhmain[j]/pow(10.0, 9), 1.1)*pow(accrat[j]/0.01, 1.2)*pow(a_spin, 2.0);
                    nuLnu_td[j] = pow(bhmain[j]/pow(10.0, 9), 0.32)*pow(accrat[j]/0.01, -1.2)*radio_lum_td[j];}
                
                bh_L1375_new[j] = pow((1375.0/4450.0), 0.44) * bhenergy[j] * 4450.0 / pow(1375.0, 2) * pow(10.0, (-0.8 + 0.067 * (log10(bhenergy[j]/lum_sun) - 12.0) - 0.017 * pow((log10(bhenergy[j]/lum_sun) - 12.0), 2.0) + 0.0023 * pow((log10(bhenergy[j]/lum_sun) - 12.0), 3.0))); // erg/s/A
                bh_L1375[j] += bh_L1375_new[j];}}

        else
            mgf_hot[j] = mgsn_hot[j];

// end of bhs implementation

        if (mgsn_cold[j] < pow(10.0,-5)) {
            mgout_agn_new[j] = 0.0;
            met_cold[j] = 0.0;}
        else {
            mgout_agn_new[j] = min((mgsn_cold[j]-bhvaraccmass) * ((bh_coupled_f * bhenergy[j] + jet_couplingk * jet_power[j]) * tstep / mass_sun / ((mgsn_cold[j]-bhvaraccmass) * vrot2 * pow(Km, 2))), mgsn_cold[j]-bhvaraccmass); 
            
            met_res[j] += met_cold[j]*mgout_agn_new[j]/(mgsn_cold[j]-bhvaraccmass);
            met_cold[j] -= met_cold[j]*mgout_agn_new[j]/(mgsn_cold[j]-bhvaraccmass);     

            mgout_agn_int[j] += mgout_agn_new[j];
            mg_res[j] += mgout_agn_new[j];

            if (mgsn_cold[j] - bhvaraccmass - mgout_agn_new[j] > 0.0) {
                mg_heated[j] = min(2.0 * jet_power[j]*(0.01) / (vrot2*pow(Km, 2)) * tstep / mass_sun, mgsn_cold[j]-bhvaraccmass-mgout_agn_new[j]);
                met_hot[j] += met_cold[j]*mg_heated[j]/(mgsn_cold[j] - bhvaraccmass - mgout_agn_new[j]);
                met_cold[j] -= met_cold[j]*mg_heated[j]/(mgsn_cold[j] - bhvaraccmass - mgout_agn_new[j]);}

            mgf_cold[j] = mgsn_cold[j] - bhvaraccmass - mgout_agn_new[j] - mg_heated[j]; 
            if (mgf_cold[j] < 0.0) 
            mgf_cold[j] = 0.0;                   
            mgf_hot[j] += mg_heated[j];}

        if (mgf_cold[j] < 0.0) {
            cout << "WARNING: mgf_cold < 0!" << endl;      
            mgf_cold[j] = 0.0;}
            
        if (mgf_hot[j] == 0.0)
            met_hot[j] = 0.0;
        if (mgf_cold[j] == 0.0)
            met_cold[j] = 0.0;
        if (mg_res[j] == 0.0)
            met_res[j] = 0.0;
            
        lcontsnew_int[j] = log10(mstarfb_new[j] * pow(10.0,coeff_uv));
        lconts_int[j] += pow(10.0, lcontsnew_int[j]);      
        
        uvflux_rat[j] = bh_L1375[j] / lconts_int[j];

        if (main_branch_flag[j] == 1)
            cout << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << bhmain[j] << "\t" << mgi_cold[j] << "\t" << mgf_cold[j] << "\t" << mgi_hot[j] << "\t" << mgf_hot[j] << "\t" << mgcool_new[j] << "\t" << tvir[i] << "\t" << log10(met_star[j]/mstarfb_int[j]/met_sun) << endl;

    }}

        
double fac_hz = 1375.0*1375.0/(c*pow(10.0,8.0));
double fac_all = fac_hz/(4.0*pi*100.0*3.086*3.086*pow(10.0,18.0)*pow(10.0,18.0));
double red_lum = 2.0;

for (i = 0; i < nomt; i++) {

    if (lconts_int[i]>0) 
        uv_obsfb[i] = (-2.5*log10((lconts_int[i])*fac_all/red_lum))-48.6; //2.9

    else        
        uv_obsfb[i] = -100.0;
        
    if (bh_L1375[i] > 0.0)
        uv_obsbh[i] = (-2.5*log10((bh_L1375[i])*fac_all))-48.6; //2.9
        
    else
        uv_obsbh[i] = -100.0;

        
    if ((lconts_int[i]>0.0) || (bh_L1375[i] > 0.0))
        uv_obstot[i] = (-2.5*log10((lconts_int[i] + (bh_L1375[i])*red_lum)*fac_all/red_lum))-48.6; //2.9
        
    else
        uv_obstot[i] = -100.0;}

// printing output files 
      
const char* name_list2[13] = {"phyprop_all_z4start_nofb_fstareff_mergertree_z1.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z2.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z3.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z4.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z5.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z6.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z7.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z8.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z9.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z10.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z11.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z12.dat"};

double lenzlist2[13] = {z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12};

ofstream foutmacchianera("main_branch_allgalaxies.dat");

for (i = galnmin; i < galnmax; i++) {

    for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
    
    	if (main_branch_flag[j] == 1) {
	
        foutmacchianera	<< zp_mt[j] << "\t" << mhp_mt[j] << "\t" << 1.0-(mh_int[j]/pow(10.0, mhp_mt[j])) << "\t" << mgacc_int[j] << "\t" << mstarfb_int[j] << "\t" << mgacc_new[j] << "\t" << mstarfb_new[j] << "\t" << nprog[j] << "\t" << mgout_sf_int[j] << "\t" << mgout_agn_int[j] << "\t" << bhenergy[j] << "\t" << bhmain[j] << "\t" << nd_gal[i] << "\t" << accrat[j] << "\t" << bhaccmass[j] << "\t" << bh_L1375[j] << "\t" << mgout_sf_new[j] << "\t" << mgout_agn_new[j] << "\t" << bondi_acc_mass_int[j] << "\t" << bondi_acc_mass_new[j] << "\t" << mgi_cold[j] << "\t" << mgi_hot[j] << "\t" << mgcool_new[j] << "\t" << mgcool_int[j] << "\t" << jet_power[j] << "\t" << mg_heated[j] << "\t" << mg_res[j] << endl;}}}

foutmacchianera.close();

for (int w = 0; w < 13; w++) {

    ofstream foutpippo(name_list2[w]);

    for (i = galnmin; i < galnmax; i++) {
    
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == lenzlist2[w]) {
            
		foutpippo << << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << 1.0-(mh_int[j]/pow(10.0, mhp_mt[j])) << "\t" << mgacc_int[j] << "\t" << mstarfb_int[j] << "\t" << mgacc_new[j] << "\t" << mstarfb_new[j] << "\t" << nprog[j] << "\t" << mgout_sf_int[j] << "\t" << mgout_agn_int[j] << "\t" << bhenergy[j] << "\t" << bhmain[j] << "\t" << nd_gal[i] << "\t" << accrat[j] << "\t" << bhaccmass[j] << "\t" << bh_L1375[j] << "\t" << mgout_sf_new[j] << "\t" << mgout_agn_new[j] << "\t" << bondi_acc_mass_int[j] << "\t" << bondi_acc_mass_new[j] << "\t" << mgi_cold[j] << "\t" << mgi_hot[j] << "\t" << mgcool_new[j] << "\t" << mgcool_int[j] << "\t" << jet_power[j] << "\t" << mg_heated[j] << "\t" << mg_res[j] << endl;}}}
                 
     foutpippo.close();}

return 0;} 

