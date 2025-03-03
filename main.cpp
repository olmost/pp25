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
// define the hot gas and dark matter density profiles
// *********************************************************

double corer_to_vir = 0.1;
double rcore_to_virc = 0.1;
double dm_conc = 3.5;

//double density_hot(double mg_hot, double rvir, double radius) {

//double rho_hot = 1.1725 * mg_hot/rvir/4.0/pi * 1.0/(pow(radius, 2) + pow(corer_to_vir*rvir, 2));
//return(rho_hot);}

double density_hot(double mg_hot, double rvir, double radius) {

double rho_hot = mg_hot/rvir/4.0/pi * 1.0/pow(radius, 2);
return(rho_hot);}

double density_cold(double mg_cold, double rvir, double radius) {

double rho_cold = mg_cold/4.0/pi * pow(rvir - rcore_to_virc * rvir * atan(1.0/rcore_to_virc), -1) / (pow(radius,2) + pow(rcore_to_virc*rvir,2));
return (rho_cold);}

double rho_dm(double mdm_tot, double rvir, double radius) {

double rho_dm = mdm_tot/(4.0 * pi * pow(rvir/dm_conc, 3)) / (log10(dm_conc+1.0)-(dm_conc/(dm_conc+1.0))) / (radius/rvir*dm_conc*pow((1.0+radius/rvir*dm_conc), 2));
return(rho_dm);}

// *********************************************************
// functions to calculate the cooling timescale and the cooling radius, as the minimum between the cooling radius and the free-fall radius
// assuming that the total matter density profile follows the hot gas distribution
// *********************************************************

//double r_coolt(double mg_hot, double rvir, double vvir, double tvir, double lambda_cooling) {  //cooling radius

//double coolr = pow(0.2472 / mean_molw / mprot / boltzk * mg_hot / vvir * lambda_cooling/tvir - 0.01 * pow(rvir, 2), 0.5);
//return(coolr);} // Kpc

double r_coolt(double mg_hot, double rvir, double vvir, double tvir, double lambda_cooling) { //cooling radius

double coolr = pow(1.0/6.0 * mg_hot * lambda_cooling / mean_molw / mprot / pi / boltzk / tvir / vvir, 0.5);
return(coolr);}

//double r_fft(double f_hot, double mg_hot, double rvir, double vvir) {  //free-fall radius

//double ffr = pow(1.2590 * pow(rvir/vvir, 2) * G * mg_hot / rvir / f_hot - 0.01 * pow(rvir, 2), 0.5);
//return(ffr);} // Kpc

double r_fft(double f_hot, double mg_hot, double rvir, double vvir) {  //free-fall radius

double ffr = pow(8.0/3.0 * rvir/pow(vvir, 2) * G * mg_hot / f_hot / pow(pi, 2), 0.5);
return(ffr);} // Kpc

double r_fft_cold(double f_cold, double mg_cold, double rvir, double vvir) {

double ffr = pow(8.0/3.0 * rvir/pow(vvir, 2) * G * mg_cold / f_cold / pow(pi, 2), 0.5);
return(ffr);} // Kpc

double r_cooling(double f_hot, double mg_hot, double rvir, double vvir, double tvir, double lambda_cooling) { //accretion radius

//double coolr = pow(0.2472 / mean_molw / mprot / boltzk * mg_hot / vvir * lambda_cooling/tvir - 0.01 * pow(rvir, 2), 0.5);
double coolr = pow(1.0/6.0 * mg_hot * lambda_cooling / mean_molw / pi / mprot / boltzk / tvir /vvir, 0.5);
//double ffr = pow(1.2590 * pow(rvir/vvir, 2) * G * mg_hot / rvir / f_hot - 0.01 * pow(rvir, 2), 0.5);
double ffr = pow(8.0/3.0 * rvir/pow(vvir, 2) * G * mg_hot / f_hot / pow(pi, 2), 0.5);

double racc = min(coolr, ffr);
double racc_true = min(racc, rvir);

return(racc_true);} // Kpc

//assuming that the total matter density profile follows the dark matter distribution
//double ff_radius_dmprofile(double age_h, double f_dm, double mdm_tot, double vir_r) {

//double ffr = pow(1.2590 * pow(age_h, 2) * G * mg_hot / rvir / f_hot - 0.01 * pow(rvir, 2), 0.5);
//return(ffrdm);}

// *********************************************************
// Runge-Kutta solver for cooling flows
// *********************************************************

double tstep = 20.0 * Myr;

//double cool_flow(double mg_hot, double rvir, double vvir, double rcool) {

//double mg_cool_rate = 1.1725 * mg_hot / pow(rvir, 2) * pow(rcool,3)/(pow(rcool, 2) + pow(0.1*rvir, 2)) * vvir;
//return(mg_cool_rate);}

double cool_flow(double mg_hot, double rvir, double vvir, double rcool) {

double mg_cool_rate = 0.5 * mg_hot / pow(rvir, 2) * rcool * vvir;
return(mg_cool_rate);}


double rk_routine(double mg_hot0, double rvir, double vvir, double f_hot, double tvir, double lambda) {

//double k1 = 1.1725 * mg_hot0 / pow(rvir, 2) * pow(r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda), 3)/(pow(r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda), 2) + pow(0.1*rvir, 2)) * vvir;
//double k2 = 1.1725 * (mg_hot0 + tstep/2.0*k1)  / pow(rvir, 2) * pow(r_cooling(f_hot, mg_hot0 + tstep/2.0*k1, rvir, vvir, tvir, lambda), 3)/(pow(r_cooling(f_hot, mg_hot0 + tstep/2.0*k1, rvir, vvir, tvir, lambda), 2) + pow(0.1*rvir, 2)) * vvir;
//double k3 = 1.1725 * (mg_hot0 + tstep/2.0*k2)  / pow(rvir, 2) * pow(r_cooling(f_hot, mg_hot0 + tstep/2.0*k2, rvir, vvir, tvir, lambda), 3)/(pow(r_cooling(f_hot, mg_hot0 + tstep/2.0*k2, rvir, vvir, tvir, lambda), 2) + pow(0.1*rvir, 2)) * vvir;
//double k4 = 1.1725 * (mg_hot0 + tstep*k3)  / pow(rvir, 2) * pow(r_cooling(f_hot, mg_hot0 + tstep*k3, rvir, vvir, tvir, lambda), 3)/(pow(r_cooling(f_hot, mg_hot0 + tstep*k3, rvir, vvir, tvir, lambda), 2) + pow(0.1*rvir, 2)) * vvir;

double k1 = 0.5 * mg_hot0 / pow(rvir, 2) * r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda) * vvir;
double k2 = 0.5 * (mg_hot0 + tstep/2.0*k1)  / pow(rvir, 2) * r_cooling(f_hot, mg_hot0 + tstep/2.0*k1, rvir, vvir, tvir, lambda) * vvir;
double k3 = 0.5 * (mg_hot0 + tstep/2.0*k2)  / pow(rvir, 2) * r_cooling(f_hot, mg_hot0 + tstep/2.0*k2, rvir, vvir, tvir, lambda) * vvir;
double k4 = 0.5 * (mg_hot0 + tstep*k3)  / pow(rvir, 2) * r_cooling(f_hot, mg_hot0 + tstep*k3, rvir, vvir, tvir, lambda) * vvir;

double mg_hot1 = mg_hot0 - 1.0/6.0 * tstep * (k1 + 2.0*k2 + 2.0*k3 + k4);
double mg_hotf;

if (r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda) < rvir)
	mg_hotf = mg_hot1;

else 
	mg_hotf = 0.0;
	
return(mg_hotf);}

double k1_routine(double mg_hot0, double rvir, double vvir, double f_hot, double tvir, double lambda) {

//double k1 = 1.1725 * mg_hot0 / pow(rvir, 2) * pow(r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda), 3)/(pow(r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda), 2) + pow(0.1*rvir, 2)) * vvir;
double k1 = 0.5 * mg_hot0 / pow(rvir, 2) * r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda) * vvir;
double k2;
if (r_cooling(f_hot, mg_hot0, rvir, vvir, tvir, lambda) < rvir)
	k2 = mg_hot0 - k1;
else 
	k2 = 0.0;
	
return(k2);}

// *********************************************************
// defining omega_m(z) and deltac(z)
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

double itime, ftime, exec_time;
itime = omp_get_wtime();

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
// rk solver test
// ****************************************************************

double mg0_test = 1.62 * pow(10.0, 7);
double rvir_test = 1.559;
double vrot_test = 17.05;
double tvir_test = 10387.7;
double lambda_test = 1.0;
double fhot_test = 0.128;

double rk_test = rk_routine(mg0_test*mass_sun, rvir_test*Kpc, vrot_test*Km, fhot_test, tvir_test, lambda_test)/mass_sun;
double k1_test = k1_routine(mg0_test*mass_sun, rvir_test*Kpc, vrot_test*Km, fhot_test, tvir_test, lambda_test)/mass_sun;
cout << "rk solver test:" << endl;
cout << "initial hot gas mass:" << "\t" << mg0_test << endl;
cout << "virial radius:" << "\t" << rvir_test << endl;
cout << "cooling radius:" << "\t" << r_coolt(mg0_test*mass_sun, rvir_test*Kpc, vrot_test*Km, tvir_test, lambda_test)/Kpc << endl;
cout << "free-fall radius:" << "\t" << r_fft(fhot_test, mg0_test*mass_sun, rvir_test*Kpc, vrot_test*Km)/Kpc << endl;
cout << "cooling radius:" << "\t" << r_cooling(fhot_test, mg0_test*mass_sun, rvir_test*Kpc, vrot_test*Km, tvir_test, lambda_test)/Kpc << endl;
cout << "final hot gas mass:" << "\t" << rk_test << "\t" << k1_test << endl; 

// ***************************************************************
// st mass fn from durham
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
 nd_genmf[i] = nd_genmf[i]; //fab calculates/dln m
  }

lenth_file=0;
int za=286;

double *zsteps_mt = new double[286];
 
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

double *mh_tvir4 = new double[nosteps_mt];
 
for(i=0;i<nosteps_mt;i++) {

double term1 = omegamz(zsteps_mt[i])*18.0*pi*pi/(omega_m0*deltacz(zsteps_mt[i]));
double term2 = (1/1.98)*(10.0/(1.0+zsteps_mt[i]))*pow(term1,0.3333);
mh_tvir4[i] = log10(pow(term2,1.5)*pow(10.0,8.0)/hubble_val);}

ofstream fouta("mh_tvir4_fnz.dat");
  for(i=0;i<nosteps_mt;i++)  
   fouta<<zsteps_mt[i]<<"\t"<<mh_tvir4[i]<<endl;
  fouta.close();

// ****************************************************************
// importing values of cooling function as a function of temperature
// ****************************************************************

double metlist[35] = {-3, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, 2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4};

int nlines = 89;
int ncols = 36;

double coolf_matrix[90][36] {};

ifstream fin1k;
fin1k.open("lambda_cooling_tot.dat",ios::in); 

for (i=0; i<nlines; i++) {
    for (j=0; j<ncols; j++) {
        fin1k >> coolf_matrix[i][j];}}

fin1k.close();
	
// ****************************************************************
// read fof file and length of file will give no of emm
// ****************************************************************
int nomt=50000000;
cout << nomt << endl;
double *zp_mt = new double[nomt];
double *mhp_mt = new double[nomt];
int *idp_mt = new int[nomt];
int *id_parent_mt = new int[nomt];
int *flag = new int[nomt];
lenth_file=0;

ifstream fin1a("mt_parkinson_z20_z1_step20myr_20dex_9to15.dat",ios::in); //  //output_cdm_1e7.dat
    while(!fin1a.eof())
     {
	      fin1a>>zp_mt[lenth_file] 
            >>mhp_mt[lenth_file] 
            >>idp_mt[lenth_file]
            >>id_parent_mt[lenth_file]
            >>flag[lenth_file]; 
            lenth_file=lenth_file+1;
            //cout<<"lenthis is "<<lenth_file<<endl;
     }      
 fin1a.close();

nomt = lenth_file-1;

cout<<"no of lines is "<<nomt<<endl;

int no_gal=0;

for(i=0;i<nomt;i++) {
   zp_mt[i] = round(zp_mt[i]*1000.0)/1000.0;
   mhp_mt[i] = log10(mhp_mt[i]);
   if(id_parent_mt[i]==-1)
      no_gal = no_gal+1;}
      
double normbin = 6. / no_gal;
 
cout<<"no of gal are "<<no_gal<<" "<<normbin<<endl;

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
// finding begin & end of each galaxy merger tree
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

// ***************************************************************
// lum as fn of age
// ****************************************************************
int nol_spec = 500;
double *age_spec = new double[nol_spec];
double *lc_spec = new double[nol_spec];
double *llw_spec = new double[nol_spec];
double *nion_spec = new double[nol_spec];

lenth_file=0;
ifstream fin4("sal_pt1_100_lc_lw_q.txt",ios::in);
    
while(!fin4.eof()) {

    fin4>>age_spec[lenth_file]
    >>lc_spec[lenth_file]
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

//******************************************
// reading DCBH file with z, id and number density
//******************************************

//int lengthfiledcbh = 10000;
//int nlines = 0;
//double *zdcbh = new double[lengthfiledcbh];
//double *iddcbh = new double[lengthfiledcbh];
//double *nddcbh = new double[lengthfiledcbh];

//ifstream fin8("dcbh_id_lw30.dat", ios::in);
//   while(!fin8.eof())
//     {
//	      fin8 >> zdcbh[nlines]
//           >> iddcbh[nlines]
//           >> nddcbh[nlines];
//           nlines = nlines+1;
//     }      
// fin8.close();

//lengthfiledcbh = nlines - 1;

//for (i = 0; i < lengthfiledcbh; i++) {

//    zdcbh[i] = zdcbh[i];
//    iddcbh[i] = iddcbh[i];
//    nddcbh[i] = nddcbh[i];}

// ****************************************
// Finding the main branch of each galaxy
// ****************************************

int *main_branch_flag = new int[nomt];
int *nprog = new int[nomt];

for (i = 0; i < nomt; i++)
	main_branch_flag[i] = 0;


for (i = 0; i < no_gal; i++) {

	int mainpos = plcb_gal[i];
	main_branch_flag[mainpos] = 1;

//	cout << "ciao" << "\t" << mainpos << "\t" << idp_mt[mainpos] << "\t" << i << "\t" << zp_mt[mainpos] << "\t" << plce_gal[i] << endl;

	while (flag[mainpos] != -1) {
	
		for (j = mainpos; j < plce_gal[i]; j++) {
		
			if (flag[mainpos] == idp_mt[j]) {

//				cout << i << "\t" << j << "\t" << mainpos << "\t" << zp_mt[mainpos] << endl;
		
				main_branch_flag[j] = 1;
				mainpos = j;}}}}
	   
// ****************************************
// we need halo mass, initial/final gas mass, eff SFE
// stellar mass and cont lum for each gal
// ****************************************
int galnmin = 0;
int galnmax = 121;

double radeff = 0.1;
double radeff_cold0 = 0.1;
double radeff_hot0 = 0.1;
double bhar_to_sfr = 0.005;

double sf_eff = 0.02;
double fw = 0.1;
double vsup = 610.847*sqrt(fw); //km/s
double age0 = log10(2.0*1000000.0);
double mmin = 0.0;
double vmin = 0.0; //km/s
double bh_coupled_f = 0.003;
double bh_f_av = 0.00003;
double alpha_res = 1.0;
double stunted_acc_mass = 0.0001;
double mhcrit0 = pow(10.0,11.25);
double rcold = 0.1;
//double age8 = 5.0*pow(10.0,7.0);
double a_df = 0.216;
double b_df = 1.3;
double c_df = 1.9;
double d_df = 1.0;

double gasfraclim = 0.6;
double majmratio = 0.1;
double bondif = 0.13;
double jet_couplingk = 0.003;
double jet_uvfac = 0.01;

double fix_spin = 0.5;
double f_spin = pow(fix_spin, 2) * pow((1.+sqrt(1.-pow(fix_spin, 2))), -2);
double a_spin = pow(0.9663 - 0.9292*fix_spin, -0.5639);
double b_spin = pow(4.627 - 4.445*fix_spin, -0.5524);
double c_spin = pow(827.3 - 718.1*fix_spin, -0.7060);
double eddf = 1.0;
double f_hot = 0.6;
double pz = 0.018;
double rz = 0.301;

double fc = 1.0;
double obs_uv_mag=-1.0;

double *hmass = new double[nomt];
double *mh_int = new double[nomt];
double *frac_mh = new double[nomt];
double *mgacc_int = new double[nomt];
double *mgacc_new = new double[nomt];
double *mdmacc_int = new double[nomt];
double *mdmacc_new = new double[nomt];
double *mgmerge_int = new double[nomt];
double *mgmerge_new = new double[nomt];
double *mgout_sf_int = new double[nomt];
double *mgout_agn_int = new double[nomt];
double *mgout_sf_new = new double[nomt];
double *mgout_agn_new = new double[nomt];
double *mstarfb_int = new double[nomt];
double *mstarfb_new = new double[nomt];
double *mgf_int = new double[nomt];
double *lconts_int = new double[nomt];
double *lcontsnew_int = new double[nomt];
double *la_int = new double[nomt];
double *lanew_int = new double[nomt];
double *q_int = new double[nomt];
double *qnew_int = new double[nomt];
double *llw_int = new double[nomt];
double *llwnew_int = new double[nomt];
double *vrot = new double[nomt];
double *fstarej = new double[nomt];
double *fstareff = new double[nomt];
double *uv_obsfb = new double[nomt];
double *uv_obsbh = new double[nomt];
double *uv_obstot = new double[nomt];
double *tvir = new double[nomt];
double *baryon_frac = new double[nomt];
double *age = new double[nomt];
int *plc_age = new int[nomt];
double *rvir = new double[nomt];
double *rmet = new double[nomt];
double *bh_L1375 = new double[nomt];
double *bh_L1375_new = new double[nomt];
double *bhenergy = new double[nomt];
double *mdm_merged_new = new double[nomt];
double *mdm_merged_int = new double[nomt];
double *mstar_merged_new = new double[nomt];
double *mstar_merged_int = new double[nomt];

double *uv_obsbh_ueda = new double[nomt];
double *uv_obstot_ueda = new double[nomt];
double *bh_L1375_ueda = new double[nomt];
double *bhenergy_ueda = new double[nomt];
double *uvflux_rat_ueda = new double[nomt];

double *bhmain = new double[nomt];
double *bhaccmass = new double[nomt];
int *bhmainflag = new int[nomt];
int *bhflag = new int[nomt];
int *nprog_bh = new int[nomt];
int *flag_accr = new int[nomt];
double *accrat = new double[nomt];
double *mass_bh_merged = new double[nomt];
double *mass_bh_merged_new = new double[nomt];
double *uvflux_rat = new double[nomt];
int *accr_flag = new int[nomt];
int nbhmergers_tot = 0;
double *bh_mag_flux = new double[nomt];

double *mgi_cold = new double[nomt];
double *mgi_hot = new double[nomt];
double *mgf_cold = new double[nomt];
double *mgf_hot = new double[nomt];
double *mgsn_cold = new double[nomt];
double *mgsn_hot = new double[nomt];
double *mg_res = new double[nomt];
double *mgcool_new = new double[nomt];
double *mgcool_int = new double[nomt];
double *bondi_acc_mass_new = new double[nomt];
double *bondi_acc_mass_int = new double[nomt];
double *radio_lum = new double[nomt];
double *mg_ret = new double[nomt];
double *jet_power = new double[nomt];
double *jet_fmass = new double[nomt];
double *jet_uvlum = new double[nomt];
double *mg_heated = new double[nomt];

int *flag_majmerg = new int[nomt];
double *gasfrac_i = new double[nomt];
double *time_since_mm = new double[nomt];
double *gasfrac_last_mm = new double[nomt];

int *id_main_prog = new int[nomt];
double *id_df_merge = new double[nomt];

double *radio_lum_td = new double[nomt];
double *radio_lum_adaf = new double[nomt];
double *nuLnu_td = new double[nomt];
double *nuLnu_adaf = new double[nomt];

double *met_hot = new double[nomt];
double *met_cold = new double[nomt];
double *met_star = new double[nomt];
double *met_res = new double[nomt];
double *z_hot_i = new double[nomt];
double *z_cold_i = new double[nomt];
double *z_star_i = new double[nomt];
double *z_res_i = new double[nomt];

//******************************************
// computing the age of each halo. age is reset after a major merger
//******************************************

//double *tsteps_since_form = new double[nomt];
//double *hmass_at_form = new double[nomt];
//int *mainpos = new int[nomt];
//double *age_h = new double[nomt];

//for(i = 0; i < nomt; i++) {
//	tsteps_since_form[i] = 0.0; hmass_at_form[i] = 0.0; mainpos[i] = -1; age_h[i] = 0.0;}

//for (i = galnmin; i < galnmax; i++) {
//	cout << plce_gal[i] - 1 << "\t" << plcb_gal[i] << endl;
	
//  for (j = plce_gal[i] - 1; j >= plcb_gal[i]; j--) {
    	
//     	if (flag[j] != 1) {
       		
//       		double main_m = 1.0;
//       		double smain_m = 0.1;

//   			for (k = plce_gal[i] - 1; k >= j; k--) {
//				if (id_parent_mt[k] == idp_mt[j]) {
//					if (mhp_mt[k] > main_m) {

//						smain_m = main_m;
//						main_m = mhp_mt[k]; 
//						mainpos[j] = k;}}}
			
//			if (mainpos[j] == -1) {
			
//    			cout << "main progenitor not found" << "\t" << j << "\t" << zp_mt[j] << "\t" << idp_mt[j] << "\t" << id_parent_mt[j] << "\t" << flag[j] << "\t" << id_main_prog[j] << "\t" << nprog[j] << "\t" << mhp_mt[j] << endl;
//    			flag[j] = 1;}
  			
//			if (pow(10.0, main_m)/pow(10.0, smain_m) < 5.0)												
//				tsteps_since_form[j] = 0.0;
//				cout << pow(10.0, main_m)/pow(10.0, smain_m) << endl;}
//			else 
//				tsteps_since_form[j] = tsteps_since_form[mainpos[j]] + 1.0;
			
//			age_h[j] = tsteps_since_form[j] * tstep;}}}
				
//		cout << j << "\t" << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << idp_mt[j] << "\t" << id_parent_mt[j] << "\t" << nprog[j] << "\t" << tsteps_since_form[j] << "\t" << mainpos[j] << endl;}}

//double *fstartot = new double[nomt];

for(i = 0; i < nomt; i++) {

    hmass[i] = pow(10.0, mhp_mt[i]); mh_int[i]=0.0; frac_mh[i]=0.0; mgacc_int[i]=0.0; mgacc_new[i]=0.0; mdmacc_int[i]=0.0; mdmacc_new[i]=0.0; mgmerge_int[i]=0.0; mgmerge_new[i]=0.0; mgout_sf_int[i]=0.0; mgout_agn_int[i]=0.0; mgout_sf_new[i]=0.0; mgout_agn_new[i]=0.0; mstarfb_int[i]=0.0; mstarfb_new[i]=0.0; mgf_int[i]=0.0; lconts_int[i]=0.0; lcontsnew_int[i]=0.0; la_int[i]=0.0; lanew_int[i]=0.0; q_int[i]=0.0; qnew_int[i]=0.0; llw_int[i]=0.0; llwnew_int[i]=0.0; vrot[i]=0.0; fstarej[i]=0.0; fstareff[i]=0.0; uv_obsfb[i]=0.0; uv_obsbh[i]=0.0; uv_obstot[i]=0.0; tvir[i]=0.0; baryon_frac[i]=0.0; plc_age[i]=0.0; age[i]=0.0; rvir[i]=0.0; rmet[i]=0.0; bhaccmass[i] = 0.0; bhflag[i] = 0; bhmainflag[i] = 0.0; bhmain[i] = 0.0; mgsn_cold[i] = 0.0; mgsn_hot[i] = 0.0; bh_L1375[i] = 0.0; bh_L1375_new[i] = 0.0; nprog_bh[i] = 0; flag_accr[i] = 9; accrat[i] = 0.0; mass_bh_merged[i] = 0.0; mass_bh_merged_new[i] = 0.0; bhenergy[i] = 0.0; accr_flag[i] = -1; mdm_merged_new[i] = 0.0; mdm_merged_int[i] = 0.0; mstar_merged_new[i] = 0.0; mstar_merged_int[i] = 0.0; mgi_cold[i] = 0.0; mgi_hot[i] = 0.0; mgf_cold[i] = 0.0; mgf_hot[i] = 0.0; mg_res[i] = 0.0; mgcool_new[i] = 0.0; mgcool_int[i] = 0.0; bondi_acc_mass_new[i] = 0.0; bondi_acc_mass_int[i] = 0.0; radio_lum[i] = 0.0; mg_ret[i] = 0.0; flag_majmerg[i] = 0; gasfrac_i[i] = 0.0; time_since_mm[i] = 0.0; gasfrac_last_mm[i] = 0.0; bh_mag_flux[i] = 0.0; jet_power[i] = 0.0; jet_fmass[i] = 0.0; jet_uvlum[i] = 0.0; mg_heated[i] = 0.0; nprog[i] = 0;  id_main_prog[i] = 0; id_df_merge[i] = 0.0; met_hot[i] = 0.0; met_cold[i] = 0.0; met_star[i] = 0.0; met_res[i] = 0.0; z_hot_i[i] = 0.0; z_cold_i[i] = 0.0; z_star_i[i] = 0.0; z_res_i[i] = 0.0; radio_lum_td[i] = 0.0; radio_lum_adaf[i] = 0.0; nuLnu_td[i] = 0.0; nuLnu_adaf[i] = 0.0;}

double coeff_age = 2.0/(365.0*3600.0*24.0*3.0*hubble_k0);

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

//        for (j = 0; j < lengthfiledcbh; j++) {

//            if (idp_mt[i] == iddcbh[j]) {

//                bhmain[i] = (rand() % 89999 + 10001)/10.0; // dcbh mass 10**3 - 10**4
//                bhflag[i] = 3;
//                bhmainflag[i] = 3;
//                bh_mag_flux[i] = (rand() % 4901 + 100)/100.0;}}

        if ((zp_mt[i] >= 13.0) && (bhmain[i] < 1000.0)  && (mhp_mt[i] >= 7.2))  { 

            bhmain[i] = 150.0; // dcbh mass 10**2 - 10**3
            bhflag[i] = 1;
            bhmainflag[i] = 1;
            bh_mag_flux[i] = (rand() % 4901 + 100)/100.0;}}

int num_threads = 4;
//cout << omp_get_max_threads() << endl;
omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(num_threads); // Use 4 threads for all consecutive parallel regions

int *new_ilist = new int[no_gal];

//for (int j = 0; j < no_gal; j++){
//    new_ilist[j] = floor(j/32) + 4*(j%32);
//    cout << floor(j/32) << "\t" << j%32 << "\t" << new_ilist[j] << endl;}
for (int j = 0; j < no_gal-1; j++){
    new_ilist[j] = floor(j/30) + 4*(j%30);
    cout << floor(j/30) << "\t" << j%30 << "\t" << new_ilist[j] << endl;}
new_ilist[no_gal-1] = 120;

#pragma omp parallel 
{
//cout << omp_get_num_threads() << endl;

#pragma omp for
for (int t = galnmin; t < galnmax; t++) {

    int i = new_ilist[t];
    cout<<"gal no is " << i << "\t" << omp_get_thread_num() <<endl;

// for each halo - counting the number of progenitors
    
    for (int j = plce_gal[i] - 1; j >= plcb_gal[i]; j--) {
        
        double mhcrit = mhcrit0 * pow(zfactor(zp_mt[j]), -3.0/8.0);
        double second_prog = 0.0;
        int pos_mainprog = 0;        
        double vrot2 = vrot[j] * vrot[j];

        if (flag[j] != -1) {
        
            for (int k = j+1; k <= plce_gal[i]; k++) {

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
                
                    if (bhmain[k] > 99.0) {         
                        nprog_bh[j] = nprog_bh[j] + 1;            
                        if (bhflag[j] == 0)
                            bhflag[j] = bhflag[k];                    
                        else if((bhflag[j]==1) && (bhflag[k]==1))
                            bhflag[j]=1;
                        else if ((bhflag[j]==1) && (bhflag[k]>1))
                            bhflag[j]=2;
                        else if ((bhflag[j]>1) && (bhflag[k]==1))
                            bhflag[j]=2;
                        else if ((bhflag[j]>1) && (bhflag[k]>1))
                            bhflag[j]=3;}
            
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
                    lconts_int[j] += pow(10.0,templum);
                    q_int[j] += pow(10.0, nion_spec[plc_age[j]]);
                    llw_int[j] += (pow(10.0, llw_spec[plc_age[j]]));}}

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

// assuuming that the accreted gas mass is shock-heated to the halo virial temperature
   
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
        double r_ff = 0.0;

        if (mgi_hot[j] > 0.0) {
            zhot_to_zsun = log10(met_hot[j]/mgi_hot[j]/met_sun);

            int lambda_tempbin = int(round((log10(tvir[j]) - coolf_matrix[0][0])/0.05));
            int lambda_metbin = int(round((metlist[0]-zhot_to_zsun)/(-0.1))+1);

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
            r_ff = r_fft(f_ghot, mgi_hot[j]*mass_sun, rvir[j]*Kpc, vrot[j]*Km)/Kpc;
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
                radeff = bhenergy[j]/bhaccmass[j]/mass_sun*tstep/pow(c,2);
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
        qnew_int[j] = log10(mstarfb_new[j] * pow(10.0, coeff_nion));
        q_int[j] += pow(10.0, qnew_int[j]);
        llwnew_int[j] = log10(mstarfb_new[j] * pow(10.0,coeff_lw));
        llw_int[j] += pow(10.0, llwnew_int[j]);
        
        uvflux_rat[j] = bh_L1375[j] / lconts_int[j];

        if ((i % 16 == 0) && (main_branch_flag[j] == 1)) 
            cout << zp_mt[j] << "\t" << mhp_mt[j] << "\t" << met_cold[j] << "\t" << bhmain[j] << "\t" << mgi_cold[j] << "\t" << mgf_cold[j] << "\t" << mgi_hot[j] << "\t" << mgf_hot[j] << "\t" << mgcool_new[j] << "\t" << rvir[j] << "\t" << r_ff << "\t" << r_acc << endl;

    }}
}
//		if (main_branch_flag[j] == 1)
//        	cout << zp_mt[j] << "\t" << lambda_val << "\t" << mhp_mt[j] << "\t" << mgout_sf_new[j] << "\t" << mgacc_new[j] << "\t" << mgi_hot[j] << "\t" << rho_bh << "\t" << bondi_acc_mass_new[j] << "\t" << mgsn_hot[j] << "\t" << mgcool_new[j] << "\t" << mgf_hot[j] << "\t" << mgi_cold[j] << "\t" << mgf_cold[j] << "\t" << f_ghot << "\t" << log10(bhmain[j]) << "\t" << endl;

//        cout << nprog[j] << "\t" << idp_mt[j] << "\t" << id_parent_mt[j] << "\t" << pow(10.0, mhp_mt[j]) << "\t" << rvir[j] << "\t" << vrot[j] << "\t" << tvir[j] << "\t" << lambda_val << "\t" << mgcool_new[j] << "\t" << mgi_hot[j] << "\t" << mgsn_hot[j] << "\t" << mgf_hot[j] << "\t" << mgf_cold[j] << "\t" << bondi_acc_mass[j] << "\t" << radio_lum[j]/bhenergy[j] << "\t" << f_ghot << "\t" << r_acc/rvir[j] << "\t" << bhmain[j] << endl;}}
        
double fac_hz = 1375.0*1375.0/(c*pow(10.0,8.0));
double fac_all = fac_hz/(4.0*pi*100.0*3.086*3.086*pow(10.0,18.0)*pow(10.0,18.0));
double red_lum = 2.0;

for (i = 0; i < nomt; i++) {

    if (lconts_int[i]>0) 
        uv_obsfb[i] = (-2.5*log10((lconts_int[i])*fac_all/red_lum))-48.6; //2.9

    else        
        uv_obsfb[i] = -100.0;
        
    if (bh_L1375[i] > 0.0)
        uv_obsbh[i] = (-2.5*log10((bh_L1375[i]+jet_uvlum[i])*fac_all))-48.6; //2.9
        
    else
        uv_obsbh[i] = -100.0;

        
    if ((lconts_int[i]>0.0) || (bh_L1375[i] > 0.0))
        uv_obstot[i] = (-2.5*log10((lconts_int[i] + (bh_L1375[i]+jet_uvlum[i])*red_lum)*fac_all/red_lum))-48.6; //2.9
        
    else
        uv_obstot[i] = -100.0;}

    

for (i = 0; i < nomt; i++) {

    q_int[i] = q_int[i]/1.0;
    llw_int[i] = llw_int[i]/1.0;

    double temp3 = 0.021*pow((1.0 + zp_mt[i]), 3.0) / (121.0 * 11.0);
    rmet[i] = 0.03 * pow(mstarfb_int[i], 0.2) * pow(20.0,0.4) / pow(temp3,0.2);
    //cout<<"val are "<<rmet[i]<<"\t"<<rvir[i]<<endl;
    }
    
// correcting bh luminosity and magnitude with Ueda corrections

double bh_bollum_bin[19] = {5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75};
double ueda_fac[18] = {0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.08695653723183, 0.1422303222811, 0.2055488904526, 0.2751277076365, 0.3521595865757, 0.3521595865757, 0.3521595865757, 0.3521595865757};
int shut_on_frac[18] = {(int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.08695653723183), (int)round(1.0/0.1422303222811), (int)round(1.0/0.2055488904526), (int)round(1.0/0.2751277076365), (int)round(1.0/0.3521595865757), (int)round(1.0/0.3521595865757), (int)round(1.0/0.3521595865757), (int)round(1.0/0.3521595865757)};

int counter_agn_bollum_bin[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

for(j=0;j<no_gal;j++) {

	for(i=plcb_gal[j];i<plce_gal[j];i++) {
	
		uv_obsbh_ueda[i] = uv_obsbh[i];
		uv_obstot_ueda[i] = uv_obstot[i];
		bh_L1375_ueda[i] = bh_L1375[i];
		bhenergy_ueda[i] = bhenergy[i];
		uvflux_rat_ueda[i] = uvflux_rat[i];
		
    	for (int v = 0; v < 19; v++) {
        	
    		if ((log10(bhenergy[i]/lum_sun) >= bh_bollum_bin[v]) && (log10(bhenergy[i]/lum_sun) < bh_bollum_bin[v+1])) {
    		
    			counter_agn_bollum_bin[v] = counter_agn_bollum_bin[v] + 1;
    			
				if (counter_agn_bollum_bin[v] % shut_on_frac[v] != (int)round(shut_on_frac[v]/2.)) {
				
					uv_obstot_ueda[i] = uv_obsfb[i];
					uv_obsbh_ueda[i] = 0.0;
					bh_L1375_ueda[i] = 0.0;
					bhenergy_ueda[i] = 0.0;
					uvflux_rat_ueda[i] = 0.0;}}}}}

//**********************
// SMD,SSFR ETC.
// *********************
// binning UV LF

double minuv = 0.0;
double maxuv = 30;
double stepuv = 0.5;
int nouvsteps = ceil((maxuv-minuv)/stepuv) + 1;

double *uvbin = new double[nouvsteps];
double *mid_uvbin = new double[nouvsteps-1];

uvbin[0] = minuv;

for (i = 1; i < nouvsteps; i++)
    uvbin[i] = uvbin[i-1] + stepuv;
 
for (i = 0; i < nouvsteps - 1; i++)
    mid_uvbin[i] = (uvbin[i] + uvbin[i+1]) / 2.0;
 
 // binning bh bolometric LF   

double minbhbol = 40.5;
double maxbhbol = 49.5;
double stepbhbol = 0.5;
int nobhbolsteps = ceil((maxbhbol-minbhbol)/stepbhbol) + 1;

double *bhbolbin = new double[nouvsteps];
double *mid_bhbolbin = new double[nouvsteps-1];

bhbolbin[0] = minbhbol;

for (i = 1; i < nobhbolsteps; i++)
    bhbolbin[i] = bhbolbin[i-1] + stepbhbol;
 
for (i = 0; i < nobhbolsteps - 1; i++)
    mid_bhbolbin[i] = (bhbolbin[i] + bhbolbin[i+1]) / 2.0;


int nozint = 9;
double z_int[12] = {z12, z11, z10, z9, z8, z7, z6, z5, z4, z3, z2, z1};

double msall[nosteps_mt]={};
double ms15[nosteps_mt]={};
double ms17[nosteps_mt]={};
double ms20[nosteps_mt]={};
double uvall[nosteps_mt]={};
double uv15[nosteps_mt]={};
double uv17[nosteps_mt]={};
double uv20[nosteps_mt]={};
double ms15_uncorr[nosteps_mt] = {};
double ms17_uncorr[nosteps_mt] = {};

double sfrda[nosteps_mt] = {};
double sfrd15[nosteps_mt] = {};
double sfrd17[nosteps_mt] = {};
double sfrd20[nosteps_mt] = {};
double ssfrda[nosteps_mt] = {};
double ssfrd15[nosteps_mt] = {};
double ssfrd17[nosteps_mt] = {};
double ssfrd20[nosteps_mt] = {};

for (i = 0; i < nosteps_mt; i++) {

	msall[i] = 0.0; ms15[i] = 0.0; ms17[i] = 0.0; ms20[i] = 0.0; uvall[i] = 0.0; uv15[i] = 0.0; uv17[i] = 0.0; ms15_uncorr[i] = 0.0; ms17_uncorr[i] = 0.0; sfrda[i] = 0.0; sfrd15[i] = 0.0; sfrd17[i] = 0.0; sfrd20[i] = 0.0; ssfrda[i] = 0.0; ssfrd15[i] = 0.0; ssfrd17[i] = 0.0; ssfrd20[i] = 0.0;}
	
	
for(k=0;k<nosteps_mt;k++) {

	for(j=0;j<no_gal;j++) {

		for(i=plcb_gal[j];i<plce_gal[j];i++) {

			if(zp_mt[i] == zsteps_mt[k]) {			
				msall[k] = msall[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
				uvall[k] = uvall[k] + (lconts_int[i]*fac_hz/red_lum*normbin*pow(10.0,nd_gal[j]));}
								
			if ((zp_mt[i] == zsteps_mt[k]) && (abs(uv_obsfb[i])>=15.0))
				ms15_uncorr[k] = ms15_uncorr[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
				
			if ((zp_mt[i] == zsteps_mt[k]) && (abs(uv_obsfb[i])>=17.7))
				ms17_uncorr[k] = ms17_uncorr[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
      
			if((zp_mt[i] == zsteps_mt[k]) && (abs(uv_obsfb[i])>=15.0)) {
				ms15[k] = ms15[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
				uv15[k] = uv15[k] + (lconts_int[i]*fac_hz/red_lum*normbin*pow(10.0,nd_gal[j]));}

			if((zp_mt[i] == zsteps_mt[k]) && (abs(uv_obsfb[i])>=17.7)) {
				ms17[k] = ms17[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
				uv17[k] = uv17[k] + (lconts_int[i]*fac_hz/red_lum*normbin*pow(10.0,nd_gal[j]));}

			if((zp_mt[i] == zsteps_mt[k]) && (abs(uv_obsfb[i])>=20.)) {
				ms20[k] = ms20[k] + (mstarfb_int[i]*normbin*pow(10.0,nd_gal[j]));
				uv20[k] = uv20[k] + (lconts_int[i]*fac_hz/red_lum*normbin*pow(10.0,nd_gal[j]));}

// ****************************************
// computing the star formation rate density
// ****************************************

			if (zp_mt[i] == zsteps_mt[k]) {

				sfrda[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7;
				ssfrda[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7/mstarfb_int[i];

				if (abs(uv_obsfb[i])>=15.0) {
					sfrd15[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7;
					ssfrd15[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7/mstarfb_int[i];}

				if (abs(uv_obsfb[i])>=17.7) {
					sfrd17[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7;
					ssfrd17[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7/mstarfb_int[i];}

				if (abs(uv_obsfb[i])>=20.0) {
					sfrd20[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7;
					ssfrd20[k] += mstarfb_new[i]*normbin*pow(10.0,nd_gal[j])/2.0e7/mstarfb_int[i];}}}}}

       
ofstream foutsmd("smd_op.dat");

for (k = 0; k < nosteps_mt; k++)
    foutsmd << zsteps_mt[k] << "\t" << log10(msall[k]) << "\t" << log10(ms15[k]) << "\t" << log10(ms17[k]) << "\t" << log10(ms20[k]) << "\t" << log10(ms15_uncorr[k]) << "\t" << log10(ms17_uncorr[k]) << endl;
    foutsmd.close();

ofstream foutsfrd("sfrd_all1518_z528.dat");
for (k = 0; k < nosteps_mt; k++)
    foutsfrd << zsteps_mt[k] << "\t" << log10(sfrda[k]) << "\t" << log10(sfrd15[k]) << "\t" << log10(sfrd17[k]) << "\t" << log10(sfrd20[k]) << endl;    
foutsfrd.close(); 

ofstream foutssfrd("ssfrd_all1518_z528.dat");
for (k = 0; k < nosteps_mt; k++)
    foutssfrd << zsteps_mt[k] << log10(ssfrda[k]) << "\t" << log10(ssfrd15[k]) << "\t" << log10(ssfrd17[k]) << "\t" << log10(ssfrd20[k]) << endl;    
foutssfrd.close(); 

// computing the madau star formation rate density
                          
double sfrdmadall[nosteps_mt]={};
double sfrdmad15[nosteps_mt]={};
double sfrdmad17[nosteps_mt]={};
double sfrdmad20[nosteps_mt]={};
double madau_lum = 1.15 * pow(10.0,-28);

for (i = 0; i < nosteps_mt; i++) {
    sfrdmadall[i] = 0.0; sfrdmad15[i] = 0.0; sfrdmad17[i] = 0.0; sfrdmad20[i] = 0.0;}


for(k = 0; k < nosteps_mt; k++) {
    sfrdmadall[k] = uvall[k]*madau_lum;
    sfrdmad15[k] = uv15[k]*madau_lum;
    sfrdmad17[k] = uv17[k]*madau_lum;
    sfrdmad20[k] = uv20[k]*madau_lum;}

ofstream foutsfrdm("sfrdmadau_op.dat");

for (k = 0; k < nosteps_mt; k++)
    foutsfrdm<<zsteps_mt[k]<<"\t"<<log10(sfrdmadall[k])<<"\t"<<log10(sfrdmad15[k])<<"\t"<<log10(sfrdmad17[k])<<"\t"<<log10(sfrdmad20[k])<<endl;
    foutsfrdm.close();       
 
// ****************************************    
// computing the uv luminosity density
// ****************************************

double rho_uvall[nosteps_mt] = {};
double rho_uv15[nosteps_mt] = {};
double rho_uv17[nosteps_mt] = {};
double rho_uv20[nosteps_mt] = {};
double agn_uvall_rho[nosteps_mt] = {};
double agn_uv15_rho[nosteps_mt] = {};
double agn_uv17_rho[nosteps_mt] = {};
double agn_uv20_rho[nosteps_mt] = {};

for (i = 0; i < nosteps_mt; i++) {
    rho_uvall[i] = 0.0; rho_uv15[i] = 0.0; rho_uv17[i] = 0.0; rho_uv20[i] = 0.0; agn_uvall_rho[i] = 0.0; agn_uv15_rho[i] = 0.0; agn_uv17_rho[i] = 0.0; agn_uv20_rho[i] = 0.0;}

for(k=0;k<nosteps_mt;k++) {

	for(j=0;j<no_gal;j++) {

		for(i=plcb_gal[j];i<plce_gal[j];i++) {

            if(zp_mt[i] == zsteps_mt[k])
				agn_uvall_rho[k] = agn_uvall_rho[k] + bh_L1375_ueda[i]*pow(1375.0, 2)/c/pow(10.0, 8)*normbin*pow(10.0,nd_gal[j]);

			if ((zp_mt[i] == z_int[k]) && (abs(uv_obstot_ueda[i])>=15.0))
				agn_uv15_rho[k] = agn_uv15_rho[k] + bh_L1375_ueda[i]*pow(1375.0, 2)/c/pow(10.0, 8)*normbin*pow(10.0,nd_gal[j]);

			if ((zp_mt[i] == z_int[k]) && (abs(uv_obstot_ueda[i])>=17.7))
				agn_uv17_rho[k] = agn_uv17_rho[k] + bh_L1375_ueda[i]*pow(1375.0, 2)/c/pow(10.0, 8)*normbin*pow(10.0,nd_gal[j]);				

			if ((zp_mt[i] == z_int[k]) && (abs(uv_obstot_ueda[i])>=20))
				agn_uv20_rho[k] = agn_uv20_rho[k] + bh_L1375_ueda[i]*pow(1375.0, 2)/c/pow(10.0, 8)*normbin*pow(10.0,nd_gal[j]);}}}				

ofstream foutrho("rho_uv_all1517.dat");

for (i = 0; i < nosteps_mt; i++) {
	
    rho_uvall[i] = sfrda[i]/madau_lum + agn_uvall_rho[i];
    rho_uv15[i] = sfrd15[i]/madau_lum + agn_uv15_rho[i];
    rho_uv17[i] = sfrd17[i]/madau_lum + agn_uv17_rho[i];
    rho_uv20[i] = sfrd20[i]/madau_lum + agn_uv20_rho[i];
    
    foutrho << zsteps_mt[i] << "\t"  << log10(rho_uvall[i]) << "\t" << log10(rho_uv15[i]) << "\t" << log10(rho_uv17[i]) << "\t" << sfrda[i]/madau_lum << "\t" << sfrd15[i]/madau_lum << "\t" << sfrd17[i]/madau_lum << "\t" << agn_uvall_rho[i] << "\t" << agn_uv15_rho[i] << "\t" << agn_uv17_rho[i] << "\t" << log10(rho_uv20[i]) << "\t" << sfrd20[i]/madau_lum << "\t" << agn_uv20_rho[i] << "\t" << uvall[i] << "\t" << uv15[i] << "\t" << uv17[i] << "\t" << uv20[i] << endl;}
    
foutrho.close();

// ****************************************
// UV LF from z = 5 to z = 20
// ****************************************

double z_list[12] = {5.0, 6.06, 7.0, 7.94, 8.97, 10.04, 11.05, 11.88, 13.42, 15.52, 17.42, 18.6};
const char* name_list_AGN[12] = {"AGN_UVLF_z5.dat", "AGN_UVLF_z6.dat", "AGN_UVLF_z7.dat", "AGN_UVLF_z8.dat", "AGN_UVLF_z9.dat", "AGN_UVLF_z10.dat", "AGN_UVLF_z11.dat", "AGN_UVLF_z12.dat", "AGN_UVLF_z13.dat", "AGN_UVLF_z15.dat", "AGN_UVLF_z17.dat", "AGN_UVLF_z19.dat"};
const char* name_list_QSO[12] = {"QSO_UVLF_z5.dat", "QSO_UVLF_z6.dat", "QSO_UVLF_z7.dat", "QSO_UVLF_z8.dat", "QSO_UVLF_z9.dat", "QSO_UVLF_z10.dat", "QSO_UVLF_z11.dat", "QSO_UVLF_z12.dat", "QSO_UVLF_z13.dat", "QSO_UVLF_z15.dat", "QSO_UVLF_z17.dat", "QSO_UVLF_z19.dat"};
const char* name_list_QSO_bhbol[12] = {"QSO_UVLF_bhbol_z5.dat", "QSO_UVLF_bhbol_z6.dat", "QSO_UVLF_bhbol_z7.dat", "QSO_UVLF_bhbol_z8.dat", "QSO_UVLF_bhbol_z9.dat", "QSO_UVLF_bhbol_z10.dat", "QSO_UVLF_bhbol_z11.dat", "QSO_UVLF_bhbol_z12.dat", "QSO_UVLF_bhbol_z13.dat", "QSO_UVLF_bhbol_z15.dat", "QSO_UVLF_bhbol_z17.dat", "QSO_UVLF_bhbol_z19.dat"};
const char* name_list_QSO_bhbol_nonstoch[12] = {"QSO_UVLF_bhbol_nonstoch_z5.dat", "QSO_UVLF_bhbol_nonstoch_z6.dat", "QSO_UVLF_bhbol_nonstoch_z7.dat", "QSO_UVLF_bhbol_nonstoch_z8.dat", "QSO_UVLF_bhbol_nonstoch_z9.dat", "QSO_UVLF_bhbol_nonstoch_z10.dat", "QSO_UVLF_bhbol_nonstoch_z11.dat", "QSO_UVLF_bhbol_nonstoch_z12.dat", "QSO_UVLF_bhbol_nonstoch_z13.dat", "QSO_UVLF_bhbol_nonstoch_z15.dat", "QSO_UVLF_bhbol_nonstoch_z17.dat", "QSO_UVLF_bhbol_nonstoch_z19.dat"};
const char* name_list_gal[12] = {"gal_UVLF_z5.dat", "gal_UVLF_z6.dat", "gal_UVLF_z7.dat", "gal_UVLF_z8.dat", "gal_UVLF_z9.dat", "gal_UVLF_z10.dat", "gal_UVLF_z11.dat", "gal_UVLF_z12.dat", "gal_UVLF_z13.dat", "gal_UVLF_z15.dat", "gal_UVLF_z17.dat", "gal_UVLF_z19.dat"};
const char* name_list_uvlf_ueda1[12] = {"ueda1_UVLF_z5.dat", "ueda1_UVLF_z6.dat", "ueda1_UVLF_z7.dat", "ueda1_UVLF_z8.dat", "ueda1_UVLF_z9.dat", "ueda1_UVLF_z10.dat", "ueda1_UVLF_z11.dat", "ueda1_UVLF_z12.dat", "ueda1_UVLF_z13.dat", "ueda1_UVLF_z15.dat", "ueda1_UVLF_z17.dat", "ueda1_UVLF_z19.dat"};
const char* name_list_uvlf_ueda2[12] = {"ueda2_UVLF_z5.dat", "ueda2_UVLF_z6.dat", "ueda2_UVLF_z7.dat", "ueda2_UVLF_z8.dat", "ueda2_UVLF_z9.dat", "ueda2_UVLF_z10.dat", "ueda2_UVLF_z11.dat", "ueda2_UVLF_z12.dat", "ueda2_UVLF_z13.dat", "ueda2_UVLF_z15.dat", "ueda2_UVLF_z17.dat", "ueda2_UVLF_z19.dat"};
const char* name_list_uvlf_ueda3[12] = {"ueda3_UVLF_z5.dat", "ueda3_UVLF_z6.dat", "ueda3_UVLF_z7.dat", "ueda3_UVLF_z8.dat", "ueda3_UVLF_z9.dat", "ueda3_UVLF_z10.dat", "ueda3_UVLF_z11.dat", "ueda3_UVLF_z12.dat", "ueda3_UVLF_z13.dat", "ueda3_UVLF_z15.dat", "ueda3_UVLF_z17.dat", "ueda3_UVLF_z19.dat"};


int lenzlist = 12; 

// total uvlf with ueda's corrections - 1

for (int w = 0; w < lenzlist; w++) {

    double *nd12_fbd = new double[nouvsteps-1];
    double *err12 = new double[nouvsteps-1];
    double *ngal12 = new double[nouvsteps-1];
    double *erruo12_fbd = new double[nouvsteps-1];
    double *errdo12_fbd= new double[nouvsteps-1];
    
    for (int p = 0; p < nouvsteps-1; p++) {
        nd12_fbd[p] = 0.0; err12[p] = 0; ngal12[p] = 0; erruo12_fbd[p] = 0.0; errdo12_fbd[p] = 0.0;}

    for (i = 0; i < no_gal; i++) {
        
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == z_list[w]) {
            							
                for (k = 0; k < nouvsteps-1; k++) {
         
                    if ((abs(uv_obstot_ueda[j]) >= uvbin[k]) && (abs(uv_obstot_ueda[j]) < uvbin[k+1]) && (uv_obstot_ueda[j] < 0.0)) {
            
                        double nd12 = normbin * pow(10.0, nd_gal[i]);
                        nd12_fbd[k] = nd12_fbd[k] + nd12;
                        ngal12[k] = ngal12[k] + 1;
                        double temp = nd12 - nd12_fbd[k];
                        err12[k] = err12[k] + (temp * temp);}}}}}
                        
    ofstream fout1g(name_list_uvlf_ueda1[w]);        
    
    for (int x = 0; x < nouvsteps-1; x++) {

        nd12_fbd[x] = nd12_fbd[x] / stepuv;
        
        if (ngal12[x] > 0) 
            err12[x] = sqrt(err12[x]/ngal12[x]);
        else
            err12[x] = 0;
        
        if (ngal12[x] != 0) {
            erruo12_fbd[x] = log10(abs(nd12_fbd[x] + err12[x]));
            errdo12_fbd[x] = log10(abs(nd12_fbd[x] - err12[x]));}
        
        else {
            erruo12_fbd[x] = -100.0;
            errdo12_fbd[x] = -100.0;}
        
        if (nd12_fbd[x] > 0)
            nd12_fbd[x] = log10(nd12_fbd[x]);
        else
            nd12_fbd[x] = -100.0;
 
        fout1g<<(-1.0*mid_uvbin[x])<<"\t"<<nd12_fbd[x]<<"\t"<<erruo12_fbd[x]<<"\t"<<errdo12_fbd[x]<<"\t"<<err12[x]<<endl;}
            
    fout1g.close();
    
    delete []nd12_fbd;
    delete []err12;
    delete []ngal12;
    delete []erruo12_fbd;
    delete []errdo12_fbd;}

//************************************************************************
// AGN, QSO and galaxies UVLF
//************************************************************************

for (int w = 0; w < lenzlist; w++) {

    cout << z_list[w] << endl;

    double *nd12_fbd = new double[nouvsteps-1]; // without ueda correction
    double *nd13_fbd = new double[nouvsteps-1]; // with ueda correction
    double *err13 = new double[nouvsteps-1];
    double *ngal13 = new double[nouvsteps-1];
    double *erruo13_fbd = new double[nouvsteps-1];
    double *errdo13_fbd= new double[nouvsteps-1];
    int *count0 = new int[nouvsteps-1];
    int *count1 = new int[nouvsteps-1];
    
    for (int p = 0; p < nouvsteps-1; p++) {
        nd12_fbd[p] = 0.0; nd13_fbd[p] = 0.0; err13[p] = 0; ngal13[p] = 0; erruo13_fbd[p] = 0.0; errdo13_fbd[p] = 0.0; count0[p] = 0; count1[p] = 0;}

    for (i = 0; i < no_gal; i++) {
        
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == z_list[w]) {
         
                for (k = 0; k < nouvsteps-1; k++) {
         
                    if ((abs(uv_obsbh[j]) >= uvbin[k]) && (abs(uv_obsbh[j]) < uvbin[k+1]) && (uv_obsbh[j] < 0.0)) {
                    
                    	if (accr_flag[j] == 0)
                    		count0[k] = count0[k] + 1;
                		
                		else if (accr_flag[j] == 1)
                			count1[k] = count1[k] + 1;
            
                        double nd12 = normbin * pow(10.0, nd_gal[i]);
                        nd12_fbd[k] = nd12_fbd[k] + nd12;}

                    if ((abs(uv_obsbh_ueda[j]) >= uvbin[k]) && (abs(uv_obsbh_ueda[j]) < uvbin[k+1]) && (uv_obsbh_ueda[j] < 0.0)) {
                    
                    	if (accr_flag[j] == 0)
                    		count0[k] = count0[k] + 1;
                		
                		else if (accr_flag[j] == 1)
                			count1[k] = count1[k] + 1;

                        double nd13 = normbin * pow(10.0, nd_gal[i]);
                        nd13_fbd[k] = nd13_fbd[k] + nd13;
                        ngal13[k] = ngal13[k] + 1;
                        double temp = nd13 - nd13_fbd[k];
                        err13[k] = err13[k] + (temp * temp);}}}}}
                            
    ofstream fout2e(name_list_AGN[w]);        
    
    for (int x = 0; x < nouvsteps-1; x++) {

        nd12_fbd[x] = nd12_fbd[x] / stepuv;
        nd13_fbd[x] = nd13_fbd[x] / stepuv;
        
        if (ngal13[x] > 0) 
            err13[x] = sqrt(err13[x]/ngal13[x]);
        else
            err13[x] = 0;
        
        if (ngal13[x] != 0) {
            erruo13_fbd[x] = log10(abs(nd13_fbd[x] + err13[x]));
            errdo13_fbd[x] = log10(abs(nd13_fbd[x] - err13[x]));}
        
        else {
            erruo13_fbd[x] = -100.0;
            errdo13_fbd[x] = -100.0;}
        
        if (nd12_fbd[x] > 0)
            nd12_fbd[x] = log10(nd12_fbd[x]);
        else
            nd12_fbd[x] = -100.0;

        if (nd13_fbd[x] > 0)
            nd13_fbd[x] = log10(nd13_fbd[x]);
        else
            nd13_fbd[x] = -100.0;
		fout2e << (-1.0*mid_uvbin[x]) << "\t" << nd12_fbd[x] << "\t" << erruo13_fbd[x] << "\t" << errdo13_fbd[x] << "\t" << err13[x] << "\t" << count0[x] << "\t" << count1[x] << "\t" << nd13_fbd[x] << endl;}
            
    fout2e.close();
    
    delete []nd12_fbd;
    delete []nd13_fbd;
    delete []err13;
    delete []ngal13;
    delete []erruo13_fbd;
    delete []errdo13_fbd;
    delete []count0;
    delete []count1;}

// QSO UVLF selected by bolometric luminosity - stochastic obscuration correction

double mock_vol = pow(1000.0, 3);

for (int w = 0; w < lenzlist; w++) {

    cout << z_list[w] << endl;

    double *nd12_fbd_44 = new double[nouvsteps-1]; // without ueda correction
    double *nd12_fbd_45 = new double[nouvsteps-1]; // without ueda correction
    double *nd12_fbd_43 = new double[nouvsteps-1]; // without ueda correction
    double *nd13_fbd_44 = new double[nouvsteps-1]; // with ueda correction
    double *nd13_fbd_45 = new double[nouvsteps-1]; // with ueda correction
    double *nd13_fbd_43 = new double[nouvsteps-1]; // with ueda correction
    double *err13_44 = new double[nouvsteps-1];
    double *err13_45 = new double[nouvsteps-1];
    double *err13_43 = new double[nouvsteps-1];
    double *ngal13_44 = new double[nouvsteps-1];
    double *ngal13_45 = new double[nouvsteps-1];
    double *ngal13_43 = new double[nouvsteps-1];
    double *erruo13_fbd_44 = new double[nouvsteps-1];
    double *erruo13_fbd_45 = new double[nouvsteps-1];
    double *erruo13_fbd_43 = new double[nouvsteps-1];
    double *errdo13_fbd_44 = new double[nouvsteps-1];
    double *errdo13_fbd_45 = new double[nouvsteps-1];
    double *errdo13_fbd_43 = new double[nouvsteps-1];
    double *err12_44 = new double[nouvsteps-1];
    double *err12_45 = new double[nouvsteps-1];
    double *err12_43 = new double[nouvsteps-1];
    double *ngal12_44 = new double[nouvsteps-1];
    double *ngal12_45 = new double[nouvsteps-1];
    double *ngal12_43 = new double[nouvsteps-1];
    double *erruo12_fbd_44 = new double[nouvsteps-1];
    double *erruo12_fbd_45 = new double[nouvsteps-1];
    double *erruo12_fbd_43 = new double[nouvsteps-1];
    double *errdo12_fbd_44 = new double[nouvsteps-1];
    double *errdo12_fbd_45 = new double[nouvsteps-1];
    double *errdo12_fbd_43 = new double[nouvsteps-1];
    int *count0 = new int[nouvsteps-1];
    int *count1 = new int[nouvsteps-1];
    
    for (int p = 0; p < nouvsteps-1; p++) {
        nd12_fbd_44[p] = 0.0; nd12_fbd_45[p] = 0.0; nd12_fbd_43[p] = 0.0; nd13_fbd_44[p] = 0.0; nd13_fbd_45[p] = 0.0; nd13_fbd_43[p] = 0.0; err13_44[p] = 0; err13_45[p] = 0; err13_43[p] = 0; ngal13_44[p] = 0; ngal13_45[p] = 0; ngal13_43[p] = 0; erruo13_fbd_44[p] = 0.0; erruo13_fbd_45[p] = 0.0; erruo13_fbd_43[p] = 0.0; errdo13_fbd_44[p] = 0.0; errdo13_fbd_45[p] = 0.0; errdo13_fbd_43[p] = 0.0; err12_44[p] = 0; err12_45[p] = 0; err12_43[p] = 0; ngal12_44[p] = 0; ngal12_45[p] = 0; ngal12_43[p] = 0; erruo12_fbd_44[p] = 0.0; erruo12_fbd_45[p] = 0.0; erruo12_fbd_43[p] = 0.0; errdo12_fbd_44[p] = 0.0; errdo12_fbd_45[p] = 0.0; errdo12_fbd_43[p] = 0.0; count0[p] = 0; count1[p] = 0;}

    for (i = 0; i < no_gal; i++) {
        
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == z_list[w]) {
         
                for (k = 0; k < nouvsteps-1; k++) {
         
                    if ((abs(uv_obstot[j]) >= uvbin[k]) && (abs(uv_obstot[j]) < uvbin[k+1]) && (uv_obstot[j] < 0.0)) {

                        double nd12 = normbin * pow(10.0, nd_gal[i]) / stepuv;

						if (bhenergy[j] >= pow(10.0, 44)){                    
                        	nd12_fbd_44[k] = nd12_fbd_44[k] + nd12;
		                    ngal12_44[k] = ngal12_44[k] + nd12*mock_vol;}
                    	if (bhenergy[j] >= pow(10.0, 45)){   
                        	nd12_fbd_45[k] = nd12_fbd_45[k] + nd12;
		                    ngal12_45[k] = ngal12_45[k] + nd12*mock_vol;}
                    	if (bhenergy[j] >= pow(10.0, 43)){   
                        	nd12_fbd_43[k] = nd12_fbd_43[k] + nd12;
		                    ngal12_43[k] = ngal12_43[k] + nd12*mock_vol;}}

                    if ((abs(uv_obstot_ueda[j]) >= uvbin[k]) && (abs(uv_obstot_ueda[j]) < uvbin[k+1]) && (uv_obstot_ueda[j] < 0.0)) {
                    
                        double nd13 = normbin * pow(10.0, nd_gal[i]) / stepuv;
						if (bhenergy_ueda[j] >= pow(10.0, 44)) {                   
                        	nd13_fbd_44[k] = nd13_fbd_44[k] + nd13;
		                    ngal13_44[k] = ngal13_44[k] + nd13*mock_vol;}
                       	if (bhenergy_ueda[j] >= pow(10.0, 45)) { 
                        	nd13_fbd_45[k] = nd13_fbd_45[k] + nd13;
		                    ngal13_45[k] = ngal13_45[k] + nd13*mock_vol;}
                    	if (bhenergy_ueda[j] >= pow(10.0, 43)) {
                        	nd13_fbd_43[k] = nd13_fbd_43[k] + nd13;
		                    ngal13_43[k] = ngal13_43[k] + nd13*mock_vol;}}}}}}
                            
    ofstream fout2e(name_list_QSO_bhbol[w]);        
    
    for (int x = 0; x < nouvsteps-1; x++) {

        if (ngal12_44[x] > 0) 
            err12_44[x] = sqrt(ngal12_44[x])/mock_vol;
        else
            err12_44[x] = 0;
        if (ngal12_45[x] > 0) 
            err12_45[x] = sqrt(ngal12_45[x])/mock_vol;
        else
            err12_45[x] = 0;
        if (ngal12_43[x] > 0) 
            err12_43[x] = sqrt(ngal12_43[x])/mock_vol;
        else
            err12_43[x] = 0;
        
        if (ngal12_44[x] != 0) {
            erruo12_fbd_44[x] = log10(abs(nd12_fbd_44[x] + err12_44[x]));
            errdo12_fbd_44[x] = log10(abs(nd12_fbd_44[x] - err12_44[x]));}
        else {
            erruo12_fbd_44[x] = -100.0;
            errdo12_fbd_44[x] = -100.0;}
        if (ngal12_45[x] != 0) {
            erruo12_fbd_45[x] = log10(abs(nd12_fbd_45[x] + err12_45[x]));
            errdo12_fbd_45[x] = log10(abs(nd12_fbd_45[x] - err12_45[x]));}
        else {
            erruo12_fbd_45[x] = -100.0;
            errdo12_fbd_45[x] = -100.0;}
        if (ngal12_43[x] != 0) {
            erruo12_fbd_43[x] = log10(abs(nd12_fbd_43[x] + err12_43[x]));
            errdo12_fbd_43[x] = log10(abs(nd12_fbd_43[x] - err12_43[x]));}
        else {
            erruo12_fbd_43[x] = -100.0;
            errdo12_fbd_43[x] = -100.0;}


        if (ngal13_44[x] > 0) 
            err13_44[x] = sqrt(err13_44[x]/ngal13_44[x]);
        else
            err13_44[x] = 0;
        if (ngal13_45[x] > 0) 
            err13_45[x] = sqrt(err13_45[x]/ngal13_45[x]);
        else
            err13_45[x] = 0;
        if (ngal13_43[x] > 0) 
            err13_43[x] = sqrt(err13_43[x]/ngal13_43[x]);
        else
            err13_43[x] = 0;
        
        if (ngal13_44[x] != 0) {
            erruo13_fbd_44[x] = log10(abs(nd13_fbd_44[x] + err13_44[x]));
            errdo13_fbd_44[x] = log10(abs(nd13_fbd_44[x] - err13_44[x]));}
        else {
            erruo13_fbd_44[x] = -100.0;
            errdo13_fbd_44[x] = -100.0;}
        if (ngal13_45[x] != 0) {
            erruo13_fbd_45[x] = log10(abs(nd13_fbd_45[x] + err13_45[x]));
            errdo13_fbd_45[x] = log10(abs(nd13_fbd_45[x] - err13_45[x]));}
        else {
            erruo13_fbd_45[x] = -100.0;
            errdo13_fbd_45[x] = -100.0;}
        if (ngal13_43[x] != 0) {
            erruo13_fbd_43[x] = log10(abs(nd13_fbd_43[x] + err13_43[x]));
            errdo13_fbd_43[x] = log10(abs(nd13_fbd_43[x] - err13_43[x]));}
        else {
            erruo13_fbd_43[x] = -100.0;
            errdo13_fbd_43[x] = -100.0;}
        
        
        if (nd12_fbd_44[x] > 0)
            nd12_fbd_44[x] = log10(nd12_fbd_44[x]);
        else
            nd12_fbd_44[x] = -100.0;
        if (nd12_fbd_45[x] > 0)
            nd12_fbd_45[x] = log10(nd12_fbd_45[x]);
        else
            nd12_fbd_45[x] = -100.0;
        if (nd12_fbd_43[x] > 0)
            nd12_fbd_43[x] = log10(nd12_fbd_43[x]);
        else
            nd12_fbd_43[x] = -100.0;

        if (nd13_fbd_44[x] > 0)
            nd13_fbd_44[x] = log10(nd13_fbd_44[x]);
        else
            nd13_fbd_44[x] = -100.0;
        if (nd13_fbd_45[x] > 0)
            nd13_fbd_45[x] = log10(nd13_fbd_45[x]);
        else
            nd13_fbd_45[x] = -100.0;
        if (nd13_fbd_43[x] > 0)
            nd13_fbd_43[x] = log10(nd13_fbd_43[x]);
        else
            nd13_fbd_43[x] = -100.0;
 
        fout2e << (-1.0*mid_uvbin[x]) << "\t" << nd12_fbd_44[x] << "\t" << nd12_fbd_45[x] << "\t" << nd12_fbd_43[x] << "\t" <<erruo12_fbd_44[x] << "\t" <<erruo12_fbd_45[x] << "\t" <<erruo12_fbd_43[x] << "\t" << errdo12_fbd_44[x] << "\t" << errdo12_fbd_45[x] << "\t" << errdo12_fbd_43[x] << "\t" << err12_44[x] << "\t" << err12_45[x] << "\t" << err12_43[x] << "\t" <<erruo13_fbd_44[x] << "\t" <<erruo13_fbd_45[x] << "\t" <<erruo13_fbd_43[x] << "\t" << errdo13_fbd_44[x] << "\t" << errdo13_fbd_45[x] << "\t" << errdo13_fbd_43[x] << "\t" << err13_44[x] << "\t" << err13_45[x] << "\t" << err13_43[x] << "\t" << count0[x] << "\t" << count1[x] << "\t" << nd13_fbd_44[x] << "\t" << nd13_fbd_45[x] << "\t" << nd13_fbd_43[x] << endl;}
            
    fout2e.close();
    
    delete []nd12_fbd_44;
    delete []nd12_fbd_45;
    delete []nd12_fbd_43;
    delete []nd13_fbd_44;
    delete []nd13_fbd_45;
    delete []nd13_fbd_43;
    delete []err13_44;
    delete []err13_45;
    delete []err13_43;
    delete []ngal13_44;
    delete []ngal13_45;
    delete []ngal13_43;
    delete []erruo13_fbd_44;
    delete []erruo13_fbd_45;
    delete []erruo13_fbd_43;
    delete []errdo13_fbd_44;
    delete []errdo13_fbd_45;
    delete []errdo13_fbd_43;
    delete []err12_44;
    delete []err12_45;
    delete []err12_43;
    delete []ngal12_44;
    delete []ngal12_45;
    delete []ngal12_43;
    delete []erruo12_fbd_44;
    delete []erruo12_fbd_45;
    delete []erruo12_fbd_43;
    delete []errdo12_fbd_44;
    delete []errdo12_fbd_45;
    delete []errdo12_fbd_43;
    delete []count0;
    delete []count1;}

// QSO UVLF selected by bolometric luminosity - non-corrected

for (int w = 0; w < lenzlist; w++) {

    cout << z_list[w] << endl;

    double *nd12_fbd_43 = new double[nobhbolsteps-1]; // without ueda correction
    double *nd13_fbd_43 = new double[nobhbolsteps-1]; // with ueda correction
    double *err13_43 = new double[nobhbolsteps-1];
    double *ngal13_43 = new double[nobhbolsteps-1];
    double *erruo13_fbd_43 = new double[nobhbolsteps-1];
    double *errdo13_fbd_43 = new double[nobhbolsteps-1];
    double *err12_43 = new double[nobhbolsteps-1];
    double *ngal12_43 = new double[nobhbolsteps-1];
    double *erruo12_fbd_43 = new double[nobhbolsteps-1];
    double *errdo12_fbd_43 = new double[nobhbolsteps-1];
    int *count0 = new int[nobhbolsteps-1];
    int *count1 = new int[nobhbolsteps-1];
    
    for (int p = 0; p < nobhbolsteps-1; p++) {
        nd12_fbd_43[p] = 0.0; nd13_fbd_43[p] = 0.0; err13_43[p] = 0; ngal13_43[p] = 0; erruo13_fbd_43[p] = 0.0; errdo13_fbd_43[p] = 0.0; err12_43[p] = 0; ngal12_43[p] = 0; erruo12_fbd_43[p] = 0.0; errdo12_fbd_43[p] = 0.0; count0[p] = 0; count1[p] = 0;}

    for (i = 0; i < no_gal; i++) {
        
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == z_list[w]) {
         
                for (k = 0; k < nobhbolsteps-1; k++) {
         
                    if ((log10(abs(bhenergy[j])) >= bhbolbin[k]) && (log10(abs(bhenergy[j])) < bhbolbin[k+1]) && (uv_obstot[j] < 0.0)) {

                        double nd12 = normbin * pow(10.0, nd_gal[i]) / stepbhbol;

                    	nd12_fbd_43[k] = nd12_fbd_43[k] + nd12;
	                    ngal12_43[k] = ngal12_43[k] + nd12*mock_vol;}

                    if ((log10(abs(bhenergy_ueda[j])) >= bhbolbin[k]) && (log10(abs(bhenergy_ueda[j])) < bhbolbin[k+1]) && (uv_obstot_ueda[j] < 0.0)) {
                    
                        double nd13 = normbin * pow(10.0, nd_gal[i]) / stepbhbol;
                    	nd13_fbd_43[k] = nd13_fbd_43[k] + nd13;
	                    ngal13_43[k] = ngal13_43[k] + nd13*mock_vol;}}}}}
                                                        
    ofstream fout2e(name_list_QSO_bhbol_nonstoch[w]);        
    
    for (int x = 0; x < nobhbolsteps-1; x++) {

        if (ngal12_43[x] > 0) 
            err12_43[x] = sqrt(ngal12_43[x])/mock_vol;
        else
            err12_43[x] = 0;
        
        if (ngal12_43[x] != 0) {
            erruo12_fbd_43[x] = log10(abs(nd12_fbd_43[x] + err12_43[x]));
            errdo12_fbd_43[x] = log10(abs(nd12_fbd_43[x] - err12_43[x]));}
        else {
            erruo12_fbd_43[x] = -100.0;
            errdo12_fbd_43[x] = -100.0;}


        if (ngal13_43[x] > 0) 
            err13_43[x] = sqrt(err13_43[x]/ngal13_43[x]);
        else
            err13_43[x] = 0;
        
        if (ngal13_43[x] != 0) {
            erruo13_fbd_43[x] = log10(abs(nd13_fbd_43[x] + err13_43[x]));
            errdo13_fbd_43[x] = log10(abs(nd13_fbd_43[x] - err13_43[x]));}
        else {
            erruo13_fbd_43[x] = -100.0;
            errdo13_fbd_43[x] = -100.0;}
        
        
        if (nd12_fbd_43[x] > 0)
            nd12_fbd_43[x] = log10(nd12_fbd_43[x]);
        else
            nd12_fbd_43[x] = -100.0;

        if (nd13_fbd_43[x] > 0)
            nd13_fbd_43[x] = log10(nd13_fbd_43[x]);
        else
            nd13_fbd_43[x] = -100.0;
 
        fout2e << mid_bhbolbin[x] << "\t" << nd12_fbd_43[x] << "\t" << erruo12_fbd_43[x] << "\t" << errdo12_fbd_43[x] << "\t" << err12_43[x] << "\t" << erruo13_fbd_43[x] << "\t" << errdo13_fbd_43[x] << "\t" << err13_43[x] << "\t" << count0[x] << "\t" << count1[x] << "\t" << nd13_fbd_43[x] << endl;}
            
    fout2e.close();
    
    delete []nd12_fbd_43;
    delete []nd13_fbd_43;
    delete []err13_43;
    delete []ngal13_43;
    delete []erruo13_fbd_43;
    delete []errdo13_fbd_43;
    delete []err12_43;
    delete []ngal12_43;
    delete []erruo12_fbd_43;
    delete []errdo12_fbd_43;
    delete []count0;
    delete []count1;}


for (int w = 0; w < lenzlist; w++) {

    cout << z_list[w] << endl;

    double *nd12_fbd = new double[nouvsteps-1];
    double *err12 = new double[nouvsteps-1];
    double *ngal12 = new double[nouvsteps-1];
    double *erruo12_fbd = new double[nouvsteps-1];
    double *errdo12_fbd= new double[nouvsteps-1];
    

    for (int p = 0; p < nouvsteps-1; p++) {
        nd12_fbd[p] = 0.0; err12[p] = 0; ngal12[p] = 0; erruo12_fbd[p] = 0.0; errdo12_fbd[p] = 0.0;}

    for (i = 0; i < no_gal; i++) {
        
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == z_list[w]) {
         
                for (k = 0; k < nouvsteps-1; k++) {
         
                    if ((abs(uv_obsfb[j]) >= uvbin[k]) && (abs(uv_obsfb[j]) < uvbin[k+1]) && (uv_obsfb[j] < 0.0)) {
            
                        double nd12 = normbin * pow(10.0, nd_gal[i]);
                        nd12_fbd[k] = nd12_fbd[k] + nd12;
                        ngal12[k] = ngal12[k] + 1;
                        double temp = nd12 - nd12_fbd[k];
                        err12[k] = err12[k] + (temp * temp);}}}}}
                        
    
    ofstream fout3e(name_list_gal[w]);        
    
    for (int x = 0; x < nouvsteps-1; x++) {

        nd12_fbd[x] = nd12_fbd[x] / stepuv;
        
        if (ngal12[x] > 0) 
            err12[x] = sqrt(err12[x]/ngal12[x]);
        else
            err12[x] = 0;
        
        if (ngal12[x] != 0) {
            erruo12_fbd[x] = log10(abs(nd12_fbd[x] + err12[x]));
            errdo12_fbd[x] = log10(abs(nd12_fbd[x] - err12[x]));}
        
        else {
            erruo12_fbd[x] = -100.0;
            errdo12_fbd[x] = -100.0;}
        
        if (nd12_fbd[x] > 0)
            nd12_fbd[x] = log10(nd12_fbd[x]);
        else
            nd12_fbd[x] = -100.0;
 
        fout3e<<(-1.0*mid_uvbin[x])<<"\t"<<nd12_fbd[x]<<"\t"<<erruo12_fbd[x]<<"\t"<<errdo12_fbd[x]<<"\t"<<err12[x]<<endl;}
            
    fout3e.close();
    
    delete []nd12_fbd;
    delete []err12;
    delete []ngal12;
    delete []erruo12_fbd;
    delete []errdo12_fbd;}


// ***********************************************************************
// printing all the parameters of all the halos at the selected redshifts
// ***********************************************************************

for (i = 0; i < no_gal; i++)

    for (j = plcb_gal[i]; j < plce_gal[i]; j++)
        
        if (lconts_int[j] > 0)
            
            lconts_int[j] = log10(lconts_int[j]);
        
        else
        
            lconts_int[j] = -100.0;
      
const char* name_list2[13] = {"phyprop_all_z4start_nofb_fstareff_mergertree_z1.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z2.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z3.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z4.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z5.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z6.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z7.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z8.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z9.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z10.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z11.dat", "phyprop_all_z4start_nofb_fstareff_mergertree_z12.dat"};
const char* name_list3[13] = {"bhprop_all_mres10to9_mmax10to15_sedd_z1.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z2.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z3.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z4.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z5.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z6.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z7.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z8.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z9.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z10.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z11.dat", "bhprop_all_mres10to9_mmax10to15_sedd_z12.dat"};


double lenzlist2[13] = {z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12};

ofstream foutmacchianera("main_branch_allgalaxies.dat");

for (i = galnmin; i < galnmax; i++) {

    for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
    
    	if (main_branch_flag[j] == 1) {
	
        foutmacchianera	<< zp_mt[j] << "\t" << log10(mh_int[j]) << "\t" << mhp_mt[j] << "\t" << mh_int[j]/pow(10.0, mhp_mt[j]) << "\t" << mgacc_int[j] << "\t" << mgmerge_int[j] << "\t" << mstarfb_int[j] << "\t" << mgf_int[j] << "\t" << mgacc_new[j] << "\t" << mgmerge_new[j] << "\t" << mstarfb_new[j] << "\t" << fstarej[j] << "\t" << sf_eff << "\t" << fstareff[j] << "\t" << lconts_int[j] << "\t" << uv_obstot[j] << "\t" << nprog[j] << "\t" << bhflag[j] << "\t" << mgout_sf_int[j] << "\t" << mgout_agn_int[j] << "\t" << bhenergy[j] << "\t" << idp_mt[j] << "\t" << id_parent_mt[j] << "\t" << id_main_prog[j] << "\t" << bhmain[j] << "\t" << nprog_bh[j] << "\t" << flag_accr[j] << "\t" << nd_gal[i] << "\t" << accrat[j] << "\t" << mass_bh_merged_new[j] << "\t" << uvflux_rat[j] << "\t" << bhaccmass[j] << "\t" << bh_L1375[j] << "\t" << mgout_sf_new[j] << "\t" << mgout_agn_new[j] << "\t" << id_df_merge[j] << "\t" << uv_obsbh[j] << "\t" << uv_obsfb[j] << "\t" << accr_flag[j] << "\t" << mass_bh_merged[j] << "\t" << mdm_merged_int[j] << "\t" << mdm_merged_new[j] << "\t" << mstar_merged_int[j] << "\t" << mstar_merged_new[j] << "\t" << bhenergy_ueda[j] << "\t" << uvflux_rat_ueda[j] << "\t" << bh_L1375_ueda[j] << "\t" << uv_obstot_ueda[j] << "\t" << uv_obsbh_ueda[j] << "\t" << bondi_acc_mass_int[j] << "\t" << bondi_acc_mass_new[j] << "\t" << radio_lum[j] << "\t" << mgi_cold[j] << "\t" << mgf_cold[j] << "\t" << mgi_hot[j] << "\t" << mgf_hot[j] << "\t" << mgcool_new[j] << "\t" << mgcool_int[j] << "\t" << time_since_mm[j] << "\t" << gasfrac_i[j] << "\t" << jet_power[j] << "\t" << jet_uvlum[j] << "\t" << mg_heated[j] << "\t" << met_cold[j] << "\t" << met_hot[j] << "\t" << met_res[j] << "\t" << met_star[j] << "\t" << z_cold_i[j] << "\t" << z_hot_i[j] << "\t" << z_res_i[j] << "\t" << mg_res[j] << endl;}}}

foutmacchianera.close();

for (int w = 0; w < 13; w++) {

    ofstream foutpippo(name_list2[w]);

    for (i = galnmin; i < galnmax; i++) {
    
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if (zp_mt[j] == lenzlist2[w]) {
            
		foutpippo << zp_mt[j] << "\t" << log10(mh_int[j]) << "\t" << mhp_mt[j] << "\t" << mh_int[j]/pow(10.0, mhp_mt[j]) << "\t" << mgacc_int[j] << "\t" << mgmerge_int[j] << "\t" << mstarfb_int[j] << "\t" << mgf_int[j] << "\t" << mgacc_new[j] << "\t" << mgmerge_new[j] << "\t" << mstarfb_new[j] << "\t" << fstarej[j] << "\t" << sf_eff << "\t" << fstareff[j] << "\t" << lconts_int[j] << "\t" << uv_obstot[j] << "\t" << nprog[j] << "\t" << bhflag[j] << "\t" << mgout_sf_int[j] << "\t" << mgout_agn_int[j] << "\t" << bhenergy[j] << "\t" << idp_mt[j] << "\t" << id_parent_mt[j] << "\t" << id_main_prog[j] << "\t" << bhmain[j] << "\t" << nprog_bh[j] << "\t" << flag_accr[j] << "\t" << nd_gal[i] << "\t" << accrat[j] << "\t" << mass_bh_merged_new[j] << "\t" << uvflux_rat[j] << "\t" << bhaccmass[j] << "\t" << bh_L1375[j] << "\t" << mgout_sf_new[j] << "\t" << mgout_agn_new[j] << "\t" << id_df_merge[j] << "\t" << uv_obsbh[j] << "\t" << uv_obsfb[j] << "\t" << accr_flag[j] << "\t" << mass_bh_merged[j] << "\t" << mdm_merged_int[j] << "\t" << mdm_merged_new[j] << "\t" << mstar_merged_int[j] << "\t" << mstar_merged_new[j] << "\t" << bhenergy_ueda[j] << "\t" << uvflux_rat_ueda[j] << "\t" << bh_L1375_ueda[j] << "\t" << uv_obstot_ueda[j] << "\t" << uv_obsbh_ueda[j] << "\t" << flag[j] << "\t" << rvir[j] << "\t" << bondi_acc_mass_int[j] << "\t" << bondi_acc_mass_new[j] << "\t" << radio_lum[j] << "\t" << mgi_cold[j] << "\t" << mgf_cold[j] << "\t" << mgi_hot[j] << "\t" << mgf_hot[j] << "\t" << mgcool_new[j] << "\t" << mgcool_int[j] << "\t" << time_since_mm[j] << "\t" << gasfrac_i[j] << "\t" << jet_power[j] << "\t" << jet_uvlum[j] << "\t" << mg_heated[j] << "\t" << met_cold[j] << "\t" << met_hot[j] << "\t" << met_res[j] << "\t" << met_star[j] << "\t" << z_cold_i[j] << "\t" << z_hot_i[j] << "\t" << z_res_i[j]  << "\t" << nuLnu_adaf[j] << "\t" << nuLnu_td[j] << "\t" << radio_lum_adaf[j] << "\t" << radio_lum_td[j] << "\t" << mg_res[j] << endl;}}}
                 
     foutpippo.close();}

for (int w = 0; w < 13; w++) {

    ofstream foutpluto(name_list3[w]);

    for (i = galnmin; i < galnmax; i++) {
    
        for (j = plcb_gal[i]; j < plce_gal[i]; j++) {
        
            if ((zp_mt[j] == lenzlist2[w]) && (bhmain[j] > 0.0)) {
            
        foutpluto << mhp_mt[j] << "\t" << nd_gal[i] << "\t" << bhmain[j]-bhaccmass[j] << "\t" << bhaccmass[j] << "\t" << accrat[j] << "\t" << jet_power[j] << "\t" << bhenergy[j]/bhaccmass[j]/mass_sun*tstep/pow(c,2) << endl;}}}
                 
     foutpluto.close();}

ftime = omp_get_wtime();
exec_time = ftime - itime;
printf("\n\nTime taken is %f", exec_time);     

return 0;} 

