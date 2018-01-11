#ifndef PARAMETER_H
#define PARAMETER_H

#include <math.h>

const double e0 = 8.8541878e-21;	//Vacuum Permitivity [F/nm]
const double q0 = 1.60217733e-19;	//Elementary Charge [Coulomb]
const double Avogadro = 6.0221409e23;

const double nm = 1e-9;		//Nanometer[m]

const double hbar = 1.0545718 * pow(10,-34); // m^2*kg/s or J*S
const double hbar_eV = 6.582119514 * pow(10,-16); // eV*s
const double h = 6.62606957 * pow(10,-34); // m^2*kg/s or J*S
const double h_eV = 4.135667516 * pow(10,-15); //eV*S
const double kb=1.380650E-23; //Boltzmann constant (J/K)
const double Tamb=300; //temperatureK
const double VT=kb*Tamb/q0; // V ~0.0259 v
const double C_speed=3 * pow(10,8); // speed of light m/s

const double m0 = 9.109 * pow(10,-31); // electron mass kg
const double m0_eV = 0.511e6; // eV
const double EBs = 3.34; // Barrier high from Si (eV)
const double EBe = 3.5; // Barrier high from Electrode (eV)

// K value
const double SiO2_permi=3.9;
const double Si_permi=11.7;
const double ZnO_permi=8.5;
const double Water_permi=80;
const double Ta2O5_permi=25;

const double eta = 1e-21;       // N*s/nm^2, Viscousity in nm
const double eta_m=eta*1e18;
const double Drho= 1e-24;       // kg/nm^3 , Density. 1g/cm3=1e-3g/m3=1e-6kg/m3=1e-24kg/m3
const double Si_Eg=1.12;
const double Si_chi=4.05;
const double Si_me=1.09;
const double Si_mh=1.15;
const double ZnO_Eg=3.37;
const double ZnO_me=0.318;
const double ZnO_mh=0.5;
const double ZnO_chi=3.37;

#endif
