/*
 * SimSettings.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#ifndef SRC_SIM_SIMSETTINGS_HPP_
#define SRC_SIM_SIMSETTINGS_HPP_

#include "../common/utils.hpp"
#include "../parsers/SimFileParser.hpp"

namespace MB{

struct SimSettings{

public:



	//simulattion name
	std::string simname;
	//shall I include spatial hole burning or not
	int shb;
	//characteristic time... e.g. seconds per picosecond and length meters per millimeter
	double tch, lch;
	//grid size in x direction
	int N;


	//diffusion constant in [cm^2/sec]
	double D;
	//dispersion constant in [ps^2/mm]!
	double disp;
	//cavity length (mm)
	double Ltot;
	//period length [nm]
	double Lp;

	//doping density in [cm^-3]
	double dN;
	//doping region thickness [nm] (set equal to Lp if average doping density is known)
	double Ld;
	//mode overlap factor [dimensionless]
	double Overlap;

	//refractive index of THz and GHz waves !
	double nTHz, nRF;
	//modulation amplitude factor
	double modA;
	//modulation frequency factor as a fraction of the RT freq.
	double modF;
	//initial bias value!
	double bias;
	//initial current value (A/mm)
	double current;
	//voltage applied to the laser (in volts)
	double voltage;
	//pure dephasing on/off
	int deph = 1;
	//pure dephasing  time for the inj -> ull transition[ps] if any
	double Tdeph_1;
	//pure dephasing  times for the inj -> lll and ull->lll transition[ps] if any
	double Tdeph_2;
	//pure dephasing  times for the inj -> lll and ull->lll transition[ps] if any
	double Tdeph_3;
	//cavity loss in 1/cm;  (from optica paper)
	double loss;

	//################# simuation parameters #################
	//number of round trips to be simulated (dimensionless)
	int simRT;
	int plotCtr;
	//number of round trips to record the data for (should be smaller than the total simRT);
	int recordRT;
	//number of steps for the MS solver (choose 5):
	int nrSteps;

	double zUL; // dipole matrix element
	std::vector<double> HTB;// TB-hamiltonian
	std::vector<double> Wmtx;// scattering rates mtx



	char* setFileName;
	char* simFileName;
	SimSettings(char* simFName,char* setFName);
	void initSimSettings();
	~SimSettings();
};



}

#endif /* SRC_SIM_SIMSETTINGS_HPP_ */
