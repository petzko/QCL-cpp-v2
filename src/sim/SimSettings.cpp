/*
 * SimSettings.cpp
 *
 *  Created on: Apr 11, 2016
 *      Author: petzko
 */

#include "SimSettings.hpp"


MB::SimSettings::SimSettings(char* simFName,char* setFName){

	this->simFileName = simFName;
	this->setFileName = setFName;
}

MB::SimSettings::~SimSettings(){

}

void MB::SimSettings::initSimSettings(){

	std::vector<std::string> res(0);

	res =  getOptionFromFile(this->setFileName,"scenario");
//	std::stringstream(res[0]) >> simname;

	//shall I include spatial hole burning or not
	res =  getOptionFromFile(this->setFileName,"shb");
	std::stringstream(res[0]) >> shb;
	//characteristic time... e.g. seconds per picosecond and length meters per millimeter
	res =  getOptionFromFile(this->setFileName,"tch");
	std::stringstream(res[0]) >> tch;

	res =  getOptionFromFile(this->setFileName,"lch");
	std::stringstream(res[0]) >> lch;

	//grid size in x direction
	res =  getOptionFromFile(this->setFileName,"N");
	std::stringstream(res[0]) >> N;
	//diffusion constant in [cm^2/sec]
	res =  getOptionFromFile(this->setFileName,"D");
	std::stringstream(res[0]) >> D;
	//dispersion constant in [ps^2/mm]!
	res =  getOptionFromFile(this->setFileName,"disp");
	std::stringstream(res[0]) >> disp;


	//cavity length (mm)
	res =  getOptionFromFile(this->setFileName,"Ltot");
	std::stringstream(res[0]) >> Ltot;
	//period length [nm]
	res =  getOptionFromFile(this->setFileName,"Lp");
	std::stringstream(res[0]) >> Lp;

	//doping density in [cm^-3]
	res =  getOptionFromFile(this->setFileName,"dN");
	std::stringstream(res[0]) >> dN;

	//doping region thickness [nm] (set equal to Lp if average doping density is known)
	res =  getOptionFromFile(this->setFileName,"Ld");
	std::stringstream(res[0]) >> Ld;

	//mode overlap factor [dimensionless]
	res =  getOptionFromFile(this->setFileName,"Overlap");
	std::stringstream(res[0]) >> Overlap;

	//			//refractive index of THz and GHz waves !
	res =  getOptionFromFile(this->setFileName,"nTHz");
	std::stringstream(res[0]) >> nTHz;
	res =  getOptionFromFile(this->setFileName,"nRF");
	std::stringstream(res[0]) >> nRF;


	//modulation amplitude factor
	res =  getOptionFromFile(this->setFileName,"modA");
	std::stringstream(res[0]) >> modA;

	//modulation frequency factor as a fraction of the RT freq.

	res =  getOptionFromFile(this->setFileName,"modF");
	std::stringstream(res[0]) >> modF;
	//initial bias value!
	res =  getOptionFromFile(this->setFileName,"bias");
	std::stringstream(res[0]) >> bias;

	//initial current value (A/mm)
	res =  getOptionFromFile(this->setFileName,"current");
	std::stringstream(res[0]) >> current;

	//voltage applied to the laser (in volts)
	res =  getOptionFromFile(this->setFileName,"voltage");
	std::stringstream(res[0]) >> voltage;


	//pure dephasing on/off
	res =  getOptionFromFile(this->setFileName,"deph");
	std::stringstream(res[0]) >> deph;

	//pure dephasing  time for the inj -> ull transition[ps] if any
	res =  getOptionFromFile(this->setFileName,"Tdeph_1");
	std::stringstream(res[0]) >> Tdeph_1;

	//pure dephasing  times for the inj -> lll and ull->lll transition[ps] if any
	res =  getOptionFromFile(this->setFileName,"Tdeph_2");
	std::stringstream(res[0]) >> Tdeph_2;
	res =  getOptionFromFile(this->setFileName,"Tdeph_3");
	std::stringstream(res[0]) >> Tdeph_3;

	//cavity loss in 1/cm;  (from optica paper)
	res =  getOptionFromFile(this->setFileName,"loss");
	std::stringstream(res[0]) >> loss;

	//################# simuation parameters #################
	res =  getOptionFromFile(this->setFileName,"simRT");
	std::stringstream(res[0]) >> simRT;
	res =  getOptionFromFile(this->setFileName,"plotCtr");
	std::stringstream(res[0]) >> plotCtr;
	res =  getOptionFromFile(this->setFileName,"recordRT");
	std::stringstream(res[0]) >> recordRT;
	res =  getOptionFromFile(this->setFileName,"nr_steps");
	std::stringstream(res[0]) >> nrSteps;


	res =  getOptionFromFile(this->simFileName,"zUL");
	std::stringstream(res[0]) >> zUL;

	res =  getOptionFromFile(this->simFileName,"HTB");
	HTB.resize(res.size());
	for (unsigned int i = 0; i < res.size();i++)
		std::stringstream(res[i]) >> HTB[i];

	res =  getOptionFromFile(this->simFileName,"Wmtx");
	Wmtx.resize(res.size());
	for (unsigned int i = 0; i < res.size();i++)
		std::stringstream(res[i]) >> Wmtx[i];



}
