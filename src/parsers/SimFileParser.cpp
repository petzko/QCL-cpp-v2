/*
 * SimFileParser.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */


#include "SimFileParser.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

#define MAX_CHARS_PER_LINE 512

int parseFile(char * filename,SimSettings* simset){


	//		simset->simname = *get_option_from_file(filename,"simname");
	//		//shall I include spatial hole burning or not
	//		simset->shb = atoi(get_option_from_file(filename,"shb")[0]);
	//		//characteristic time... e.g. seconds per picosecond and length meters per millimeter
	//		simset->tch = atof(get_option_from_file(filename,"tch")[0]);
	//		simset->lch = atof(get_option_from_file(filename,"lch")[0]);
	//		//grid size in x direction
	//		simset->N = atoi(get_option_from_file(filename,"N")[0]);
	//		//diffusion constant in [cm^2/sec]
	//		simset->D = atof(get_option_from_file(filename,"D")[0]);
	//		//dispersion constant in [ps^2/mm]!
	//		simset->disp = atof(get_option_from_file(filename,"disp")[0]);
	//
	//		//cavity length (mm)
	//		simset->Ltot = atof(get_option_from_file(filename,"Ltot")[0]);
	//		//period length [nm]
	//		simset->Lp = atof(get_option_from_file(filename,"Lp")[0]);
	//
	//		//doping density in [cm^-3]
	//		simset->dN = atof(get_option_from_file(filename,"dN")[0]);
	//		//doping region thickness [nm] (set equal to Lp if average doping density is known)
	//		simset->Ld = atof(get_option_from_file(filename,"Ld")[0]);
	//		//mode overlap factor [dimensionless]
	//		simset->Overlap = atof(get_option_from_file(filename,"Overlap")[0]);
	//		//refractive index of THz and GHz waves !
	//		simset->nTHz = atof(get_option_from_file(filename,"nTHz")[0]);
	//		simset->nRF = atof(get_option_from_file(filename,"nRF")[0]);
	//		//modulation amplitude factor
	//		simset->modA = atof(get_option_from_file(filename,"modA")[0]);
	//		//modulation frequency factor as a fraction of the RT freq.
	//		simset->modF = atof(get_option_from_file(filename,"modF")[0]);
	//		//initial bias value!
	//		simset->bias = atof(get_option_from_file(filename,"bias")[0]);
	//		//initial current value (A/mm)
	//		simset->current = atof(get_option_from_file(filename,"current")[0]);
	//		//voltage applied to the laser (in volts)
	//		simset->voltage = atof(get_option_from_file(filename,"voltage")[0]);
	//		//pure dephasing on/off
	//		simset->deph = atoi(get_option_from_file(filename,"deph")[0]);
	//
	//		//pure dephasing  time for the inj -> ull transition[ps] if any
	//		simset->Tdeph_1 = atof(get_option_from_file(filename,"Tdeph_1")[0]);
	//		//pure dephasing  times for the inj -> lll and ull->lll transition[ps] if any
	//		simset->Tdeph_2 = atof(get_option_from_file(filename,"Tdeph_2")[0]);
	//		//pure dephasing  times for the inj -> lll and ull->lll transition[ps] if any
	//		simset->Tdeph_3 = atof(get_option_from_file(filename,"Tdeph_3")[0]);
	//		//cavity loss in 1/cm;  (from optica paper)
	//		simset->loss = atof(get_option_from_file(filename,"loss")[0]);
	//
	//
	//		//################# simuation parameters #################
	//		simset->simRT = atoi(get_option_from_file(filename,"simRT")[0]);
	//		simset->plotCtr = atoi(get_option_from_file(filename,"plotCtr")[0]);
	//		simset->recordRT = atoi(get_option_from_file(filename,"recordRT")[0]);
	//		simset->nrSteps = atoi(get_option_from_file(filename,"nrSteps")[0]);
	//
	//		simset->zUL = atoi(get_option_from_file(filename,"zUL")[0]);

	double* HTB;// TB-hamiltonian
	double* Wmtx;// scattering rates mtx

}


std::string normalizeStr(std::string str){

	std::string result(str);
	result.erase(0,str.find_first_not_of(" \n\r\t"));
	result.erase(str.find_last_not_of(" \n\r\t")+1);

	std::locale loc;
	for (std::string::size_type i = 0; i < str.length(); ++i)
		result[i] = std::tolower(str[i], loc);

	return result;

}

std::vector<std::string> tokenize(std::string str, std::string DELIMITER){

	std::vector<std::string> v;
	int ctr_0 =0;
	bool tokenfound = false;

	//trim and convert to lower case!
	//        str = normalizeStr(str);

	for(int i = 0; i< str.length(); i++){
		if(DELIMITER.find(str[i]) != std::string::npos){
			if (!tokenfound){
				v.push_back(str.substr(ctr_0,i-ctr_0));
				tokenfound = true;
			}
			ctr_0 = i+1;
		}else{tokenfound = false;}
	}
	v.push_back(str.substr(ctr_0));

	return v;

}


std::vector<std::string>  get_option_from_file(char* filename, std::string option){

	std::ifstream file(filename);
	std::string line;


	std::vector<std::string> result;
	bool found  = false;

	// skip through the file until the first non commented or empty line
	if (file.is_open()){

		std::string DELIMITER = " =\n\t\r";

		while(file.good()){

			getline(file, line);
			while ((line.empty() || line.find("%#") !=std::string::npos) && file.good())
				getline(file, line);


			//remove trailing and leding spaces. convert to lower case!
			line = normalizeStr(line);
			// now parse line looking for your option ...
			//tokenize the string removing empty spaces

			std::vector<std::string> tokens = tokenize(line,DELIMITER);
			if(tokens.size() !=0)
				if( tokens[0].compare(option) ){
					found = true;
					for (int i=1; i<tokens.size();i++)
						result.push_back(tokens[i]);
					std::cout << "Option: " << option << " found. Value is:\n";
					for_each(result.begin(),result.end(),[](std::string str){std::cout << str << " "; });
					break;
				}
		}
		file.close();

	}
}


