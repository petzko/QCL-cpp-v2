///*
// * SimLoop.hpp
// *
// *  Created on: Apr 11, 2016
// *      Author: petzko
// */
//
#ifndef SRC_SIM_SIMLOOP_HPP_
#define SRC_SIM_SIMLOOP_HPP_
#include "SimSettings.hpp"
#include "SimData.hpp"
#include "mat.h"

void write2Matlab(mxArray *var, const char *varname,const char *filename);


void makeMaxwellVars(MB::SimSettings& set, MB::SimData& dat);
void makeBlochVars(MB::SimSettings& set, MB::SimData& dat);

void stepMaxwellVars(MB::SimSettings& set, MB::SimData& dat);
void stepBlochVars(MB::SimSettings& set, MB::SimData& dat);
void updateSolvers(MB::SimSettings& set, MB::SimData& dat);

void startSim(char* simFile,char* setFile);


#endif /* SRC_SIM_SIMLOOP_HPP_ */
