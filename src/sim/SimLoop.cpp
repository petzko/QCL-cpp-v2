/*
 * SimLoop.cpp
 *
 *  Created on: Apr 11, 2016
 *      Author: petzko
 */



#include "SimLoop.hpp"
#include "SimSettings.hpp"
#include "SimData.hpp"

#include "../solvers/RNFDSolver.hpp"
#include "../solvers/MultistepDMSolver.hpp"

#include "../common/utils.hpp"
#include "../common/CONSTANTS.hpp"

void makeMaxwellVars(MB::SimSettings&  set, MB::SimData& dat){
	dat.U.init(1,set.N);
	dat.V.init(1,set.N);

	dat.U_solver = new MB::RNFDSolver<_TYPE_>(set.N,dat.dx,1,dat.c,dat.U);
	dat.V_solver = new MB::RNFDSolver<_TYPE_>(set.N,dat.dx,-1,dat.c,dat.V);

	dat.losses = -dat.c*dat.l0;
	dat.K.init(1,set.N);
	dat.K = dat.losses*dat.K;
	dat.n32p_t.init(1,set.N);
	dat.n32m_t.init(1,set.N);

	dat.factor = -std::complex<double>(0,1)*dat.c*dat.trace_rho;


}
void makeBlochVars(MB::SimSettings& set, MB::SimData& dat){

	//  populations
	dat.r110.init(1,set.N); dat.r110 = 1./3.*ones<_TYPE_>(1,set.N);
	dat.r330.init(1,set.N); dat.r330 = 1./3.*ones<_TYPE_>(1,set.N);
	dat.r220.init(1,set.N); dat.r220 = 1./3.*ones<_TYPE_>(1,set.N);
	dat.rRES.init(1,set.N);

	// coherences
	dat.r130.init(1,set.N); dat.r130 = 1e-15*(randInit<_TYPE_>(1,set.N)+_TYPE_(0,1)*randInit<_TYPE_>(1,set.N,time(NULL)));

	dat.n32p.init(1,set.N); dat.n32p = 1e-15*(randInit<_TYPE_>(1,set.N)+_TYPE_(0,1)*randInit<_TYPE_>(1,set.N,time(NULL)));
	dat.n32m.init(1,set.N); dat.n32m = 1e-15*(randInit<_TYPE_>(1,set.N)+_TYPE_(0,1)*randInit<_TYPE_>(1,set.N,time(NULL)));

	dat.n12p.init(1,set.N);
	dat.n12m.init(1,set.N);

	// solvers
	dat.r110_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r110);
	dat.r330_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r330);
	dat.r220_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r220);
	dat.rRES_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.rRES);

	dat.r130_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r130);
	dat.n32p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.n32p);
	dat.n32m_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.n32m);
	dat.n12p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.n12p);
	dat.n12m_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.n12m);


	if(set.shb > 0){
		dat.r11p.init(1,set.N);
		dat.r11p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r11p);
		dat.r33p.init(1,set.N);
		dat.r33p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r33p);
		dat.r22p.init(1,set.N);
		dat.r22p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r22p);
		dat.r13p.init(1,set.N);
		dat.r13p_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r13p);
		dat.r13m.init(1,set.N);
		dat.r13m_solver = new MB::MultistepDMSolver<_TYPE_>(set.nrSteps,set.N,dat.r13m);
	}



}

void updateSolvers(MB::SimSettings& set, MB::SimData& dat){

	dat.r110 = dat.r110_solver->getLatestSolution();
	dat.r330 = dat.r330_solver->getLatestSolution();
	dat.r220 = dat.r220_solver->getLatestSolution();
	dat.rRES = dat.rRES_solver->getLatestSolution();


	dat.r130 = dat.r130_solver->getLatestSolution();

	dat.n32p = dat.n32p_solver->getLatestSolution();
	dat.n32m = dat.n32m_solver->getLatestSolution();

	dat.n12p = dat.n12p_solver->getLatestSolution();
	dat.n12m = dat.n12m_solver->getLatestSolution();


	if(set.shb > 0 ){

		dat.r11p = dat.r11p_solver->getLatestSolution();
		dat.r33p = dat.r33p_solver->getLatestSolution();
		dat.r22p = dat.r22p_solver->getLatestSolution();

		dat.r13p = dat.r13p_solver->getLatestSolution();
		dat.r13m = dat.r13m_solver->getLatestSolution();
	}



}
void stepMaxwellVars(MB::SimSettings&  set, MB::SimData& dat){

	_TYPE_ factor = dat.factor*dat.dipR;

	_TYPE_ bdry4V = dat.U_solver->makeStep(factor*dat.n32p,factor*dat.n32p_t,dat.K,dat.dt);
	_TYPE_ bdry4U = dat.V_solver->makeStep(factor*dat.n32m,factor*dat.n32m_t,dat.K,dat.dt);

	dat.U_solver->setBdry(bdry4U,0);
	dat.V_solver->setBdry(bdry4V,set.N-1);
}
void stepBlochVars(MB::SimSettings& set, MB::SimData& dat){
	/**
	 * TODO major stuff!
	 *
	 *
	 */
}





void startSim(char* simFile,char* setFile){

	MB::SimSettings set(simFile,setFile);
	set.initSimSettings();
	MB::SimData dat;

	dat.c = C0*1./(set.lch/set.tch)/set.nTHz;
	std::cout <<"phase velocity: " << dat.c << "\n";

	dat.T_R = 2*set.Ltot/dat.c; dat.f_R = 1/dat.T_R;

	dat.hbar = HBAR/Q0/set.tch;

	dat.INJ = 0; dat.ULL = 1;dat.LLL = 2; dat.RES = 3; dat.DEPOP = 4;
	dat.NLVL = 5;

	dat.E0 = (set.HTB[dat.ULL+dat.NLVL*dat.ULL]-set.HTB[dat.LLL+dat.NLVL*dat.LLL])/dat.hbar;
	std::cout <<"Central Frequency (THz): " << dat.E0/2./M_PI << "\n";

	dat.l0 = set.loss*100/(1/set.lch);

	dat.dx = set.Ltot/(set.N-1);
	dat.dt =  dat.dx/dat.c;
	dat.diffusion =4*dat.E0/(dat.c*dat.c)*set.D*1e2/(1/set.tch);
	dat.zUL = set.zUL;
	dat.dipR = 1;

	double zUL2 = (dat.zUL*1e-9*Q0)*(dat.zUL*1e-9*Q0);
	dat.Ncarriers = set.dN*1e3*set.Ld/set.Lp;
	std::cout << "Carrier density! " << dat.Ncarriers << "\n";
	dat.trace_rho = (dat.E0*1e12*dat.Ncarriers*set.Overlap*zUL2)/(EPS0*set.nTHz*C0*HBAR)/(1./set.tch/set.lch);
	std::cout << "Trace: " << dat.trace_rho << "\n";

	dat.G.resize(dat.NLVL);
	dat.W.resize(dat.NLVL);
	for(int i = 0 ; i< dat.NLVL;i++)
		dat.W[i].resize(dat.NLVL);

	for(int i = 0 ; i< dat.NLVL;i++)
		for(int j = 0 ; j< dat.NLVL;j++)
			dat.W[i][j] = set.Wmtx[j+i*dat.NLVL];

	for(int i = 0 ; i< dat.NLVL;i++){
		for(int j = 0 ; j< dat.NLVL;j++)
			std::cout << dat.W[i][j]<<" ";
		std::cout << "\n";
	}

	dat.t = dat.dt;
	dat.N_t = int(set.simRT*dat.T_R);
	makeMaxwellVars(set,dat);
	makeBlochVars(set,dat);

//	while (dat.iter_ctr < dat.N_t){
//
//		/**
//		 * store some data?
//		 *
//		 */
//
//		stepBlochVars(set,dat);
//		stepMaxwellVars(set,dat);
//
//		dat.iter_ctr ++ ;
//		dat.t += dat.dt;
//	}

}

