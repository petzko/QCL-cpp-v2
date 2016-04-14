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






// allocating memory for a matlab array:
// path to MATLAB LIBRARY -> do not forget to include the matlab so lib.
// LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/MATLAB/R2016a/bin/glnxa64/
// declaring & defining a matlab array: mxArray *e_field =  mxCreateDoubleMatrix(s.N_x, s.N_t/s.ct_plot, mxREAL);


void write2Matlab(mxArray *var, const char *varname,const char *filename){

	MATFile *pmat;
	pmat = matOpen(filename,"u");

	if(pmat == NULL){
		pmat = matOpen(filename, "w");
	}



	if (pmat == NULL) {
		printf("Error creating file %s\n", filename);
		return;
	}

	matPutVariable(pmat, varname, var);
	matClose(pmat);
}



void makeMaxwellVars(MB::SimSettings&  set, MB::SimData& dat){

	dat.U.init(1,set.N);
	// init u as a gaussian !
	double tp = 1;
	double x0 = 6*set.Ltot/7;
	double argument = 0.f ;
	for (int i = 0 ; i < set.N; i++ ){
		argument = ((dat.dx*i-x0)/dat.c)/tp;
		dat.U(0,i) = exp(-argument*argument);
	}

	dat.V.init(1,set.N);

	dat.U_solver = new MB::RNFDSolver<_TYPE_>(set.N,dat.dx,true,dat.c,dat.U);
	dat.V_solver = new MB::RNFDSolver<_TYPE_>(set.N,dat.dx,false,dat.c,dat.V);

	dat.losses = -dat.c*dat.l0;
	dat.K.init(1,set.N);
	dat.K = dat.losses*ones<_TYPE_>(1,set.N);

	dat.n32p_t.init(1,set.N);
	dat.n32m_t.init(1,set.N);

	dat.factor = -_i*dat.c*dat.trace_rho;


}
void makeBlochVars(MB::SimSettings& set, MB::SimData& dat){

	//  populations
	dat.r110.init(1,set.N); dat.r110 = 1./3*ones<_TYPE_>(1,set.N);
	dat.r330.init(1,set.N); dat.r330 = 1./3*ones<_TYPE_>(1,set.N);
	dat.r220.init(1,set.N); dat.r220 = 1./3*ones<_TYPE_>(1,set.N);
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

	dat.U = dat.U_solver->getLatestSolution();
	dat.V = dat.V_solver->getLatestSolution();



}
void stepMaxwellVars(MB::SimSettings&  set, MB::SimData& dat){

	_TYPE_ factor = dat.factor*dat.dipR;
	MB::Matrix<_TYPE_> tmp = factor*dat.n32p;
	MB::Matrix<_TYPE_> tmp_t = factor*dat.n32p_t;

	_TYPE_ bdry4V = dat.U_solver->makeStep(tmp,tmp_t,dat.K,dat.dt);
	_TYPE_ bdry4U = dat.V_solver->makeStep(tmp,tmp_t,dat.K,dat.dt);

	dat.U_solver->setBdry(bdry4U,0);
	dat.V_solver->setBdry(bdry4V,set.N-1);



}
void stepBlochVars(MB::SimSettings& set, MB::SimData& dat){

	size_t INJ = dat.INJ; size_t ULL = dat.ULL; size_t LLL = dat.LLL; size_t RES = dat.RES;
	size_t DEPOP = dat.DEPOP;





	MB::Matrix<_TYPE_> r110_t = _i*dat.HTB[INJ][ULL]*(dat.r130-conjVal(dat.r130)) +(dat.W[ULL][INJ] + dat.W[ULL][DEPOP])*dat.r330 + (dat.W[LLL][INJ]+dat.W[LLL][DEPOP])*dat.r220 + (dat.W[RES][INJ]+dat.W[RES][DEPOP])*dat.rRES- dat.G[INJ]*dat.r110;
	dat.r110_solver->makeStep(r110_t,dat.dt);


	MB::Matrix<_TYPE_>  lmInteraction = dat.dipR*(mMult(conjVal(dat.U),dat.n32p) + mMult(conjVal(dat.V),dat.n32m));

	MB::Matrix<_TYPE_> r330_t = _i*1.0*dat.HTB[INJ][ULL]*(conjVal(dat.r130) - dat.r130) +_i/2.0*(lmInteraction-conjVal(lmInteraction)) +  dat.r110*dat.W[INJ][ULL] + dat.r220*dat.W[LLL][ULL] + dat.rRES*dat.W[RES][ULL]- dat.G[ULL]*dat.r330;
	dat.r330_solver->makeStep(r330_t,dat.dt);

	MB::Matrix<_TYPE_> r220_t = -_i/2.0*(lmInteraction-conjVal(lmInteraction)) + dat.r110*dat.W[INJ][LLL]+ dat.r330*dat.W[ULL][LLL] +dat.rRES*dat.W[RES][LLL]- dat.G[LLL]*dat.r220;
	dat.r220_solver->makeStep(r220_t,dat.dt);

	MB::Matrix<_TYPE_> rRES_t =  dat.r110*dat.W[INJ][RES]+ dat.r330*dat.W[ULL][RES]+ dat.r220*dat.W[LLL][RES] - dat.G[RES]*dat.rRES;
	dat.rRES_solver->makeStep(rRES_t,dat.dt);

	// COHERENCES
	MB::Matrix<_TYPE_> r130_t = dat.dE13*dat.r130 + _i*dat.HTB[INJ][ULL]*(dat.r110-dat.r330) +_i/2.0*dat.dipR*(mMult(conjVal(dat.U),dat.n12p) + mMult(conjVal(dat.V),dat.n12m));
	dat.r130_solver->makeStep(r130_t,dat.dt);

	dat.n32p_t = dat.dE32*dat.n32p + _i/2.0*dat.dipR*(mMult(dat.U,(dat.r330-dat.r220)) + mMult(dat.V,(dat.r33p-dat.r22p))) - _i*dat.HTB[INJ][ULL]*dat.n12p;
	dat.n32p_solver->makeStep(dat.n32p_t,dat.dt);

	dat.n32m_t = dat.dE32*dat.n32m + _i/2.0*dat.dipR*(mMult(dat.V,(dat.r330-dat.r220)) + mMult(dat.U,conjVal(dat.r33p-dat.r22p))) - _i*dat.HTB[INJ][ULL]*dat.n12m;
	dat.n32m_solver->makeStep(dat.n32m_t,dat.dt);

	MB::Matrix<_TYPE_> n12p_t = dat.dE12*dat.n12p +_i/2.0*dat.dipR*(mMult(dat.U,dat.r130) + mMult(dat.V,dat.r13p)) - _i*dat.HTB[INJ][ULL]*dat.n32p;
	dat.n12p_solver->makeStep(n12p_t,dat.dt);

	MB::Matrix<_TYPE_> n12m_t = dat.dE12*dat.n12m + _i/2.0*dat.dipR*(mMult(dat.V,dat.r130) +mMult(dat.U,dat.r13m)) - _i*dat.HTB[INJ][ULL]*dat.n32m;
	dat.n12m_solver->makeStep(n12m_t,dat.dt);

	if(set.shb > 0 ){

		//r11+
		MB::Matrix<_TYPE_> r11p_t = _i*dat.HTB[INJ][ULL]*(dat.r13p-conjVal(dat.r13m)) + (dat.W[ULL][INJ]+dat.W[ULL][DEPOP])*dat.r33p+ (dat.W[LLL][INJ]+dat.W[LLL][DEPOP])*dat.r22p - (dat.G[INJ]+dat.diffusion)*dat.r11p;
		dat.r11p_solver->makeStep(r11p_t,dat.dt);

		/// r33+
		MB::Matrix<_TYPE_>  r33p_t = _i*dat.HTB[INJ][ULL]*(conjVal(dat.r13m)-dat.r13p)+_i/2.0*dat.dipR*(mMult(conjVal(dat.V),(dat.n32p)) -mMult((dat.U),conjVal(dat.n32m))) +  dat.W[INJ][ULL]*dat.r11p + dat.W[LLL][ULL]*dat.r22p - (dat.G[ULL]+dat.diffusion)*dat.r33p;
		dat.r33p_solver->makeStep(r33p_t,dat.dt);

		// r22+
		MB::Matrix<_TYPE_> r22p_t = -_i/2.0*dat.dipR*(mMult(conjVal(dat.V),(dat.n32p)) -mMult((dat.U),conjVal(dat.n32m))) + dat.W[INJ][LLL]*dat.r11p + dat.W[ULL][LLL]*dat.r33p - (dat.G[LLL]+dat.diffusion)*dat.r22p;
		dat.r22p_solver->makeStep(r22p_t,dat.dt);

		// r13+
		MB::Matrix<_TYPE_>  r13p_t = (dat.dE13-dat.diffusion)*dat.r13p + _i*dat.HTB[INJ][ULL]*(dat.r11p-dat.r33p) +_i/2.0*dat.dipR*(mMult(conjVal(dat.V),dat.n12p));
		dat.r13p_solver->makeStep(r13p_t,dat.dt);

		// r13-
		MB::Matrix<_TYPE_>  r13m_t = (dat.dE13-dat.diffusion)*dat.r13m + _i*dat.HTB[INJ][ULL]*conjVal(dat.r11p-dat.r33p) + _i/2.0*dat.dipR*(mMult(conjVal(dat.U),dat.n12m));
		dat.r13m_solver->makeStep(r13m_t,dat.dt);
	}













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

	dat.l0 = set.loss*100.0/(1/set.lch);
	std::cout << "losses: " << dat.l0  << "\n";

	dat.dx = set.Ltot/(set.N-1);
	dat.dt =  dat.dx/dat.c;
	std::cout << "timestep: " << dat.dt << "\n" ;
	std::cout << "grid size: " << dat.dx << "\n" ;

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
	dat.HTB.resize(dat.NLVL);

	for(size_t i = 0 ; i< dat.NLVL;i++){
		dat.W[i].resize(dat.NLVL);
		dat.HTB[i].resize(dat.NLVL);
		for(size_t j = 0 ; j< dat.NLVL;j++){
			dat.W[i][j] = set.Wmtx[j+i*dat.NLVL];
			dat.HTB[i][j] = set.HTB[j+i*dat.NLVL];
		}
	}

	std::cout<<"Wmtx \n";
	for(size_t i = 0 ; i< dat.NLVL;i++){
		for(size_t j = 0 ; j< dat.NLVL;j++)
			std::cout << dat.W[i][j]<<" ";
		std::cout << "\n";
	}

	std::cout<<"HTB \n";
	for(size_t i = 0 ; i< dat.NLVL;i++){
		for(size_t j = 0 ; j< dat.NLVL;j++)
			std::cout << dat.HTB[i][j]<<" ";
		std::cout << "\n";
	}

	dat.t = 0;
	dat.N_t = 3200;
	dat.iter_ctr = 0;
	makeMaxwellVars(set,dat);
	makeBlochVars(set,dat);

	mxArray *record_U = mxCreateDoubleMatrix(1,dat.N_t,mxCOMPLEX);
	mxArray *record_V = mxCreateDoubleMatrix(1,dat.N_t,mxCOMPLEX);




	int idx =  100;

	while (dat.iter_ctr < dat.N_t){

		/**
		 * store some data?
		 */

		stepBlochVars(set,dat);
		stepMaxwellVars(set,dat);
		updateSolvers(set, dat);

		*(mxGetPr(record_U)+dat.iter_ctr) = REAL(dat.U(0,idx));
		*(mxGetPi(record_U)+dat.iter_ctr) = IMAG(dat.U(0,idx));

		*(mxGetPr(record_V)+dat.iter_ctr) = REAL(dat.V(0,idx));
		*(mxGetPi(record_V)+dat.iter_ctr) = IMAG(dat.V(0,idx));



		if (dat.iter_ctr % 1000 == 0)
			std::cout << "Time @ iteration " << dat.iter_ctr << ": " << dat.t << "\n";

		dat.iter_ctr++ ;
		dat.t += dat.dt;
	}

	write2Matlab(record_U, "record_U","out_test.mat");
	write2Matlab(record_V, "record_V","out_test.mat");

	mxDestroyArray(record_U);
	mxDestroyArray(record_V);




}

