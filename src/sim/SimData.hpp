///*
// * SimData.hpp
// *
// *  Created on: Apr 8, 2016
// *      Author: petzko
// */
//
#ifndef SRC_SIM_SIMDATA_HPP_
#define SRC_SIM_SIMDATA_HPP_

#include <string>
#include "../matrix/matrix.hpp"
#include "../common/utils.hpp"

#include "../solvers/PropagationEqnSolver.hpp"
#include "../solvers/DMSolver.hpp"



namespace MB{

struct SimData{
	double c,T_R,f_R,hbar;
	size_t INJ,ULL,LLL,RES,DEPOP;
	size_t NLVL;

	double E0,l0;
	double dx,dt;

	double diffusion;
	double zUL, Ncarriers,trace_rho;
	double t;
	size_t iter_ctr;
	size_t N_t;

	std::vector< std::vector<double> > W;
	std::vector< std::vector<double> > HTB;
	std::vector<double>  G;

	double losses;
	double dipR = 1;
	COMPLEXDOUBLE factor;

	COMPLEXDOUBLE dE32;
	COMPLEXDOUBLE dE13;
	COMPLEXDOUBLE dE12;



	// Maxwell Variables
	MB::Matrix<_TYPE_> U;
	MB::Matrix<_TYPE_> V;

	//wave eqn. solvers!
	MB::PropagationEqnSolver<_TYPE_>* U_solver;
	MB::PropagationEqnSolver<_TYPE_>* V_solver;


	// Bloch Variables

	//populations
	MB::Matrix<_TYPE_> r110;
	MB::Matrix<_TYPE_> r330;
	MB::Matrix<_TYPE_> r220;
	MB::Matrix<_TYPE_> rRES;

	//SHB POPS
	MB::Matrix<_TYPE_> r11p;
	MB::Matrix<_TYPE_> r33p;
	MB::Matrix<_TYPE_> r22p;

	//COHERENCES
	MB::Matrix<_TYPE_> r130;
	MB::Matrix<_TYPE_> n32p;
	MB::Matrix<_TYPE_> n32p_t;

	MB::Matrix<_TYPE_> n32m;
	MB::Matrix<_TYPE_> n32m_t;

	MB::Matrix<_TYPE_> K;


	MB::Matrix<_TYPE_> n12p;
	MB::Matrix<_TYPE_> n12m;

	// coherences SHB
	MB::Matrix<_TYPE_> r13p;
	MB::Matrix<_TYPE_> r13m;

	////	 bloch solvers!
	MB::DMSolver<_TYPE_>* r110_solver = NULL;
	MB::DMSolver<_TYPE_>* r330_solver = NULL;
	MB::DMSolver<_TYPE_>* r220_solver = NULL;
	MB::DMSolver<_TYPE_>* rRES_solver = NULL;

	MB::DMSolver<_TYPE_>* r11p_solver = NULL;
	MB::DMSolver<_TYPE_>* r33p_solver = NULL;
	MB::DMSolver<_TYPE_>* r22p_solver = NULL;

	MB::DMSolver<_TYPE_>* r130_solver = NULL;
	MB::DMSolver<_TYPE_>* n32p_solver = NULL;
	MB::DMSolver<_TYPE_>* n32m_solver = NULL;
	MB::DMSolver<_TYPE_>* n12p_solver = NULL;
	MB::DMSolver<_TYPE_>* n12m_solver = NULL;

	MB::DMSolver<_TYPE_>* r13p_solver = NULL;
	MB::DMSolver<_TYPE_>* r13m_solver = NULL;


	~SimData(){
		if(r110_solver != NULL)
			delete r110_solver;
		if(r330_solver != NULL)
			delete r330_solver;
		if(r220_solver != NULL)
			delete r220_solver;
		if(rRES_solver != NULL)
			delete rRES_solver;

		if(r11p_solver != NULL)
			delete r11p_solver;
		if(r33p_solver != NULL)
			delete r33p_solver;
		if(r22p_solver != NULL)
			delete r22p_solver;

		if(r130_solver != NULL)
			delete r130_solver;
		if(n32p_solver != NULL)
			delete n32p_solver;
		if(n32m_solver != NULL)
			delete n32m_solver;
		if(n12p_solver != NULL)
			delete n12p_solver;
		if(n12m_solver != NULL)
			delete n12m_solver;
		if(r13p_solver != NULL)
			delete r13p_solver;
		if(r13m_solver != NULL)
			delete r13m_solver;
		if(U_solver != NULL)
			delete U_solver;
		if(V_solver != NULL)
			delete V_solver;


	}
};


}
#endif /* SRC_SIM_SIMDATA_HPP_ */
