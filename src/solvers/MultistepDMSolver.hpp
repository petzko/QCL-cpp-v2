/*
 * MultistepDMSolver.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#ifndef SRC_SOLVERS_MULTISTEPDMSOLVER_HPP_
#define SRC_SOLVERS_MULTISTEPDMSOLVER_HPP_

#include "DMSolver.hpp"
namespace MB{

template<typename _Tp>
class MultistepDMSolver: public DMSolver<_Tp>{

private:

	unsigned int const _MAXSTEPS;
	// coeffs are stored in the following manner " newest  <--->  oldest  <--->  zeros padding"
	std::vector<double> _coefs; // = [ 1 0 0 0 0 ;  3/2 -1/2 0 0 0; 23/12 -4/3  5/12  0  0 ;55/24 -59/24  37/24 -3/8  0  ; 1901/720  -1387/360  109/30 -637/360  251/720];
	unsigned int _m; // nr of steps...
	unsigned int _N; // nr of points on the spatial grid..
	//data container storing the rhs of the previous "m" steps;
	 // data is a R^(Nxm) matrix ...
	MB::Matrix<_Tp> _data;
	// a vector of size N which stores the solution at the previous time step
	MB::Matrix<_Tp> _sol;
	//counts how many iterations have been made so far...
	unsigned int  _iterCtr;
	unsigned int  _storageIdx;
	std::vector<double> getCoeffs(unsigned int step,unsigned int mStep);

public:

	MultistepDMSolver(unsigned int nrSteps,unsigned int nrPts,MB::Matrix<_Tp> initRhsDat,MB::Matrix<_Tp> initSol);
	MB::Matrix<_Tp> makeStep(MB::Matrix<_Tp> rhs, double dt);
	MB::Matrix<_Tp> getLatestSolution();

	void setLatestSolution(MB::Matrix<_Tp> newsol);
	virtual ~MultistepDMSolver();

};
}

#endif /* SRC_SOLVERS_MULTISTEPDMSOLVER_HPP_ */
