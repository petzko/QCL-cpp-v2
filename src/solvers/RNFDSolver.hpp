/*
 * RNFDSolver.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: petzko
 */

#ifndef SRC_SOLVERS_RNFDSOLVER_HPP_
#define SRC_SOLVERS_RNFDSOLVER_HPP_

#include "PropagationEqnSolver.hpp"


namespace MB{
template<typename _Tp>
class RNFDSolver: public PropagationEqnSolver<_Tp>{

private:
	int _N;
	double _dx;
	bool p_direction;
	double _c;
	MB::Matrix<_Tp> _sol;

public:

	RNFDSolver(unsigned int N, double dx, bool direction, double velocity, MB::Matrix<_Tp> U_0);
	_Tp makeStep(MB::Matrix<_Tp>& F, MB::Matrix<_Tp>& F_t, MB::Matrix<_Tp>& K, double dt) ;

	/**
	 * set the boundary of the wave vector!
	 * @param val -> of type _Tp is the value to be set at the boundary.
	 * @param flag -> of type int is the index of the element at which the boundary is to be set
	 */
	void setBdry(_Tp val,unsigned int flag);
	void setLatestSolution(Matrix<_Tp> sol);
	Matrix<_Tp> getLatestSolution();

	virtual ~RNFDSolver();



};

}


#endif /* SRC_SOLVERS_RNFDSOLVER_HPP_ */
