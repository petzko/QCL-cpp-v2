/*
 * PropagationEqnSolver.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: petzko
 */

#ifndef SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_
#define SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_
#include "../matrix/matrix.hpp"
#include "../common/utils.hpp"

namespace MB{

template<typename _Tp>
class PropagationEqnSolver{

public:
	virtual _Tp makeStep(MB::Matrix<_Tp>& F, MB::Matrix<_Tp>& F_t, MB::Matrix<_Tp>& K, double dt) = 0;

	virtual void setBdry(_Tp val,unsigned int flag)= 0;
	virtual	void setLatestSolution(Matrix<_Tp> sol)= 0;
	virtual	Matrix<_Tp> getLatestSolution()= 0;
	virtual ~PropagationEqnSolver(){;}

};


}



#endif /* SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_ */
