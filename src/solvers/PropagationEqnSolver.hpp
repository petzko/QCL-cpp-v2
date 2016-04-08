/*
 * PropagationEqnSolver.hpp
 *
 *  Created on: Apr 7, 2016
 *      Author: petzko
 */

#ifndef SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_
#define SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_
#include "../matrix/matrix.hpp"

namespace MB{

template<typename _Tp>

class PropagationEqnSolver{
public:
	virtual MB::Matrix<_Tp>& makeStep(MB::Matrix<_Tp> F, MB::Matrix<_Tp> F_t, MB::Matrix<_Tp> K, double dt) = 0;
};


}



#endif /* SRC_SOLVERS_PROPAGATIONEQNSOLVER_HPP_ */
