/*
 * DMSolver.hpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#ifndef SRC_SOLVERS_DMSOLVER_HPP_
#define SRC_SOLVERS_DMSOLVER_HPP_
#include "../matrix/matrix.hpp"
#include "../common/utils.hpp"

namespace MB{

template<typename _Tp>

class DMSolver{

	virtual MB::Matrix<_Tp>& makeStep(MB::Matrix<_Tp> rhs, double dt) = 0;

};

}
#endif /* SRC_SOLVERS_DMSOLVER_HPP_ */
