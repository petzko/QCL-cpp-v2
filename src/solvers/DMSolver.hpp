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

public:
	virtual MB::Matrix<_Tp> makeStep(MB::Matrix<_Tp> rhs, double dt) = 0;
	virtual MB::Matrix<_Tp> getLatestSolution() = 0;
	virtual void setLatestSolution(MB::Matrix<_Tp> newsol)= 0;
	virtual ~DMSolver(){;} ;
};

}
#endif /* SRC_SOLVERS_DMSOLVER_HPP_ */
