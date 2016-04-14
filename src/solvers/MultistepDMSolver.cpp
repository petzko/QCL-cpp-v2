/*
 * MultiStepDMSolver.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#include "MultistepDMSolver.hpp"


template<typename _Tp>
MB::MultistepDMSolver<_Tp>::MultistepDMSolver (unsigned int nrSteps,unsigned int nrPts,MB::Matrix<_Tp> initRhsDat,MB::Matrix<_Tp> initSol) : _MAXSTEPS(5) {

	assert(nrSteps>=1);
	assert(nrSteps <= _MAXSTEPS); // check if nrSteps is in order!
	unsigned int colDim = initRhsDat.getDim_j();
	assert((colDim >= 0 && colDim <= (nrSteps-1))); // check if init rhs data is in order
	assert(initSol.getDim_j() == nrPts); // check if vector dimensions are in order

	_m = nrSteps; _N = nrPts;
	_data.init(_m,_N);

	// assign the initial values in their correct location inside the _data member matrix
	for(int k = 1 ; k <= colDim; k++){
		_data.setSlice(k,k,0,_N-1);
		initRhsDat.setSlice(colDim-k,colDim-k,0,_N-1); // <- need to invert!
		_data = initRhsDat; // copy data!
	}
	_coefs = getCoeffs(colDim,_m);

	_sol.init(1,_N);
	_sol = initSol;
	_iterCtr = colDim;
}




template<typename _Tp>
MB::MultistepDMSolver<_Tp>::MultistepDMSolver (unsigned int nrSteps,unsigned int nrPts,MB::Matrix<_Tp> initSol) : _MAXSTEPS(5) {

	assert(nrSteps>=1);
	assert(nrSteps <= _MAXSTEPS); // check if nrSteps is in order!

	assert(initSol.getDim_j() == nrPts); // check if vector dimensions are in order

	_m = nrSteps; _N = nrPts;
	_data.init(_m,_N);

	_coefs = getCoeffs(1,_m);

	_sol.init(1,_N);
	_sol = initSol;
	_iterCtr = 1;
}

template<typename _Tp>
MB::Matrix<_Tp> MB::MultistepDMSolver<_Tp>::makeStep(MB::Matrix<_Tp> rhs, double dt){


	assert(rhs.getDim_j() == _N);
	assert(rhs.getDim_i() ==  1); // we expect only a row vector!

	_data.setSlice(0,0,0,_N-1);
	_data = rhs; // assign rhs to first position in data array
	_data.resetSlice();

	for (int k = 0; k < _m ; k++){
		_data.setSlice(k,k,0,_N-1);
//		_sol = _sol +dt*_coefs[k]*_data;
	}
	_data.resetSlice();

	_iterCtr++;
	if(_iterCtr<=_m)
		_coefs = getCoeffs(_iterCtr,_m);

	/*
	 *  circularly shift the rows of the data MTX
	 */
	for(int k = 0; k < _m-1; k++){
		_data.setRowPtr(k,_data.getMtxData()[k+1]);
	}
}

template<typename _Tp>
MB::Matrix<_Tp> MB::MultistepDMSolver<_Tp>::getLatestSolution(){
	return this->_sol;
}

template<typename _Tp>
void MB::MultistepDMSolver<_Tp>::setLatestSolution(MB::Matrix<_Tp> newsol){
	this->_sol = newsol;
}
template<typename _Tp>
MB::MultistepDMSolver<_Tp>::~MultistepDMSolver(){
	; // do nothing as long as we do not cause memory leaks! the destructors of the MB::Matrix<_Tp> fields WILL be called and memory on heap WILL be freed!
}

template<typename _Tp>
std::vector<double> MB::MultistepDMSolver<_Tp>::getCoeffs(unsigned int step,unsigned int mStep){
	if (step > this->_MAXSTEPS || step <= 0 ){
		MB_OUT_ERR("getCoeffs(..)! Maximum number of steps exceeds the given step argument or is (leq) than 0! ",__FILE__,__LINE__);
		throw std::domain_error("step param outside of the allowed domain!");
	}
	assert(step <= mStep);
	std::vector<double> cfs(mStep);// automatically initialized with zeros!

	switch (step){
	case 1:
		cfs[0] =  1.f;
		break;
	case 2:
		cfs[0] =3.f/2.f; cfs[1]= -1.f/2.f;
		break;
	case 3:
		cfs[0] = 23.f/12.f; cfs[1] =  -4.f/3.f ; cfs[2] = 5.f/12.f ;
		break;
	case 4:
		cfs[0] = 55.f/24.f;cfs[1] = -59.f/24.f; cfs[2] = 37.f/24.f; cfs[3] = -3.f/8.f;
		break;
	case 5:
		cfs[0] = 1901.f/720.f; cfs[1] =  -1387.f/360.f; cfs[2] = 109.f/30.f; cfs[3] =  -637.f/360.f; cfs[4] = 251.f/720.f;
		break;
	default:
		MB_OUT_ERR("getCoeffs(..)! How did you manage to sneak in a wrong value here?",__FILE__,__LINE__);
		throw std::domain_error("step param outside of the allowed domain!");
		break;
	}
	return cfs;

}

template class MB::MultistepDMSolver<float>;
template class MB::MultistepDMSolver< COMPLEXFLOAT >;
template class MB::MultistepDMSolver<double>;
template class MB::MultistepDMSolver< COMPLEXDOUBLE >;
