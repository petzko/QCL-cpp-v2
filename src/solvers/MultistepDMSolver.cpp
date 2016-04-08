/*
 * MultiStepDMSolver.cpp
 *
 *  Created on: Apr 8, 2016
 *      Author: petzko
 */

#include "MultistepDMSolver.hpp"


template<typename _Tp>
MB::MultistepDMSolver<_Tp>::MultistepDMSolver (unsigned int nrSteps,unsigned int nrPts,MB::Matrix<_Tp> initRhsDat,MB::Matrix<_Tp> initSol) : _MAXSTEPS(5) {


	/**
	 *    MATLAB CODE !
	 *    % do some error checking ...
	 *           assert(nr_steps >=1,'please set a positive nr of steps for the MULTIstep algorithm to run correctly.');
	 *           assert(nr_steps <= 5,['unfortunately current implementation does not support more then a 5 step adams bashforth algorithm!'...
	 *               ' Please specify a new nr of steps parameter and try again,']);
	 *           [ rowdim,coldim ]  = size(init_rhs_data);
	 *           % do some assertions
	 *           assert((coldim >=0 && coldim <=(nr_steps-1)), 'Incorrectly specified initial data! Please try again');
	 *
	 *           %if all data is correctly initalized, proceed with the
	 *           %constructor...
	 *
	 *           obj.m = nr_steps;
	 *           obj.N = nr_pts;
	 *
	 *           obj.data = zeros(obj.N,obj.m);
	 *           for k = 1:coldim
	 *               obj.data(:,k+1) = init_rhs_data(:,coldim - k+1);
	 *           end
	 *
	 *          % solution at step n2
	 *           obj.coefs = get_coeffs(obj.m);
	 *          obj.prev_solution = init_solution(); % this sets - up the inital data...
	 *           obj.iter_ctr = coldim+1;
	 */

	assert(nrSteps>=1);
	assert(nrSteps <= _MAXSTEPS); // check if nrSteps is in order!
	unsigned int rowDim = initRhsDat.getDim_i();
	unsigned int colDim = initRhsDat.getDim_j();
	assert((colDim >= 0 && colDim <= (nrSteps-1))); // check if init rhs data is in order
	assert(rowDim == nrPts); // check if vector dimensions are in order

	_m = nrSteps; _N = nrPts;
	_data = MB::Matrix<_Tp>(_m,_N);

	// assign the initial values in their correct location inside the _data member matrix
	for(int k = 1 ; k <= colDim; k++){
		_data.setSlice(k,k,0,_N-1);
		initRhsDat.setSlice(colDim-k,colDim-k,0,_N-1); // <- need to invert!
		_data = initRhsDat; // copy data!
	}
	_coefs = getCoeffs(colDim,_m);

	_sol = initSol;
	_iterCtr = colDim;
}

template<typename _Tp>
MB::Matrix<_Tp> MB::MultistepDMSolver<_Tp>::makeStep(MB::Matrix<_Tp> rhs, double dt){
	/**
	 *		MATLAB CODE:
	 *		    % Make a single propagation step of size "dt"  from tn = n*dt -> to tn+1 = (n+1)*dt, using the stored
	 *            % previous data. input arguments are the solver object itself,
	 *           % the right hand side of the equation evaluated at the current timestep (tn) and the timestep size dt
	 *
	 *           [n1,n2] = size(prev_rhs);
	 *           %check if vector dimensions are consice..
	 *           err_msg = ['Cannot evolve equation ' ...
	 *               'Rhs vector dimension does not agree with solution vector dimension.' ];
	 *           assert(n1 == obj.N,err_msg);
	 */
	assert(rhs.getDim_j() == _N);
	assert(rhs.getDim_i() ==  1); // we expect only a row vector!

	/*           %%% this ensures that we get the right initial conditions
	 *
	 *           % replace the oldest rhs with the newest rhs
	 *           obj.data(:,1) = prev_rhs;  % finally store the result...
	 */
	_data.setSlice(0,0,0,_N-1);
	_data = rhs; // assign rhs to first position in data array
	_data.resetSlice();

	for (int k = 0; k < _m ; k++){
		_data.setSlice(k,k,0,_N-1);
		_sol = _sol +dt*_coefs[k]*_data;
	}
	_data.resetSlice();

	_iterCtr++;
	if(_iterCtr<=_m)
		_coefs = getCoeffs(_iterCtr,_m);

	/*
	 * circularly shift the rows of the data MTX
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
	this->_prefSolution = newsol;
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
