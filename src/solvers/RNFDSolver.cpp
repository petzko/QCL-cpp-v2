#include "RNFDSolver.hpp"


template<typename _Tp>
MB::RNFDSolver<_Tp>::RNFDSolver(unsigned int N, double dx, bool direction, double velocity, MB::Matrix<_Tp> U_0){

	if(N != U_0.getDim_i()*U_0.getDim_j()){
		MB_OUT_ERR("Input array size not equal to specified grid size N. Aborting! ",__FILE__,__LINE__);
		throw std::length_error("Length error exception thrown. Aborting! ");
	}

	this->_N = N;
	this->_dx = dx;
	this->_direction = direction;
	this->_c  = velocity;
	this->_sol.init(U_0.getDim_i(),U_0.getDim_j());
	this->_sol = U_0;
}
template<typename _Tp>
MB::Matrix<_Tp>& MB::RNFDSolver<_Tp>::makeStep(MB::Matrix<_Tp> F, MB::Matrix<_Tp> F_t, MB::Matrix<_Tp> K, double dt) {

	if(this->_c*dt != this->_dx){
		MB_OUT_ERR("RNFDsolver make step failed! Magic time step not selected!",__FILE__,__LINE__);
		throw std::domain_error("Magic time step requirement violated!");
	}

	if(F.getLen() != _N || F_t.getLen() != _N || K.getLen() != _N){
		MB_OUT_ERR("RNFDsolver makeStep!Input vector size mismatch! ",__FILE__,__LINE__);
		throw std::length_error("one of the input vectors, F, F_t or K's size does not correspond to the grid size! Aborting! ");

	}



	//
	//	 The solver assumes an equation of the form
	//	 D_t(E) = -/+ D_x(E) + f(x,t) + k E


	_sol.resetSlice();
	MB::Matrix<_Tp> E = _sol;

	double dt_2 = dt*dt/2;
	double dx_1 = 1/_dx;

	E.setSlice(0,0,1,_N-1);
	F.setSlice(0,0,1,_N-1);

	MB::Matrix<_Tp> E_x = E*dx_1;
	MB::Matrix<_Tp> F_x = F*dx_1;

	E.setSlice(0,0,0,_N-2);
	F.setSlice(0,0,2,_N-2);
	E_x = E_x-E*dx_1;
	F_x = F_x-F*dx_1;

	/** MATLAB CODE
	 *
	 *  if(obj.PROPERTIES(4) > 0)
	 *                   obj.SOLUTION(2:N) = E(1:end-1) +dt_2*( F_t(2:end) -c*F_x - 2*c*k(2:end).*E_x + ...
	 *                       k(2:end).*F(2:end)+(k(2:end).^2).*E(2:end)) +dt*(F(2:end)+k(2:end).*E(2:end));
	 *              else
	 *                   obj.SOLUTION(1:end-1) = E(2:end)+dt_2*(F_t(1:end-1) + c*F_x + 2*c*k(1:end-1).*E_x + ...
	 *                     k(1:end-1).*F(1:end-1)+(k(1:end-1).^2).*E(1:end-1)) +dt*(F(1:end-1)+k(1:end-1).*E(1:end-1));
	 *              end
	 *
	 */


	if(this->_direction > 0){
		// 1:N-1
		_sol.setSlice(0,0,1,_N-1);
		E.setSlice(0,0,0,_N-2);
		_sol = E;
		E.setSlice(0,0,1,_N-1); F_t.setSlice(0,0,1,_N-1); F.setSlice(0,0,1,_N-1);K.setSlice(0,0,1,_N-1);
		_sol = _sol+dt_2*(F_t-_c*F_x-2*_c*mMult(K,E_x) + mMult(K,F) + mMult(mMult(K,K),E)) + dt*(F+mMult(K,E));

	}else{
		//	obj.SOLUTION(0:N-2) = E(1:N-1) +dt_2*( F_t(0:N-2) + c*F_x + 2*c*k(0:N-2).*E_x + k(0:N-2).*F(0:N-2)+...
		//	(k(0:N-2).^2).*E(0:N-2)) +dt*(F(0:N-2)+k(0:N-2).*E(0:N-2));

		// 1:N-1
		_sol.setSlice(0,0,0,_N-2);
		E.setSlice(0,0,1,_N-1);
		_sol = E;
		E.setSlice(0,0,0,_N-2); F_t.setSlice(0,0,0,_N-2); F.setSlice(0,0,0,_N-2);K.setSlice(0,0,0,_N-2);
		_sol = _sol+dt_2*(F_t+_c*F_x+2*_c*mMult(K,E_x) + mMult(K,F) + mMult(mMult(K,K),E)) + dt*(F+mMult(K,E));
	}

	this->_sol.resetSlice();
	return this->_sol;

}

template<typename _Tp>
void  MB::RNFDSolver<_Tp>::setBdry(_Tp val,unsigned int flag){
	if (flag >1 ){
		MB_OUT_ERR(" setBdry(:)! Boundary location flag outside allowed domain.",__FILE__,__LINE__);
		throw std::domain_error("Check boundary flag? ") ;

	}
	this->_sol(0,flag) = val;
}

template<typename _Tp>
void  MB::RNFDSolver<_Tp>::setLatestSolution(MB::Matrix<_Tp> sol){
	if(this->_sol.getDim_i() != sol.getDim_i() || this->_sol.getDim_j() != sol.getDim_j()){
		MB_OUT_ERR("Input array size not equal to previously stored array size. Aborting! ",__FILE__,__LINE__);
		throw std::length_error("Length error exception thrown. Aborting! ");
	}else{
		this->_sol = sol;
	}
}

template<typename _Tp>
MB::Matrix<_Tp>&  MB::RNFDSolver<_Tp>::getLatestSolution(){
	return this->_sol;
}
template<typename _Tp>
MB::RNFDSolver<_Tp>::~RNFDSolver(){
	;
}


template class MB::RNFDSolver<float>;
template class MB::RNFDSolver< complex float>;
template class MB::RNFDSolver<double>;
template class MB::RNFDSolver<complex double>;
