/*
 * matrix.cpp
 *
 *  Created on: Feb 4, 2015
 *      Author: petz
 */

#include "matrix/matrix.hpp"

using namespace MB;


template<typename _Tp>
int Matrix<_Tp>::getDim_i() const {
	return _dimI;
}

template<typename _Tp>
int Matrix<_Tp>::getDim_j() const {
	return _dimJ;
}

template<typename _Tp>
_Tp** const Matrix<_Tp>::getMtxData() const {
	return _data;
}

/**
 * the standard constructor - takes a row and column dimension.
 * Initializes the matrix data in a row-major format (i.e. data in a single row is ordered contiguously in memory)
 * the param dim_i specifies the number of rows and the param dim_j the number of columns.
 */
template<typename _Tp>
_Tp* Matrix<_Tp>::setRowPtr(int fromIdx,_Tp* toPtr){

	_Tp* fromPtr = _data[fromIdx];
	_data[fromIdx] = toPtr;
	return fromPtr;
}


template<typename _Tp>
void Matrix<_Tp>::init(int dim_i, int dim_j) {

	_dimI = dim_i;
	_dimJ = dim_j;

	/**
	 * FULLY CONTIGUOUS ALLOCATION OF DATA IN MEMORY does not work
	 * with swaprows functionality!
	 *
	 * evaluate if the speedup I get from swappig pointers outweights the speedup which I get from having all data
	 * stored contiguously in memory!
	 *
	 *  _Tp* tmp_data = (_Tp*) calloc(_dimI * _dimJ, sizeof(_Tp));
	 *	_data = (_Tp**) calloc(_dimI, sizeof(_Tp*));
	 *	for (int i = 0; i < _dimI; i++)
	 *		_data[i] = tmp_data + i * _dimJ;
	 *	_Tp* tmp_data = (_Tp*) calloc(_dimI * _dimJ, sizeof(_Tp));
	 *
	 *
	 */


	_data = (_Tp**) calloc(_dimI, sizeof(_Tp*));

	for (int i = 0; i < _dimI; i++)
		_data[i] =  (_Tp*) calloc(_dimJ, sizeof(_Tp));;

	_bytes = _dimI * _dimJ * sizeof(_Tp);
	_i1 =  0; _i2 =  _dimI-1;
	_j1 =  0; _j2 = _dimJ-1 ;

	_len = _dimI*_dimJ;
}


template<typename _Tp>
Matrix<_Tp>::Matrix(){
	_data = NULL;
}


/*
 * Destructor
 */
template<typename _Tp>
Matrix<_Tp>::~Matrix() {
	// first free the data elements
	//	free(_data[0]);

	if(_data != NULL){
		for (int i = 0; i < _dimI; i++){
			if(_data[i] !=NULL)
				free(_data[i]);
		}
		free(_data);

	}
}

template<typename _Tp>
Matrix<_Tp>::Matrix(int dim_i, int dim_j) {
	init(dim_i, dim_j);
	_i1 = 0;  _i2 =   _dimI-1 ;
	_j1 = 0; _j2 = _dimJ-1;
}

template<typename _Tp >
Matrix<_Tp> Matrix<_Tp>::getSliceMtx(int i1, int i2, int j1,
		int j2){


	int rows = i2-i1+1;
	int cols = j2-j1+1;
	if(rows < 0 || rows > this->_dimI || cols < 0 || cols > this->_dimJ)
	{
		MB_OUT_ERR("Submatrix is of invalid size. Cannot proceed",
				__FILE__, __LINE__);
		throw std::length_error("Length error exception thrown");
	}

	//	copy row by row
	Matrix<_Tp> res(rows,cols);
	switch (gettype<_Tp>()) {
	case FLT:

		// copy the r+i1-th row of the current matrix into the rth row of the output.
		// the initial copy shall start from column index j1.
		for (int r = 0; r < rows; r++)
			cblas_scopy(cols, (float*)_data[r+i1]+j1, 1,
					(float*) res.getMtxData()[r], 1);
		break;

	case DBL:
		for (int r = 0; r < rows; r++)
			cblas_dcopy(cols, (double*)_data[r+i1]+j1, 1,
					(double*) res.getMtxData()[r], 1);
		break;
	case CPLX_FLOAT:
		for(int r = 0 ; r < rows ;r ++ )
			cblas_ccopy(cols, (void*)_data[r+i1]+j1, 1,
					(void*) res.getMtxData()[r], 1);
		break;
	case CPLX_DOUBLE:
		for(int r = 0; r < rows;r++ )
			cblas_zcopy(cols,(void*)_data[r+i1]+j1, 1,
					(void*) res.getMtxData()[r], 1);
		break;

	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}

	return res;

}
template<typename _Tp>
std::vector<int> MB::Matrix<_Tp>::getSliceVector() const{
	std::vector<int> res(4);
	res[0] = _i1; res[1] = _i2; res[2] = _j1; res[3] = _j2;
	return res;
}

template<typename _Tp>
void MB::Matrix<_Tp>::setSlice(int i1, int i2, int j1, int j2){


	if(i1 < 0 || i2 >= this->_dimI || j1 < 0 || j2 >= this->_dimJ){
		MB_OUT_ERR(
				"setSlice(...) ! Wrongly specified submatrix indices - dimension constraint violation",
				__FILE__, __LINE__);
		throw std::length_error("Length error exception thrown. Aborting! ");
	}
	this->_i1 = i1; this->_i2 = i2;
	this-> _j1 = j1; this->_j2 = j2;
}

template<typename _Tp>
void MB::Matrix<_Tp>::resetSlice(){
	this->setSlice(0,this->_dimI-1,0,this->_dimJ-1);
}


template<typename _Tp >
void Matrix<_Tp>::operator()(Matrix<_Tp> inMatrix, int i1, int i2, int j1,
		int j2){


	int rows = i2-i1+1;
	int cols = j2-j1+1;

	if(rows < 0 || rows > this->_dimI || cols < 0 || cols > this->_dimJ)
	{
		MB_OUT_ERR(	"operator(...)! Wrongly specified submatrix indices - dimension constraint violation",
				__FILE__, __LINE__);
		throw std::length_error("Length error exception thrown. Aborting! ");
	}
	if(rows != inMatrix.getDim_i() || cols != inMatrix.getDim_j())
	{
		MB_OUT_ERR("Input matrix dimensions do not match the requested array dimensions",
				__FILE__, __LINE__);
		throw std::length_error("Length error exception thrown. Aborting! ");
	}



	switch (gettype<_Tp>()) {
	case FLT:

		// copy the r+i1-th row of the current matrix into the rth row of the output.
		// the initial copy shall start from column index j1.
		for (int r = 0; r < rows; r++)
			cblas_scopy(cols, (float*)inMatrix.getMtxData()[r], 1,
					(float*) this->getMtxData()[r+i1]+j1, 1);
		break;

	case DBL:
		for (int r = 0; r < rows; r++)
			cblas_dcopy(cols, (double*)inMatrix.getMtxData()[r], 1,
					(double*) this->getMtxData()[r+i1]+j1, 1);
		break;
	case CPLX_FLOAT:
		for(int r = 0 ; r < rows ;r ++ )
			cblas_ccopy(cols, (void*)inMatrix.getMtxData()[r], 1,
					(void*) this->getMtxData()[r+i1]+j1, 1);
		break;
	case CPLX_DOUBLE:
		for(int r = 0; r < rows;r++ )
			cblas_zcopy(cols,(void*)inMatrix.getMtxData()[r], 1,
					(void*) this->getMtxData()[r+i1]+j1, 1);
		break;

	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}


}

/**
 *	Copy constructor
 */
template<typename _Tp>
Matrix<_Tp>::Matrix(const Matrix<_Tp>& arg) {

	init(arg._i2-arg._i1+1, arg._j2-arg._j1+1);

	int N = arg._j2-arg._j1+1;


	switch (gettype<_Tp>()) {
	case FLT:
		for (int i = arg._i1; i <= arg._i2;i++)
			cblas_scopy(N, (float*) arg.getMtxData()[i]+arg._j1, 1,	(float*) this->getMtxData()[i-arg._i1], 1);
		break;

	case DBL:
		for (int i = arg._i1; i <= arg._i2;i++)
			cblas_dcopy(N, (double*) (arg.getMtxData()[i]+arg._j1), 1,(double*) (this->getMtxData()[i-arg._i1]), 1);
		break;
	case CPLX_FLOAT:
		for (int i = arg._i1; i <= arg._i2;i++)
			cblas_ccopy(N, (void*) ((COMPLEXFLOAT*)arg.getMtxData()[i]+arg._j1), 1,(void*) (this->getMtxData()[i-arg._i1]), 1);
		break;
	case CPLX_DOUBLE:
		for (int i = arg._i1; i <= arg._i2;i++)
			cblas_zcopy(N, (void*) ((COMPLEXDOUBLE*)(arg.getMtxData()[i]+arg._j1)), 1,(void*) (this->getMtxData()[i-arg._i1]), 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");

	}
}


// make the function call operator retrive the i,j th data element
template<typename _Tp>
_Tp& Matrix<_Tp>::operator()(int const i,   int const j) const {

	if (i >= this->getDim_i() || i < 0) {
		MB_OUT_ERR(" operator(i,j)! Array index out of bounds. Cannot retrieve Matrix element. ", 	__FILE__, __LINE__);
		throw std::out_of_range("row index out of bounds.");

	}

	if (j >= this->getDim_j() || j < 0) {
		MB_OUT_ERR( "operator(i,j)! Array index out of bounds. Cannot retrieve Matrix element. ",__FILE__, __LINE__);
		throw std::out_of_range("column index out of bounds.");
	}

	return _data[i][j];

}
/* 
 * overload the assignment operator!
 *
 */
template<typename _Tp>
Matrix<_Tp>& Matrix<_Tp>::operator=(const Matrix<_Tp>& arg) {

	if(this == &arg)
		return *this;

	if (_i2-_i1 != arg._i2-arg._i1 	|| _j2-_j1 != arg._j2-arg._j1) {
		MB_OUT_ERR(" copy assignment operator! Mtx dimensions mismatch!",__FILE__,__LINE__);
		throw std::length_error("Matrix Dimensions do not aggree.");
	}
	std::vector<int> slice = arg.getSliceVector();

	int N = _j2-_j1+1;

	switch (gettype<_Tp>()){
	case FLT:
		for (int i = _i1; i<= _i2;i++)
			cblas_scopy(N, (float*) arg.getMtxData()[i+(arg._i1-_i1)] + arg._j1, 1,(float*) this->getMtxData()[i]+ _j1, 1);
		break;

	case DBL:
		for (int i = _i1; i<= _i2;i++)
			cblas_dcopy(N, (double*) arg.getMtxData()[i+(arg._i1-_i1)]+arg._j1, 1, (double*) this->getMtxData()[i]+_j1, 1);
		break;
	case CPLX_FLOAT:
		for (int i = _i1; i<= _i2;i++)
			cblas_ccopy(N, (void*) ((COMPLEXFLOAT*)arg.getMtxData()[i+(arg._i1-_i1)]+arg._j1), 1,(void*) ((COMPLEXFLOAT*)this->getMtxData()[i]+_j1), 1);
		break;
	case CPLX_DOUBLE:

		for (int i = _i1; i<= _i2;i++)
			cblas_zcopy(N, (void*) ((COMPLEXDOUBLE*) arg.getMtxData()[i+(arg._i1-_i1)]+arg._j1), 1, (void*) ((COMPLEXDOUBLE*) this->getMtxData()[i]+_j1), 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix copy operation");
	}
	return *this;

}

// overload the multiplication by matrix operator -> no slicing possible! Whole matrices are multiplied!
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator*(const Matrix<_Tp>& arg) const {

	if (this->getDim_j() != arg.getDim_i()) {
		MB_OUT_ERR(" matrix matrix multiply. Matrix dimension do not agree!",__FILE__,__LINE__);
		throw std::length_error("Matrix Dimensions do not agree.");
	}

	int L = arg.getDim_j();
	Matrix<_Tp> res(this->getDim_i(), L);

	COMPLEXFLOAT alpha_f = 1.0;
	COMPLEXDOUBLE alpha_d = 1.0;
	COMPLEXDOUBLE beta_d = 0.;
	COMPLEXFLOAT  beta_f = 0.;
	switch (gettype<_Tp>()) {
	case FLT:

		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), 1.f,
				(float*) *this->getMtxData(), this->getDim_i(),
				(float*) *arg.getMtxData(), arg.getDim_i(), 0.f,
				(float*) *res.getMtxData(), this->getDim_i());
		break;
	case DBL:
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), 1.,
				(double*) *this->getMtxData(), this->getDim_i(),
				(double*) *arg.getMtxData(), arg.getDim_i(), 0.,
				(double*) *res.getMtxData(), this->getDim_i());
		break;
	case CPLX_FLOAT:
		cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), (void*) &alpha_f,
				(void*) *this->getMtxData(), this->getDim_i(),
				(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_f,
				(void*) *res.getMtxData(), this->getDim_i());
		break;

	case CPLX_DOUBLE:
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->getDim_i(),
				arg.getDim_j(), this->getDim_j(), (void*) &alpha_d,
				(void*) *this->getMtxData(), this->getDim_i(),
				(void*) *arg.getMtxData(), arg.getDim_i(), (void*) &beta_d,
				(void*) *res.getMtxData(), this->getDim_i());
		break;
	default:
		throw std::domain_error("Unsupported data type for matrix-matrix mutiplication!");

	}

	return res;
}

//overload the addition operator -> A lot of copying is done at the end!
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator+(const Matrix<_Tp>& arg) const {

	if (_i2-_i1 != arg._i2-arg._i1
			|| _j2-_j1 != arg._j2-arg._j1) {
		MB_OUT_ERR(" copy assignment operator! Mtx dimensions mismatch!",__FILE__,__LINE__);
		throw std::length_error("Matrix Dimensions do not aggree.");
	}
	std::vector<int> slice = arg.getSliceVector();

	//the below operation invokes the copy constructor!
	// the result matrix will be of size
	// (_i2-_i1 + 1) x (_j2 - _j1 + 1) , i.e. the size of the current slice of this and arg!

	Matrix<_Tp> res = arg; // res is not sliced !

	COMPLEXDOUBLE alpha_d = 1.;
	COMPLEXFLOAT  alpha_f = 1.;
	int N = _j2-_j1+1;

	switch (gettype<_Tp>()) {
	case FLT:
		for(int i = _i1; i <= _i2; i++)
			cblas_saxpy(N, 1.0, (float*) (this->getMtxData()[i]+_j1),	1, (float*) res.getMtxData()[i-_i1], 1);

		break;
	case DBL:
		for(int i = _i1; i <= _i2; i++)
			cblas_daxpy(N, 1.0, (double*) (this->getMtxData()[i]+_j1),1, (double*) res.getMtxData()[i-_i1], 1);
		break;
	case CPLX_FLOAT:
		for(int i = _i1; i <= _i2; i++)
			cblas_caxpy(N, (void*) &alpha_f,(void*) ((COMPLEXFLOAT*) (this->getMtxData()[i]+_j1)), 1,(void*) (COMPLEXFLOAT*) res.getMtxData()[i-_i1], 1);
		break;
	case CPLX_DOUBLE:
		for(int i = _i1; i <= _i2; i++)
			cblas_zaxpy(N, (void*) &alpha_d,(void*) ( (COMPLEXDOUBLE*)(this->getMtxData()[i]+_j1)), 1,(void*) (COMPLEXDOUBLE*) res.getMtxData()[i-_i1], 1);
		break;
	default:
		throw std::domain_error(
				"Unsupported data type for matrix-matrix addition!");

	}
	return res;
}

//overload the subtraction operator
template<typename _Tp>
Matrix<_Tp> Matrix<_Tp>::operator-(const Matrix<_Tp>& arg) {

	if (_i2-_i1 != arg._i2-arg._i1
			|| _j2-_j1 != arg._j2-arg._j1) {
		MB_OUT_ERR(" copy assignment operator! Mtx dimensions mismatch!",__FILE__,__LINE__);
		throw std::length_error("Matrix Dimensions do not aggree.");
	}

	Matrix<_Tp> res = ((_Tp) -1.) * arg;
	return (*this) + res;
}
template class MB::Matrix<float>;
template class MB::Matrix< COMPLEXFLOAT>;
template class MB::Matrix<double>;
template class MB::Matrix<COMPLEXDOUBLE >;
