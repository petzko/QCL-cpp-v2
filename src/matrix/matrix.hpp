#ifndef _MATRIX_
#define _MATRIX_

#include <common/utils.hpp>

/***************************************************************
 * Usual 2D matrix of size dimI x dimJ
 ***************************************************************/
namespace MB {

template<typename _Tp>
class Matrix {

private:
	int _dimI, _dimJ, _bytes;
	int _i1, _i2,_j1,_j2; // current slice indices!
	int _len; // store the total length.
	_Tp** _data;

public:

	_Tp* setRowPtr(int fromIdx,_Tp* toPtr);
	int getLen(){return _len;}
	void init(int dim_i, int dim_j);

	Matrix();

	int getDim_i() const;
	int getDim_j() const;
	int getBytesize() const {
		return _bytes;
	}
	_Tp** const getMtxData() const;
	bool square() const {
		return _dimI == _dimJ;
	}

	//	 initialize an empty matrix
	Matrix(int dim_i, int dim_j);

	//	 copy constructor
	Matrix(const Matrix<_Tp>& arg);

	//	destructor
	virtual ~Matrix();

	/**
	 * Overload the assignment operator!! Extremely important in order to be able to correctly execute arithmetic operations on matrices with  syntax!!!
	 */
	Matrix<_Tp>& operator=( Matrix<_Tp> const &  arg);

	/**
	 * Return by reference (implicit pointer) to the i,j-th data element
	 * so that the user can actually overwrite the corresponding entry
	 */

	_Tp& operator()(int const i, int const j) const;

	/**
	 * 	Routine to retreive the submatrix locked between (inclusive) row and col indices
	 * 	i1:i2 x j1:j2
	 *
	 *	@param i1 - the row index of the upper left element  of the	submatrix
	 *  @param j1 - the col index of the upper left element  of the	submatrix
	 *  @param i2 - the row index of the lower right element of the submatrix;
	 *  @param j2 - the col index of the lower right element of the submatrix;
	 *  @return - a submatrix object of the current matrix with indices (i1,j1)-->(i2,j2)!
	 *  the returned submatrix contains a STAND ALONE INDEPENDENT copy of the data contained in the original matrix between indices (i1,j1) and (i2,j2)
	 */
	Matrix<_Tp> getSliceMtx(int i1, int i2, int j1,
			int j2);

	/**
	 *
	 * Set the slice indices to the following elements! ALl subsequent operations on the matrix will modify ONLY the
	 * values within the slice and non outside. A call of resetSlice sets the slice indices to the original whole matrix, i.e.
	 * i1 = 0 i2 = _dimI-1; j1 = 0; j2 = _dimJ-1;
	 *
	 *	@param i1 - the row index of the upper left element  of the	slice
	 *  @param j1 - the col index of the upper left element  of the	slice
	 *  @param i2 - the row index of the lower right element of the slice;
	 *  @param j2 - the col index of the lower right element of the slice;
	 */
	void setSlice(int i1, int i2, int j1, int j2);
	/**
	 * equivalent to setSlice(0,_dimI,0,_dimJ);
	 */
	void resetSlice();

	std::vector<int> getSliceVector() const;
	/**
	 * copy the data stored in inMatrix into corresponding positions specified by the pairs of indices (i1,j1) and (i2,j2). If the dimensions of inMatrix are not
	 * (i2-i1)x(j2-j1) then   the function throws a length error exception.
	 *	@param inMatrix - a matrix of size (i2-i1) x ( j2 - j1) which will be copied into the corresponding place of the original data array.
	 *	@param i1 - the row index of the upper left element  of the	submatrix
	 *  @param j1 - the col index of the upper left element  of the	submatrix
	 *  @param i2 - the row index of the lower right element of the submatrix;
	 *  @param j2 - the col index of the lower right element of the submatrix;
	 *
	 */

	void operator()(Matrix<_Tp> inMatrix, int i1, int i2, int j1,
			int j2);

	// overload the multiplication by matrix operator...
	Matrix<_Tp> operator*(const Matrix<_Tp>& arg) const;

	//overload the addition operator
	Matrix<_Tp> operator+(const Matrix<_Tp>& arg) const;

	//overload the subtraction operator
	Matrix<_Tp> operator-(const Matrix<_Tp>& arg);

	//overload the print operator -> PRINTS THE WHOLE MATRIX!!
	friend std::ostream& operator<<(std::ostream& os,
			const MB::Matrix<_Tp>& obj) {

		std::vector<int> slice = obj.getSliceVector();



		if(gettype<_Tp>() == CPLX_DOUBLE || gettype<_Tp>() ==CPLX_FLOAT )
		{
			for (int i = slice[0]; i <= slice[1]; i++) {
				for (int j = slice[2]; j <= slice[3]; j++){os << "(" << REAL((_TYPE_)obj.getMtxData()[i][j]) << "," <<IMAG((_TYPE_)obj.getMtxData()[i][j]) <<") ";}
				os << "\n";
			}

		}else{
			for (int i = slice[0]; i <= slice[1]; i++) {
				for (int j = slice[2]; j <= slice[3]; j++){os << obj.getMtxData()[i][j] << " ";}
				os << "\n";
			}
		}

		return os;
	}

	friend Matrix<_Tp> operator*(_Tp lhs, const Matrix<_Tp>& rhs) {

		MB::Matrix<_Tp> res = rhs; // res is not sliced!
		std::vector<int> slice = rhs.getSliceVector();
		int N = slice[3]-slice[2]+1;
		int M = slice[1] - slice[0]+1;
		TYPE_ID type = gettype<_Tp>();
		switch (type) {
		case FLT:
			for (int i = 0; i <M; i++)
				cblas_sscal(N, explicit_cast(float,lhs), (float*) res._data[i],1);
			break;
		case DBL:
			for (int i = 0; i <M; i++)
				cblas_dscal(N , explicit_cast(double,lhs), (double*) res._data[i],	1);
			break;
		case CPLX_FLOAT:
			for (int i = 0; i <M; i++)
				cblas_cscal(N , (void*) &lhs, (void*)  res._data[i], 1);
			break;
		case CPLX_DOUBLE:
			for (int i = 0; i <M; i++)
				cblas_zscal(N, (void*) &lhs, (void*) res._data[i], 1);
			break;
		default:
			throw std::domain_error("Unsupported matrix scale operation");
		}
		return res;
	}

	friend Matrix<_Tp> operator*(Matrix<_Tp> lhs, _Tp rhs) {
		return rhs * lhs;
	}

}
;
}

template<typename _Tp>
MB::Matrix<_Tp> eye(int M) {

	MB::Matrix<_Tp> res(M, M);
	for (int i = 0; i < M; i++)
		res(i, i) = 1.0;

	return res;
}

template<typename _Tp>
MB::Matrix<_Tp> realVal(const MB::Matrix<_Tp>& m) {

	MB::Matrix<_Tp> res = m;

	for (int i = 0; i < res.getDim_i(); i++)
		for (int j = 0; j < res.getDim_j(); j++)
			res(i, j) = REAL(res(i,j));

	return res;
}


template<typename _Tp>
MB::Matrix<_Tp> imagVal(const MB::Matrix<_Tp>& m) {

	MB::Matrix<_Tp> res = m;
	for (int i = 0; i < res.getDim_i(); i++)
		for (int j = 0; j < res.getDim_j(); j++)
			res(i, j) = IMAG(res(i,j));

	return res;

}


template<typename _Tp>
MB::Matrix<_Tp> conjVal(const MB::Matrix<_Tp>& m) {

	MB::Matrix<_Tp> res = m;
	for (int i = 0; i < res.getDim_i(); i++)
		for (int j = 0; j < res.getDim_j(); j++)
			res(i, j) = REAL(res(i,j)) - _i*IMAG(res(i,j));

	return res;

}


template<typename _Tp>
MB::Matrix<_Tp> ones(int M,int N) {

	MB::Matrix<_Tp> res(M,  N);

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			res(i, j) = 1.0;

	return res;
}

template<typename _Tp>
MB::Matrix<_Tp> randInit(int M,int N, int seed  = 0) {

	/* initialize random seed: */
	//	int time_ui = time(NULL);
	srand (time(NULL) + seed);
	MB::Matrix<_Tp> res(M, N);

	int tp = gettype<_Tp>();
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			res(i, j) = rand()/(1.0*RAND_MAX);

	return res;
}

template<typename _Tp>
MB::Matrix<_Tp> mMult(MB::Matrix<_Tp> A, MB::Matrix<_Tp> B){

	std::vector<int> Aslice = A.getSliceVector();
	std::vector<int> Bslice = B.getSliceVector();

	if (Aslice[1]-Aslice[0] != Bslice[1]-Bslice[0]
												|| Aslice[3]-Aslice[2] != Bslice[3]-Bslice[2]) {
		MB_OUT_ERR(" mMult(...) ! Mtx dimensions mismatch!",__FILE__,__LINE__);
		throw std::length_error("Argument matrix dimensions do not agree.");
	}

	MB::Matrix<_Tp> C = A;

	for(int i = 0 ; i < C.getDim_i(); i++)
		for(int j = 0; j < C.getDim_j(); j++){
			C(i,j) = A(i+Aslice[0],j+Aslice[2])*B(i+Bslice[0],j+Bslice[2]);
		}
	return C;
}

#endif
