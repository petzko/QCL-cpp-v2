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
	unsigned int _dimI, _dimJ, _bytes;
	unsigned int _i1, _i2,_j1,_j2; // current slice indices!
	unsigned int _len; // store the total length.
	_Tp** _data;

public:

	_Tp* setRowPtr(unsigned int fromIdx,_Tp* toPtr);
	unsigned int getLen(){return _len;}
	void init(unsigned int dim_i, unsigned int dim_j);

	Matrix(){
		;
	}

	unsigned int getDim_i() const;
	unsigned int getDim_j() const;
	unsigned int getBytesize() const {
		return _bytes;
	}
	_Tp** const getMtxData() const;
	bool square() const {
		return _dimI == _dimJ;
	}

	//	 initialize an empty matrix
	Matrix(unsigned int dim_i, unsigned int dim_j);

	//	 copy constructor
	Matrix(const Matrix<_Tp>& arg);

	//	destructor
	virtual ~Matrix();

	/**
	 * Overload the assignment operator!! Extremely important in order to be able to correctly execute arithmetic operations on matrices with  syntax!!!
	 */
	Matrix<_Tp> & operator=( Matrix<_Tp> const &  arg);

	/**
	 * Return by reference (implicit pointer) to the i,j-th data element
	 * so that the user can actually overwrite the corresponding entry
	 */

	_Tp& operator()(unsigned int const i, unsigned int const j) const;

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
	Matrix<_Tp> getSliceMtx(unsigned int i1, unsigned int i2, unsigned int j1,
			unsigned int j2);

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
	void setSlice(unsigned int i1, unsigned int i2, unsigned int j1, unsigned int j2);
	/**
	 * equivalent to setSlice(0,_dimI,0,_dimJ);
	 */
	void resetSlice();

	std::vector<unsigned int> getSliceVector() const;
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

	void operator()(Matrix<_Tp> inMatrix, unsigned int i1, unsigned int i2, unsigned int j1,
			unsigned int j2);

	// overload the multiplication by matrix operator...
	Matrix<_Tp> operator*(const Matrix<_Tp>& arg) const;

	//overload the addition operator
	Matrix<_Tp> operator+(const Matrix<_Tp>& arg) const;

	//overload the subtraction operator
	Matrix<_Tp> operator-(const Matrix<_Tp>& arg);

	//overload the print operator -> PRINTS THE WHOLE MATRIX!!
	friend std::ostream& operator<<(std::ostream& os,
			const MB::Matrix<_Tp>& obj) {

		std::vector<unsigned int> slice = obj.getSliceVector();

		for (int i = slice[0]; i <= slice[1]; i++) {
			for (int j = slice[2]; j <= slice[3]; j++)
				os << obj.getMtxData()[i][j] << " ";
			os << "\n";
		}
		return os;
	}

	friend Matrix<_Tp> operator*(_Tp lhs, const Matrix<_Tp>& rhs) {
		MB::Matrix<_Tp> res = rhs;

		std::vector<unsigned int> slice = rhs.getSliceVector();
		unsigned int N = slice[3]-slice[2]+1;

		//		std::cout << "slice: " << slice[0] << " " << slice[1] << " " << slice[2] << " " << slice[3]  << " \n";

		int type = gettype<_Tp>();

		switch (type) {
		case FLT:
			for (int i = slice[0]; i <= slice[1]; i++)
				cblas_sscal(N, explicit_cast(float,lhs), (float*) (res._data[i] +  slice[2]),1);
			break;
		case DBL:
			for (int i = slice[0]; i <= slice[1]; i++)
				cblas_dscal(N , explicit_cast(double,lhs), (double*) (res._data[i]+slice[2]),	1);
			break;
		case CPLXFLT:
			for (int i = slice[0]; i <= slice[1]; i++)
				cblas_cscal(N , (void*) &lhs, (void*)  (res._data[i]+slice[2]), 1);
			break;
		case CPLXDBL:
			for (int i = slice[0]; i <= slice[1]; i++)
				cblas_zscal(N, (void*) &lhs, (void*)  (res._data[i]+slice[2]), 1);
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
MB::Matrix<_Tp> eye(unsigned int M) {

	MB::Matrix<_Tp> res(M, M);
	for (int i = 0; i < M; i++)
		res(i, i) = 1.0;

	return res;
}

template<typename _Tp>
MB::Matrix<_Tp> mMult(MB::Matrix<_Tp> A, MB::Matrix<_Tp> B){

	std::vector<unsigned int> Aslice = A.getSliceVector();
	std::vector<unsigned int> Bslice = B.getSliceVector();

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
