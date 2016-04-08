#include <matrix/matrix.hpp>
#include <iostream>

using namespace std;
using namespace MB;

// placeholder for the desired type !

#define _TYPE complex float
//
//
//MB::Matrix<float> mMult(MB::Matrix<float> A, MB::Matrix<float> B){
//
//
//	if (A._i2-A._i1 != B._i2-B._i1
//			|| A._j2-A._j1 != B._j2-B._j1) {
//		MB_OUT_ERR(" mMult(...) ! Mtx dimensions mismatch!",__FILE__,__LINE__);
//		throw std::length_error("Argument matrix dimensions do not agree.");
//	}
//
//	MB::Matrix<float> C = A;
//
//	for(int i = 0 ; i < C.getDim_i(); i++)
//		for(int j = 0; j < C.getDim_j(); j++){
//			C(i,j) = A(i+A._i1,j+A._j1)*B(i+B._i1,j+B.j1);
//		}
//}

void compute(){

	MB::Matrix<_TYPE> m1 = eye<_TYPE>(4);
	MB::Matrix<_TYPE> m2 = 2*eye<_TYPE>(4);

	std::cout << "Assignment operator testing!\n";
	std::cout << "M1 before: \n";
	std::cout << m1;
	m1(1,1) = -5;
	m1.setSlice(2,3,2,3);
	m1 = 3*eye<_TYPE>(2);
	std::cout << "M1 after \n";
	std::cout<< m1;

	std::cout << "Arithmetic operations test!\n";
	m2(0,0) =-1;
	m2(0,1) =-1;
	m2(1,0) =-1;
	m2(1,1) =-1;
	m2.setSlice(0,1,0,1);
	std::cout<<"M1 slice:\n";
	std::cout<< m1 ;
	std::cout<<"M2 slice:\n";
	std::cout << m2;

	MB::Matrix<_TYPE> m3 = m1; //+m2;
	m3.resetSlice();
	std::cout<<"Matrix M3 = M1 slice:\n";
	std::cout << m3;
	m3 = m3+m2;
	std::cout<<" M3 = M3+M2 slice:\n";
	std::cout<<m3;

	m3 = m3-2*m2;
	std::cout<<" M3 = M3-2*M2 slice:\n";
	std::cout<<m3;

	m3 = mMult(m1,m2);

	std::cout<<" M3 = M1.*M2 slice:\n";
	std::cout<< m3;



}

//
//
//MB::Matrix<float> mMult(MB::Matrix<float> A, MB::Matrix<float> B){
//
//
//	if (A._i2-A._i1 != B._i2-B._i1
//			|| A._j2-A._j1 != B._j2-B._j1) {
//		MB_OUT_ERR(" mMult(...) ! Mtx dimensions mismatch!",__FILE__,__LINE__);
//		throw std::length_error("Argument matrix dimensions do not agree.");
//	}
//
//	MB::Matrix<float> C = A;
//
//	for(int i = 0 ; i < C.getDim_i(); i++)
//		for(int j = 0; j < C.getDim_j(); j++){
//			C(i,j) = A(i+A._i1,j+A._j1)*B(i+B._i1,j+B.j1);
//		}
//}

void swaprows(){

	MB::Matrix<_TYPE> m1 = eye<_TYPE>(4);
	std::cout<< "M1 before swap: \n";
	std::cout << m1;
	_TYPE* tmp1 = m1.setRowPtr(2,m1.getMtxData()[0]);
	m1.setRowPtr(0,tmp1);
	std::cout<< "M1 after swap: \n";
	std::cout << m1;

}

int main() {

//	compute();
	swaprows();
	std::cout<<"Goodbye world\n";

}

