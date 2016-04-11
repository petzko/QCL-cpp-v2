#include <matrix/matrix.hpp>
#include <sim/SimSettings.hpp>
#include <common/CONSTANTS.hpp>
#include <sim/SimLoop.hpp>
#include <iostream>

using namespace std;
using namespace MB;

// placeholder for the desired type !


void compute(){

	MB::Matrix<_TYPE_> m1 = eye<_TYPE_>(4);
	MB::Matrix<_TYPE_> m2 = 2*eye<_TYPE_>(4);

	std::cout << "Assignment operator te	sting!\n";
	std::cout << "M1 before: \n";
	std::cout << m1;
	m1(1,1) = -5;
	m1.setSlice(2,3,2,3);
	m1 = 3*eye<_TYPE_>(2);
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

	MB::Matrix<_TYPE_> m3 = m1; //+m2;
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

void parseTest(){

	char* simFname  = "/home/petzko/workspace/qcl-sim-cpp/Debug/OPTICAMAIN.sim";
	char* setFname  = "/home/petzko/workspace/qcl-sim-cpp/Debug/12-1-off.set";
	SimSettings set(simFname,setFname);
	set.initSimSettings();



}

void swaprows(){

	MB::Matrix<_TYPE_> m1 = eye<_TYPE_>(4);
	std::cout<< "M1 before swap: \n";
	std::cout << m1;
	_TYPE_* tmp1 = m1.setRowPtr(2,m1.getMtxData()[0]);
	m1.setRowPtr(0,tmp1);
	std::cout<< "M1 after swap: \n";
	std::cout << m1;
}

void initTest(){

	MB::Matrix<_TYPE_> m1 = ones<_TYPE_>(1,4);
	std::cout << "Ones:";
	std::cout << m1;

	MB::Matrix<_TYPE_> m2(1,4);
	std::cout << "zeros:";
	std::cout << m2;

	MB::Matrix<_TYPE_> m3=randInit<_TYPE_>(1,4);
	m3 = m3+ _TYPE_(0,1)*randInit<_TYPE_>(1,4,10);
	std::cout << "rand:";
	std::cout << m3;

	MB::Matrix<_TYPE_> m4(0,0);
//	m4.init(1,4);
	std::cout << "Zeros:";
	std::cout<<m4;




}

int main() {

//	compute();
//	swaprows();
//	parseTest();
//	initTest();

	char* simFname  = "/home/petzko/workspace/qcl-sim-cpp/Debug/OPTICAMAIN.sim";
	char* setFname  = "/home/petzko/workspace/qcl-sim-cpp/Debug/12-1-off.set";

	startSim(simFname,setFname);
	std::cout << "\n";
	std::cout<<"Goodbye world\n";

}

