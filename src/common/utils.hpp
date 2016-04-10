/*
 * utils.hpp
 *
 *  Created on: Jul 10, 2015
 *      Author: petzko
 */

#ifndef INCLUDE_UTILS_HPP_
#define INCLUDE_UTILS_HPP_

#include <complex>
#include <typeinfo>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>


#include <map>
#include <vector>
#include <algorithm>

#include <assert.h>


#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>


#define MB_OUT_WRN(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ":\n->WARNING: " << str << "\n";}

#define MB_OUT_ERR(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n->ERROR: " << str << "\n";}

//#define max(x,y){ (x>y ? x:y) }
//#define min(x,y){ (x<y ? x:y) }

#define explicit_cast(type,data){*((type*)&data)}
//predefine the BLAS supported datatypes -> single and double precision real and complex numbers...

const int FLT = 0;
const int DBL = 1;
const int CPLXFLT = 2;
const int CPLXDBL = 3;
const int INTGR = 4;
const int D_PTR = 5;


//std::ostream& operator<<(std::ostream& out, complex double nr);
//std::ostream& operator<<(std::ostream& out, complex float nr);

template<typename _Tp>
int gettype() {

	double dbl = 0.;
	float flt = 0.;
	std::complex<double>cdbl (0,1);
	std::complex<float> cflt (0,1);

	if (typeid(_Tp) == typeid(flt))
		return FLT;
	if (typeid(_Tp) == typeid(dbl))
		return DBL;
	if (typeid(_Tp) == typeid(cflt))
		return CPLXFLT;
	if (typeid(_Tp) == typeid(cdbl))
		return CPLXDBL;
	return -1;

}

template<typename _Tp>
int gettype(_Tp in) {

	double dbl = 0.;
	float flt = 0.;
	std::complex<double>cdbl (0,1);
	std::complex<float> cflt (0,1);

	int intr = 1;
	double* dptr;

	if (typeid(in) == typeid(flt))
		return FLT;
	if (typeid(in) == typeid(dbl))
		return DBL;
	if (typeid(in) == typeid(cflt))
		return CPLXFLT;
	if (typeid(in) == typeid(cdbl))
		return CPLXDBL;
	if (typeid(in) == typeid(intr))
			return INTGR;
	if (typeid(in) == typeid(cdbl))
			return D_PTR;
	return -1;

}




#endif /* INCLUDE_UTILS_HPP_ */
