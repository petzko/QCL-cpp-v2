/*
 * utils.hpp
 *
 *  Created on: Jul 10, 2015
 *      Author: petzko
 */

#ifndef INCLUDE_UTILS_HPP_
#define INCLUDE_UTILS_HPP_

#include <complex.h>
#include <typeinfo>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

#include <map>
#include <vector>
#include <assert.h>


#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>


#define MB_OUT_WRN(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ":\n->WARNING: " << str << "\n";}

#define MB_OUT_ERR(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n->ERROR: " << str << "\n";}

#define max(x,y){ (x>y ? x:y)}
#define min(x,y){ (x<y ? x:y)}

#define explicit_cast(type,data){*((type*)&data)}
//predefine the BLAS supported datatypes -> single and double precision real and complex numbers...

const int FLT = 0;
const int DBL = 1;
const int CPLXFLT = 2;
const int CPLXDBL = 3;

std::ostream& operator<<(std::ostream& out, complex double nr);
std::ostream& operator<<(std::ostream& out, complex float nr);

template<typename _Tp>
int gettype() {

	double dbl = 0.;
	float flt = 0.;
	complex double cdbl = 0. + I;
	complex float cflt = 0. + I;

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



#endif /* INCLUDE_UTILS_HPP_ */
