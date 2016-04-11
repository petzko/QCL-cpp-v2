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
#include <sstream>

#include <stdexcept>
#include <stdlib.h>
#include <time.h>       /* time */


#include <map>
#include <vector>
#include <algorithm>

#include <assert.h>


#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

#define _TYPE_ std::complex<double>

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
//
//template<typename _Tp>
//int gettype(_Tp in) {
//
//	double dbl = 0.;
//	float flt = 0.;
//	std::complex<double>cdbl (0,1);
//	std::complex<float> cflt (0,1);
//
//	int intr = 1;
//	double* dptr;
//
//	if (typeid(in) == typeid(flt))
//		return FLT;
//	if (typeid(in) == typeid(dbl))
//		return DBL;
//	if (typeid(in) == typeid(cflt))
//		return CPLXFLT;
//	if (typeid(in) == typeid(cdbl))
//		return CPLXDBL;
//	if (typeid(in) == typeid(intr))
//			return INTGR;
//	if (typeid(in) == typeid(cdbl))
//			return D_PTR;
//	return -1;
//
//}


/**
 *
 * std::exception <exception> interface (debatable if you should catch this)
    std::bad_alloc <new> failure to allocate storage
        std::bad_array_new_length <new> invalid array length
    std::bad_cast <typeinfo> execution of an invalid dynamic-cast
    std::bad_exception <exception> signifies an incorrect exception was thrown
    std::bad_function_call <functional> thrown by "null" std::function
    std::bad_typeid <typeinfo> using typeinfo on a null pointer
    std::bad_weak_ptr <memory> constructing a shared_ptr from a bad weak_ptr
    std::logic_error <stdexcept> errors detectable before the program executes
        std::domain_error <stdexcept> parameter outside the valid range
        std::future_error <future> violated a std::promise/std::future condition
        std::invalid_argument <stdexcept> invalid argument
        std::length_error <stdexcept> length exceeds its maximum allowable size
        std::out_of_range <stdexcept> argument value not in its expected range
    std::runtime_error <stdexcept> errors detectable when the program executes
        std::overflow_error <stdexcept> arithmetic overflow error.
        std::underflow_error <stdexcept> arithmetic underflow error.
        std::range_error <stdexcept> range errors in internal computations
        std::regex_error <regex> errors from the regular expression library.
        std::system_error <system_error> from operating system or other C API
            std::ios_base::failure <ios> Input or output error
 *
 *
 */

#endif /* INCLUDE_UTILS_HPP_ */
