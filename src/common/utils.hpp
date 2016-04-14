/*
 * utils.hpp
 *
 *  Created on: Jul 10, 2015
 *      Author: petzko
 */

#ifndef INCLUDE_UTILS_HPP_
#define INCLUDE_UTILS_HPP_

#include <typeinfo>
#include <iostream>
#include <sstream>

#include <stdexcept>
#include <stdlib.h>

/* time */
#include <time.h>

#include <map>
#include <vector>
#include <algorithm>

#include <assert.h>

#define CPP_CPLX
#ifdef CPP_CPLX

#include <complex>
#define COMPLEXDOUBLE std::complex<double>
#define COMPLEXFLOAT std::complex<float>

/**
 * redefine the imaginary unit to be on par with complex.h
 */
#define I _TYPE_(0,1)

#define REAL(x) std::real(x)
#define IMAG(x) std::imag(x)
#define ABS(x) std::abs(x)


#else

#include <complex.h>

#define COMPLEXDOUBLE    complex double
#define COMPLEXFLOAT     complex float

#define REAL(x) (creal(x))
#define IMAG(x) (cimag(x))
#define  ABS(x)  sqrt(REAL(x)*REAL(x) + IMAG(x)*IMAG(x))

#endif


#define _TYPE_  COMPLEXDOUBLE
#define _i I


enum TYPE_ID {FLT , DBL, CPLX_FLOAT, CPLX_DOUBLE , UNDEF};


#include <gsl/gsl_cblas.h>




#define MB_OUT_WRN(str,FILE,LINE){\
		std::cout<<"FILE:" << FILE << " LINE:" << LINE << ":\n->WARNING: " << str << "\n";}

#define MB_OUT_ERR(str,FILE,LINE){\
		std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n->ERROR: " << str << "\n";}

//#define max(x,y){ (x>y ? x:y) }
//#define min(x,y){ (x<y ? x:y) }

#define explicit_cast(type,data){*((type*)&data)}





template<typename _Tp>
TYPE_ID gettype() {
	double dbl = 0.;
	float flt = 0.;

	COMPLEXDOUBLE cdbl  = (COMPLEXDOUBLE)(_i);
	COMPLEXFLOAT cflt  = (COMPLEXFLOAT)(_i);


	if (typeid(_Tp) == typeid(flt))
		return FLT;
	if (typeid(_Tp) == typeid(dbl))
		return DBL;
	if (typeid(_Tp) == typeid(cflt))
		return CPLX_FLOAT;
	if (typeid(_Tp) == typeid(cdbl))
		return CPLX_DOUBLE;
	return UNDEF;

}


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
