#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#endif

#include <functional>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#ifndef  SCC_POLY_FUN_
#define  SCC_POLY_FUN_
//
//######################################################################
//                        SCC_PolyFun.h
//######################################################################
//

/*
#############################################################################
#
# Copyright  2015 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/


namespace SCC
{

class PolyFun
{

    public :

	PolyFun() : coeff()
	{
	    degree = -1;
	    zero   = 0.0;
	}

	PolyFun(const PolyFun& P)
	{
	    degree = P.degree;
	    coeff  = P.coeff;
	    zero   = P.zero;
	}

	PolyFun(int degree)
	{
	    this->degree = degree;
	    coeff.resize(degree+1,0.0);
	    zero = 0.0;
	}

	PolyFun(const std::vector<double>& coefficients)
	{
	    degree = coefficients.size()-1;
	    coeff  = coefficients;
	    zero   = 0.0;
	}

	~PolyFun()
	{}

	void initialize()
	{
	    coeff.clear();
	    degree = -1;
	    zero   = 0.0;
	}

	void initialize(const PolyFun& P)
	{
		initialize();
	    coeff  = P.coeff;
	    degree = P.degree;
	}

	void initialize(int degree)
	{
		initialize();
	    this->degree = degree;
	    coeff.resize(degree+1,0.0);
	}

	void initialize(const std::vector<double>& coefficients)
	{
		initialize();
	    degree = coefficients.size()-1;
	    coeff      = coefficients;
	}

	PolyFun integrate()
	{
	    int i;
	    PolyFun R(degree+1);

	    R.coeff[0] = 0.0;                     // pick arbitrary constant
	                                      // to be zero
	    for(i = 1; i <= degree + 1; i++)
	    {R.coeff[i] =  coeff[i-1]/double(i);}

	    return R;
	}

	double integrate(double Xmin, double Xmax)
	{
	    PolyFun R = this->integrate();
	    return (R(Xmax) - R(Xmin));
	}


	PolyFun differentiate()
	{
	    int i;
	    double rDegree;

	    if(degree > 0)
	    {rDegree = degree - 1;}
	    else
	    {rDegree = 0;}

	    PolyFun R(rDegree);

	    if(degree > 0)
	    {
	    for(i = 0; i <= degree - 1; i++)
	    {R.coeff[i] =  coeff[i+1]*double(i+1);}
	    }
	    else
	    {
	    R.coeff[0] = 0.0;
	    }

	    return R;
	}

	void operator=(const PolyFun& P)
	{
	    degree = P.degree;
	    coeff  = P.coeff;
	}

	PolyFun operator+(const PolyFun& P)
	{
	    PolyFun  R;
	    int i;

	    if(degree >= P.degree)
	    {
	      R = *this;
	      for(i = 0; i <= P.degree; i++)
	      {R.coeff[i] += P.coeff[i];}
	    }
	    else
	    {
	       R = P;
	       for(i = 0; i <= degree; i++)
	       {R.coeff[i] += coeff[i];}
	    }
	    return R;
	}

	PolyFun operator-(const PolyFun& P)
	{
	    PolyFun  R;
	    int i;

	    if(degree >= P.degree)
	    {
	      R = *this;
	      for(i = 0; i <= P.degree; i++)
	      {R.coeff[i] -= P.coeff[i];}
	    }
	    else
	    {
	       R  =    P;
	       R *= -1.0;
	       for(i = 0; i <= degree; i++)
	       {R.coeff[i] += coeff[i];}
	    }
	    return R;
	}
	PolyFun operator*(const PolyFun& P)
	{
	    int i; int k;
	    PolyFun R(degree+P.degree);

	    for(i = 0; i <= degree; i++)
	    {
	    for(k = i; k <= i+P.degree; k++)
	    {R.coeff[k] += coeff[i]*P.coeff[k-i];}
	    }

	    return R;
	}

	PolyFun operator-()
	{
	    PolyFun  R(*this);
	    int i;
	    for(i = 0; i <= degree; i++)
	    {R.coeff[i] = -coeff[i];}
	    return R;
	}

	double operator()(double x) const
	{
	//  Evaluate using Horner's method
	//
	    double result = 0.0;
	    int i;
	    for(i = degree; i>=0; i--)
	    {result = x*result + coeff[i];}
	    return result;
	}

	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return F;
	}


	const double&  operator[](long i) const
	{
	   if(i > degree) {return zero;}
	   return coeff[i];
	}

	double&  operator[](long i)
	{
	   if(i > degree)
	   {degree = i; coeff.resize(degree+1,0.0);}
	   return coeff[i];
	}

	PolyFun operator*(double alpha)
	{
	   PolyFun R(*this);
	   R *= alpha;
	   return R;
	}

	void operator*=(double alpha)
	{
	   for(int i =0; i <= degree; i++) coeff[i] *= alpha;
	}

	void operator/=(double alpha)
	{
	   for(int i =0; i <= degree; i++) coeff[i] /= alpha;
	}


	PolyFun operator/(double alpha)
	{
	   PolyFun R(*this);
	   R /= alpha;
	   return R;
	}

	friend PolyFun operator*(double alpha, const PolyFun& P)
	{
	   PolyFun R(P);
	   R *= alpha;
	   return R;
	}


	// Shift the independent variable x -> x - p

	PolyFun shift(double p)
	{
	     PolyFun R(degree);

	     PolyFun S(1);
	     S[0] = -p; S[1] = 1.0;      // S = (x - p)
	     PolyFun Q(S);

	     R[0] = coeff[0];
	     for(long i = 1; i <= degree; i++)
	     {
	     R = R + coeff[i]*Q;
	     Q = Q*S;
	     }
	     return R;
	}


	// Scale the independent variable x -> alpha*x

	PolyFun scale(double alpha)
	{
	     PolyFun R(*this);
	     double s = 1.0;
	     for(long i = 1; i <= R.degree; i++)
	     {
	     s = s*alpha;
	     R.coeff[i] = R.coeff[i]*s;
	     }
	     return R;
	}

	//

	friend std::ostream& operator <<(std::ostream& outStream, const PolyFun& P)
	{
	    if(P.degree >= 0)
	    {outStream << "1   : " << P.coeff[0] << std::endl;}
	    int i;
	    if(P.degree >= 1)
	    {
	    for(i = 1; i <= P.degree; i++)
	    {
	    outStream << "x^" << i << " : " << P.coeff[i] << std::endl;;
	    }}

	    return outStream;
	}

	/*
    friend std::ostream& operator <<(std::ostream& outStream, const PolyFun& P)
	{
		std::ios_base::fmtflags ff;
        int precisionCache;

		ff =  outStream.flags();
        precisionCache = outStream.precision();

	    outStream.setf(std::ios::fixed);
	    outStream.precision(2);
	    if(P.degree >= 0)
	    {outStream << std::setw(2) << P.coeff[0];}
	    int i;
	    if(P.degree >= 1)
	    {
	    for(i = 1; i <= P.degree; i++)
	    {
	    outStream << " + " << std::setw(2) << P.coeff[i] << "x^";
	    outStream.setf(std::ios::left);
	    outStream << std::setw(2) << i;
	    outStream.setf(std::ios::right);
	    }
	    }
	    outStream.flags(ff);
	    outStream.precision(precisionCache);

	    return outStream;
	}
    */


	int getDegree() const
	{
	     long i;
	     int deg = 0;
	     for(i = degree; i >= 0; i--)
	     {
	     if((deg == 0)&&(coeff[i] != 0.0)) deg = i;
	     }
	     return deg;
	}

	std::vector<double> getCoefficients() const
	{
	     return coeff;
	}


    bool operator==(const PolyFun& P)
	{
    	double epsTol     = 1.0e-14;
    	bool returnValue  = true;
    	int minDegree     = (coeff.size() < P.coeff.size()) ? coeff.size()-1 : P.coeff.size() -1;

    	for(int i = 0; i <= minDegree; i++)
    	{
    	if(std::abs(coeff[i] - P.coeff[i]) > epsTol) returnValue = false;
    	}

    	for(int i = minDegree+1; i < (int)coeff.size(); i++)
    	{
    	if(std::abs(coeff[i]) > epsTol) returnValue = false;
    	}

    	for(int i = minDegree+1; i < (int)P.coeff.size(); i++)
    	{
    	if(std::abs(P.coeff[i]) > epsTol) returnValue = false;
    	}

    	return returnValue;
	}

    bool operator!=(const PolyFun& P)
	{
    	return not this->operator==(P);
	}

    private :

    std::vector<double>     coeff;
    int               degree;
    double              zero;  // reference to 0 coefficient value
};

} // Namespace SCC
#endif
  
