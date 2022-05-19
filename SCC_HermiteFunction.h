
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#ifndef  SCC_HERMITE_FUNCTION_
#define  SCC_HERMITE_FUNCTION_
//
//######################################################################
//                    SCC_HermiteFunction.h
//######################################################################
//
/*

Instances of HermiteFunction are the functions

H_n(x) = C*std::exp(-((gamma*x)^2)/2)*P_n(x*gamma)

where C = std::sqrt(gamma/ [std::sqrt(pi)*(2^n)*n!] ), P_n is the nth
Hermite polynomial, and n >= 0.

These functions are orthonormal over (-oo,+oo).

If alpha and beta are real values such that alpha*beta < 0, then
H_n(x) are eigenfunctions of the operator

L = alpha*d^2/dx^2 + beta*x^2

when gamma = (beta/alpha)^(1/4). The eigenvalue is given by

lambda  = -(2*n + 1)*sign(alpha)*(|alpha*beta|)^(1/2)

The specification of the shift parameter in the constructor (or
initializer) shifts the origin of the HermiteFunction to
x = shift.

 */

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

#include "SCC_PolyFun.h"
#include "SCC_OrthoPoly.h"

namespace SCC
{

class HermiteFunction
{
    public :

//
//  Constructors
//
    HermiteFunction()
	{
    gamma         = 0.0;
	shift         = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
	}

    HermiteFunction(double gamma, double shift = 0.0)
    {
    this->gamma          = gamma;
    this->shift          = shift;
    orthoPoly.initialize(OrthoPoly::Hermite,1.0,shift);
    }

    virtual ~HermiteFunction(){};

    void initialize()
    {
    gamma         = 0.0;
	shift         = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    void initialize(double gamma, double shift = 0.0)
    {
    this->gamma          = gamma;
    this->shift          = shift;
    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    void   setGamma(double gamma) {this->gamma = gamma;}
    double getGamma()             {return this->gamma;}

    void   setShift(double shift) {this->shift = shift;}
    double getshift()             {return shift;};


    std::function< double(double) > getHermiteFunction(long n) const
    {
    double sqrtPi = std::sqrt(3.141592653589793238);
    PolyFun P     = orthoPoly.getNthOrthoPolynomial(n);

    double normConstant = gamma/(sqrtPi*pow(2.0,n)*std::exp(std::lgamma(n+1)));
    P  *= std::sqrt(normConstant);

    double gamma  = this->gamma;
    double shift  = this->shift;

    std::function< double(double) > Hn = [P, gamma, shift](double x)
	{
    	return std::exp(-0.5*gamma*(x-shift)*gamma*(x-shift))*P(gamma*(x-shift));
	};

    return Hn;
    }


	std::vector< std::function< double(double) > >  getHermiteFunctionArray(long maxIndex)
    {
	std::vector< std::function< double(double) > > H;
    for(int i = 0; i < maxIndex; i++)
    {
    	H.push_back(getHermiteFunction(i));
    }

    return H;
    }

    private :

    double                    gamma;
    double                    shift;
    OrthoPoly             orthoPoly;
};

}; // namespace SCC

#endif


