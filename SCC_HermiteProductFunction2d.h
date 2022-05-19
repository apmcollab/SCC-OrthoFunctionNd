
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>


#ifndef  SCC_HERMITE_PRODUCT_FUNCTION_2D_
#define  SCC_HERMITE_PRODUCT_FUNCTION_2D_
//
//######################################################################
//                SCC_HermiteProductFunction2d.h
//######################################################################
//
/*

The HermiteProductFunction2d creates functions from R^2 -> R that
are products of two Hermite functions.

If gammaX, and gammaY are double values greater than  0 and nx and ny
are integer values all greater than or equal to zero, then these functions
are products of functions of the form

H_n(s) = C*exp(-((gamma*s)^2)/2)*P_n(s*gamma)

where C = sqrt(gamma/ [sqrt(pi)*(2^n)*n!] ), P_n is the nth
Hermite polynomial, and n >= 0.

When s = x, gamma = gammaX and n = nx
     s = y, gamma = gammaY and n = ny

These functions are orthonormal over [-oo, oo]X[-oo, oo].

If alphaX and alphaY and betaX and betaY are real values such that
alphaX*betaX  < 0 and alphaY*betaY  < 0 then these functions
are eigenfunctions of the separable operator

L = alphaX*d^/dx^2 + alphaY*d^/dy^2 + betaX*x^2 + betaY*y^2

when gammaX = (betaX/alphaX)^(1/4) and gammaY = (betaY/alphaY)^(1/4)

The eigenvalue is given by

lambda  = [ -(2*nX + 1)*sign(alphaX)*(|alphaX*betaX|)^(1/2) ] +
          [ -(2*nY + 1)*sign(alphaY)*(|alphaY*betaY|)^(1/2) ]

e.g. the sum of the eigenvalues of each of the separable components.

By specifying shiftX and shiftY as additional parameters in the
constructor (or initializer) the center of the functions created can be
shifted to arbitrary locations (shiftX, shiftY) in R^2.

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
#include "SCC_HermiteFunction.h"

namespace SCC
{

class HermiteProductFunction2d
{
    public :

    HermiteProductFunction2d()
	{
    initialize();
	}

    HermiteProductFunction2d(double gammaX, double gammaY)
    {
    initialize(gammaX,gammaY,shiftX,shiftY);
    }

    HermiteProductFunction2d(double gammaX, double gammaY, double shiftX, double shiftY)
    {
    initialize(gammaX,gammaY,shiftX,shiftY);
    }

    virtual ~HermiteProductFunction2d(){};

    void initialize()
	{
    gammaX = gammaY = 1.0;
    shiftX = shiftY = 0.0;
    orthoPoly.initialize(OrthoPoly::Hermite);
	}

    void initialize(double gammaX, double gammaY)
    {
    shiftX = shiftY = 0.0;
    initialize(gammaX,gammaY,shiftX,shiftY);
    }


    void initialize(double gammaX, double gammaY,  double shiftX, double shiftY)
    {
    this->gammaX = gammaX;
    this->shiftX = shiftX;

    this->gammaY = gammaY;
    this->shiftY = shiftY;

    orthoPoly.initialize(OrthoPoly::Hermite);
    }

    std::function< double(double,double) > getHermiteProductFunction2d(long nX,long nY) const
    {
    double sqrtPi  = std::sqrt(3.141592653589793238);
    PolyFun PX     = orthoPoly.getNthOrthoPolynomial(nX);
    PolyFun PY     = orthoPoly.getNthOrthoPolynomial(nY);

    double normConstant;

    normConstant  = gammaX/(sqrtPi*pow(2.0,nX)*exp(std::lgamma(nX+1)));
    PX *= std::sqrt(normConstant);

    normConstant  = gammaY/(sqrtPi*pow(2.0,nY)*exp(std::lgamma(nY+1)));
    PY *= std::sqrt(normConstant);

    std::function< double(double,double) > H = [=](double x,double y)
	{
    	double val = -0.5*(gammaX*(x-shiftX)*gammaX*(x-shiftX) + gammaY*(y-shiftY)*gammaY*(y-shiftY));
    	return std::exp(val)*PX(gammaX*(x-shiftX))*PY(gammaY*(y-shiftY));
	};

    return H;
    }

    private :

    double     gammaX,gammaY;
    double     shiftX,shiftY;
    OrthoPoly      orthoPoly;
};

}; // namespace SCC

#endif


