
#include <iostream>
#include <vector>

#include "SCC_PolyFun.h"

#ifdef _MSC_VER
#include "iso646.h"          // So "and" is equivalenced to &&
typedef unsigned int uint;   // Define uint to be unsigned int
#endif

#ifndef  SCC_ORTHO_POLY_
#define  SCC_ORTHO_POLY_
//
//######################################################################
//                        SCC_OrthoPoly.h
//######################################################################
//

/*
OrthoPoly::Legendre

Polynomials orthogonal over the interval [-1,1]
with the weight w(x) = 1

Integral over (-1,1) of P[n]*P[n]  = 2/(2*n + 1)

=================================================

OrthoPoly::Hermite

Polynomials orthogonal over the interval (-oo,+oo)
with the weight w(x) = exp(-x*x)

Integral over (-oo,+oo) of exp(-x*x)*P[n]*P[n]  =  sqrt(pi)*(2^n)*n!


Note: To obtain polynomials orthonormal with respect to exp(-(x*x)/2.0)
specify a xScaleFactor of 1/sqrt(2) then multiply each by
2^(index/2), e.g. P[i] *= pow(2.0,-i/2.0); For these polynomials

Integral over (-oo,+oo) of exp(-x*x)*P[n]*P[n]  =  sqrt(pi)*sqrt(2)*n!

=================================================

OrthoPoly::Laguerre

Polynomials orthogonal over the interval [0, oo )
with the weight w(x) = exp(-x)

Integral over (0,+oo) of exp(-x)*P[n]*P[n]  =  (n!)^2

=================================================

OrthoPoly::Chebyshev

Polynomials orthogonal over the interval [-1,1]
with the weight w(x) = 1/sqrt(1 - x^2)


Integral over (-1,1) of

1/sqrt(1 - x^2)*P[0]*P[0]  =  pi

1/sqrt(1 - x^2)*P[n]*P[n]  =  pi/2  n != 0

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



namespace SCC
{

class OrthoPoly
{
    public :

	enum {Legendre, Hermite, Laguerre, Chebyshev}; // Valid PolyTypes
//
//  Constructors
//
    OrthoPoly()
	{
	shift = 0.0;
    xScaleFactor   = 1.0;
    polyType      = -1;
	}

    OrthoPoly(int polyType)
    {
    shift     = 0.0;
    xScaleFactor       = 1.0;
    this->polyType    = polyType;
    }

    OrthoPoly(int polyType, double xScaleFactor,double shift)
    {
    this->shift  = shift;
    this->xScaleFactor   = xScaleFactor;
    this->polyType       = polyType;
    }

    virtual ~OrthoPoly(){};

    void initialize()
    {
   	shift = 0.0;
    xScaleFactor   = 1.0;
    polyType      = -1;
    }

    void initialize(int polyType)
    {
    shift     = 0.0;
    xScaleFactor       = 1.0;
    this->polyType    = polyType;
    }

    void initialize(int polyType, double xScaleFactor,double shift)
    {
    this->shift  = shift;
    this->xScaleFactor    = xScaleFactor;
    this->polyType       = polyType;
    }


    // The interval of definition of the polynomial is scaled first then shifted

    void setIntervalScaleAndShift(double xScaleFactor,double shift)
    {
    this->shift  = shift;
    this->xScaleFactor    = xScaleFactor;
    }

    void setType(int polyType)
    {
    	this->polyType = polyType;
    	switch(polyType)
        {
    	case OrthoPoly::Legendre  : return;
    	case OrthoPoly::Hermite   : return;
    	case OrthoPoly::Laguerre  : return;
    	case OrthoPoly::Chebyshev : return;
    	default : std::cout << " Orthogonal Polynomial Type Index not Supported " << std::endl;
        }
    }

    double getXscaleFactor() const    {return xScaleFactor;};
    double getshift()        const    {return shift;};

    PolyFun getNthOrthoPolynomial(long index) const
    {
	PolyFun Xpoly(1);
	Xpoly[0] = 0.0; Xpoly[1] = 1.0;

	PolyFun Pnm1;
	PolyFun Pn;
	PolyFun Pnp1;

	switch(polyType)
    {
//
// Legendre Functions :
//
// P_0      = 1
// P_1      = x
// P_(n+1) = ((2*n +1)*x*P_(n) + n*P_(n-1))/(n+1)
//
    case OrthoPoly::Legendre :
    {
    Pnm1.initialize(0);
    Pnm1[0] = 1.0;

    if(index == 0) return Pnm1;

    Pn.initialize(1);
    Pn[0] = 0.0; Pn[1] = 1.0;

    if(index == 1) return Pn;

    for(int i = 1; i < index; i++)
    {
    Pnp1 =((2.0*i + 1.0)*Xpoly*Pn - i*Pnm1)/double(i+1);
    Pnm1 = Pn;
    Pn   = Pnp1;
    }
    }; break;
//
// Hermite polynomials :
//
// P_0      = 1
// P_1      = 2*x
// P_(n+1) = 2*x*P_(n) - 2*n*P_(n-1)
//
    case OrthoPoly::Hermite  :
    {
    Pnm1.initialize(0);
    Pnm1[0] = 1.0;

    if(index == 0) return Pnm1;

    Pn.initialize(1);
    Pn[0] = 0.0; Pn[1] = 2.0;

    if(index == 1) return Pn;

    for(int i = 1; i < index; i++)
    {
    Pnp1 = 2.0*Xpoly*Pn - 2.0*i*Pnm1;
    Pnm1 = Pn;
    Pn   = Pnp1;
    }
    }; break;
//
// Laguerre polynomials :
//
// P_0      = 1
// P_1      = 1 - x
// P_(n+1) = (1 + 2*n - x)*P_(n) - n*n*P_(n-1)
//
    case OrthoPoly::Laguerre :
    {
    Pnm1.initialize(0);
    Pnm1[0] = 1.0;

    if(index == 0) return Pnm1;

    Pn.initialize(1);
    Pn[0] = 1.0; Pn[1] = -1.0;

    if(index == 1) return Pn;

    for(int i = 1; i < index; i++)
    {
    Pnp1 =(1.0 + 2.0*i)*Pn  - Xpoly*Pn - i*i*Pnm1;
    Pnm1 = Pn;
    Pn   = Pnp1;
    }
    }; break;
//
// Chebyshev polynomials :
//
// P_0      = 1
// P_1      = x
// P_(n+1) = (2*x)*P_(n) - P_(n-1)
//
    case OrthoPoly::Chebyshev :
    {
    Pnm1.initialize(0);
    Pnm1[0] = 1.0;

    if(index == 0) return Pnm1;

    Pn.initialize(1);
    Pn[0] = 0.0; Pn[1] = 1.0;

    if(index == 1) return Pn;

    for(int i = 1; i < index; i++)
    {
    Pnp1 =2.0*Xpoly*Pn - Pnm1;
    Pnm1 = Pn;
    Pn   = Pnp1;
    }
    }; break;
    default : std::cout << " Orthogonal Polynomial Type index not Supported " << std::endl;
    }

	// Scale then shift if needed

	if((xScaleFactor   != 1.0))  {Pn = Pn.scale(xScaleFactor);}
    if((shift          != 0.0))  {Pn = Pn.shift(shift);}

    return Pn;
}

std::vector<PolyFun>  getOrthoPolyArray(long maxIndex) const
{
	// Create return array

	std::vector<PolyFun> P(maxIndex+1);

    PolyFun Xpoly(1);
    Xpoly[0] = 0.0; Xpoly[1] = 1.0;

    switch(polyType)
    {
//
// Legendre Functions :
//
// P_0      = 1
// P_1      = x
// P_(n+1) = ((2*n +1)*x*P_(n) + n*P_(n-1))/(n+1)
//
    case OrthoPoly::Legendre :
    {
    P[0].initialize(0);
    P[0][0] = 1.0;

    P[1].initialize(1);
    P[1][0] = 0.0; P[1][1] = 1.0;

    for(long i = 1; i < maxIndex; i++)
    {
    P[i+1] =((2.0*i + 1.0)*Xpoly*P[i] - i*P[i-1])/double(i+1);
    }
    }; break;
//
// Hermite polynomials :
//
// P_0      = 1
// P_1      = 2*x
// P_(n+1) = 2*x*P_(n) - 2*n*P_(n-1)
//
    case OrthoPoly::Hermite  :
    {
    P[0].initialize(0);
    P[0][0] = 1.0;

    P[1].initialize(1);
    P[1][0] = 0.0; P[1][1] = 2.0;

    for(long i = 1; i < maxIndex; i++)
    {
    P[i+1] = 2.0*Xpoly*P[i] - 2.0*i*P[i-1];
    }
    }; break;
//
// Laguerre polynomials :
//
// P_0      = 1
// P_1      = 1 - x
// P_(n+1) = (1 + 2*n - x)*P_(n) - n*n*P_(n-1)
//
    case OrthoPoly::Laguerre :
    {
    P[0].initialize(0);
    P[0][0] = 1.0;

    P[1].initialize(1);
    P[1][0] = 1.0; P[1][1] = -1.0;

    for(long i = 1; i < maxIndex; i++)
    {
    P[i+1] =(1.0 + 2.0*i)*P[i]  - Xpoly*P[i] - i*i*P[i-1];
    }
    }; break;
//
// Chebyshev polynomials :
//
// P_0      = 1
// P_1      = x
// P_(n+1) = (2*x)*P_(n) - P_(n-1)
//
    case OrthoPoly::Chebyshev :
    {
    P[0].initialize(0);
    P[0][0] = 1.0;

    P[1].initialize(1);
    P[1][0] = 0.0; P[1][1] = 1.0;

    for(long i = 1; i < maxIndex; i++)
    {
    P[i+1] =2.0*Xpoly*P[i] - P[i-1];
    }
    }; break;
    default : std::cout << " Orthogonal Polynomial Type Not Supported " << std::endl;
    }
//
// Scale then shift if needed
//
	if((xScaleFactor != 1.0))
    {
     for(long i = 0; i <= maxIndex; i++){P[i] = P[i].scale(xScaleFactor);}
    }
    if((shift != 0.0))
    {
     for(long i = 0; i <= maxIndex; i++){P[i] = P[i].shift(shift);}
    }

    return P;
}

    private :

    int         polyType;
    double      shift;
    double      xScaleFactor;
};

}; // namespace SCC

#endif


