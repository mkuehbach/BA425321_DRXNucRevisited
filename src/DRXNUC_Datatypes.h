/*
	Copyright Max-Planck-Institut f\"ur Eisenforschung, GmbH, D\"usseldorf
	Programmer: Markus K\"uhbach, 2017

	drxnuc --- this program implements functionality to study the effect of
	dislocation density gradients on the 1d migration of (discontinuously) dynamically
	recrystallizing nuclei.

	The source code was developed by Markus K\"uhbach during his postdoc with the Max-Planck-Institut f\"ur
	Eisenforschung GmbH in D\"usseldorf, contact should be directed to m.kuehbach at mpie.de

	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the project "Consistent physically-based modeling of dynamic
	recrystallization under hot working conditions", project no 315419526.

	This file is part of drxnuc.

	drxnuc is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	drxnuc is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with drxnuc.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __DRXNUC_DATATYPES_H__
#define __DRXNUC_DATATYPES_H__

#include "DRXNUC_Settings.h"


struct profiler
{
	double t0;
	double tstart;
	double tend;
	profiler() : t0(0.0), tstart(0.0), tend(0.0) {}
};
typedef struct profiler * profilerP;


struct runparm
{
	double msc;
	double Msc;
	double epsrate;
	double temperature;
	double rho_d0;
	double rho_da;
	unsigned int backstress;
	unsigned int rhocase;
	runparm( double _msc, double _Msc, double _epsr, double _temp, double _rhod0, double _rhoda, 
		unsigned int _bs, unsigned int _rc ) : msc(_msc), Msc(_Msc), epsrate(_epsr), temperature(_temp),
		rho_d0(_rhod0), rho_da(_rhoda), backstress(_bs), rhocase(_rc) {}
	~runparm(){}
};
typedef struct runparm * runparmP;

struct world
{
	double len;
	double dx;
	double edot;			//globally imposed macrostrain rate
	double T;
	unsigned int n;
	unsigned int i;
	unsigned int imax;
	//##MK::world() : len(15.1e-6), dx(1.0e-9), n( (unsigned int) (len/dx) ),i(1), imax(1.0e8),edot(0.0)
};
typedef struct world * worldP;

struct gb
{
	double gamma;			//GB spec. energy	
	double x0;				//initial boundary position
	double xmax;			//maximum migrated to position along +x
	double x;				//current position
	double t;				//current time
	double m0;				//mobility GB preexp
	double Q;				//activation enthalphy GB
	double m;				//mobility GB
							//##MK::add defaults
};
typedef struct gb * gbP;

struct rec
{
	double rhoGB;
	double rhoReps;
	double rhoD0;
	double rhoDa;
	double rhoDb;
	double G;
	double b;
	double halfGb;
	double C;
	double M0;
	double Q;
	double M;
	//##MK::add defaults
};
typedef struct rec * recP;

struct conf
{
	double m;
	double M;
	double rhocase;
	double edot;
	double T;
	double bs;
	double d0;
	double da;
	double stillrolling;
	double arrival;
	double totalmacrostrain;
	//##MK::add default
	conf(double _m, double _M, double _rhocase, double _edot, double _T, double _bs,
		double _d0, double _da, double _stillrolling, double _arrival, double _totalmacrostrain) :
		m(_m), M(_M), rhocase(_rhocase), edot(_edot), T(_T), bs(_bs), d0(_d0), da(_da), stillrolling(_stillrolling),
		arrival(_arrival), totalmacrostrain(_totalmacrostrain) {}
	~conf() {}
};
typedef struct conf * confP;

struct physlog
{
	double t;
	double dt;
	double x;
	double v;
	double p;
	double rhoD;
	double rhoR;
	double rhoDR;
	double dr;
	double dr_max_acc;
	double rcritBH;
	double rcritBate;
	double Fbal;
	double RadiusConsistentToBulge;
	double dEsurf;
	double dEdisl;
	double dEnet;
	physlog(double _t, double _dt, double _x, double _v, double _p, double _rhoD, double _rhoR, double _rhoDR,
		double _dr, double _drmaxacc, double _rcritBH, double _rcritBate,
		double _Fbal, double _Rad, double _dEs, double _dEd, double _dEn) :
		t(_t), dt(_dt), x(_x), v(_v), p(_p), rhoD(_rhoD), rhoR(_rhoR), rhoDR(_rhoDR),
		dr(_dr), dr_max_acc(_drmaxacc), rcritBH(_rcritBH),
		rcritBate(_rcritBate), Fbal(_Fbal), RadiusConsistentToBulge(_Rad),
		dEsurf(_dEs), dEdisl(_dEd), dEnet(_dEn) {}
	~physlog() {}
};
typedef struct physlog * physlogP;



struct runparm3
{
	//real msc;
	//real Msc;
	/*real mm0;
	real mH0;
	real Mm0;
	real MH0;
	*/
	//real epsrate;
	real aepsrateA;
	real aepsrateB;
	real temperature;
	real rho_A0;
	real rho_B0;
	real rho_cA;
	real rho_cB;

/*
	runparm3(const real _msc, const real _Msc, const real _mm0, const real _mH0, const real _Mm0, const real _MH0,
			const real _epsr, const real _epsra, const real _epsrb,
		const real _temp, const real _rhoA0, const real _rhoB0, const real _rhocA, const real _rhocB) :
		msc(_msc), Msc(_Msc), mm0(_mm0), mH0(_mH0), Mm0(_Mm0), MH0(_MH0), epsrate(_epsr),
		aepsrateA(_epsra), aepsrateB(_epsrb), temperature(_temp), rho_A0(_rhoA0), rho_B0(_rhoB0), rho_cA(_rhocA), rho_cB(_rhocB) {}
*/
	runparm3( const real _epsra, const real _epsrb,
			const real _temp, const real _rhoA0, const real _rhoB0, const real _rhocA, const real _rhocB) :
		aepsrateA(_epsra), aepsrateB(_epsrb), temperature(_temp),
		rho_A0(_rhoA0), rho_B0(_rhoB0), rho_cA(_rhocA), rho_cB(_rhocB) {}
	~runparm3() {}
};
typedef struct runparm3 * runparm3P;


struct physlog3
{
	real t;
	real x;
	real v;
	real p;
	real rhoA;
	real rhoB;
	real rhoBA;
	physlog3(const real _t, const real _x, const real _v, const real _p, 
		const real _rhoA, const real _rhoB, const real _rhoBA ) :
		t(_t), x(_x), v(_v), p(_p), rhoA(_rhoA), rhoB(_rhoB), rhoBA(_rhoBA) {}
	~physlog3() {}
};
typedef struct physlog3 * physlog3P;

struct world3
{
	real len;
	real dx;
	real edot;			//globally imposed macrostrain rate
	real edotA;			//strain rate of the two parts of the bicrystal arrangement
	real edotB;
	real T;
	unsigned int n;
	unsigned int i;
	unsigned int imax;
	//##MK::world() : len(15.1e-6), dx(1.0e-9), n( (unsigned int) (len/dx) ),i(1), imax(1.0e8),edot(0.0)
};
typedef struct world3 * world3P;

struct gb3
{
	real gamma;			//GB spec. energy	
	real x0;				//initial boundary position
	real xmax;			//maximum migrated to position along +x
	real x;				//current position
	real t;				//current time
	real m0;				//mobility GB preexp
	real Q;				//activation enthalphy GB
	real m;				//mobility GB
};
typedef struct gb3 * gb3P;

struct rec3
{
	real rhoGB;			//dislo dens in boundary
	real rhoA0;			//dislo dens base grain A (behind part)
	real rhoB0;
	real cA0;				//scaling of initial dislocation density exponential profile close to boundary in grain A
	real cB0;
	real CHoltA;			//scaling Holt grain A
	real CHoltB;			//scaling Holt grain B

	real G;				//shear modulus
	real b;				//Burgers vector
	real halfGb;			//0.5Gb^2
	real M0;				//mobility preexp dislo
	real Q;				//activation enthalpy
	real M;				//mobility dislo
};
typedef struct rec3 * rec3P;





 
#endif
