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

#ifndef __DRXNUC_NUMERICS_H__
#define __DRXNUC_NUMERICS_H__

#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
	typedef float real;
	#define	MYEPSILON							(1.0e-6)
#else
	typedef double real;
	#define MYEPSILON							(1.0e-12)
#endif


#define NUCR0									(1.0e-6);
#define CONSTANT_AHEAD							(1)
#define GRADIENT_AHEAD							(2)
#define GRADIENT_AHEAD_AND_BEHIND				(3)

#define BACKSTRESS_NONE							(1)
#define BACKSTRESS_ACCUMULATE					(2)
#define BACKSTRESS_IMMEDIATE_RELAXATION			(3)

#define	DBL_EPSILON								(1.0e-12) //0.001 of wo.dx
#endif
