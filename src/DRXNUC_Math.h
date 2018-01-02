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

#ifndef __DRXNUC_MATH_H__
#define __DRXNUC_MATH_H__

#include "DRXNUC_STLIncludes.h"
#include "DRXNUC_Parallelization.h"
#include "DRXNUC_Information.h"
#include "DRXNUC_Numerics.h"
#include "DRXNUC_Abinitio.h"

using namespace std;

struct interval
{
	double min;
	double incr;
	double max;
};
typedef struct interval * intervalP;

#define SQR(a)					((a)*(a))
#define CUBE(a)					((a)*(a)*(a))
#define MIN(X,Y)				(((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y)				(((X) > (Y)) ? (X) : (Y))

#endif
