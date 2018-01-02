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

#ifndef __DRXNUC_ABINITIO_H__
#define __DRXNUC_ABINITIO_H__

//natural constants
#define PI								(3.141592653589793238462643383279502884197169399375105820974)
#define kboltzmann						(1.3806488e-23)		//J/Kmol
#define echarge							(1.602176565e-19)	//Coulomb
#define Navogadro						(6.022140857e+23)	//1/mol
#define R								(8.31446154) //(Navogadro)*(kboltzmann)) //J/K/mol

//unit conversions
#define CELSIUS2KELVIN(T)				((273.15) + (T))
#define KELVIN2CELSIUS(T)				((T) - (273.15))

#define DEGREE2RADIANT(theta)			((theta)/(180.0)*(PI))
#define RADIANT2DEGREE(rad)				((rad)/(PI)*(180.0))

//scaling conversions
#define MICRON2METER(x)					((x)*(1.0e-6))
#define METER2MICRON(x)					((x)*(1.0e6))

#define EV2SI(EV)						((EV)*(1.602e-19))
#define JOULE2EV(J)						((J)*(6.242e+18))


//naming definitions for more intuitive code reading
#define TWO_DIMENSIONS					2
#define THREE_DIMENSIONS				3

#endif
