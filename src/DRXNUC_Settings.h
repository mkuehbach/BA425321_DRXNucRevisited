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

#ifndef __DRXNUC_SETTINGS_H__
#define __DRXNUC_SETTINGS_H__

#include "DRXNUC_Math.h"

//add thirdparty XML reader header library functionality by M. Kalicinski###
#include "thirdparty/RapidXML/rapidxml.hpp"

using namespace std;

enum E_ANALYSIS_MODE {
	E_DEBUG		//developer
	//E_1,
	//E_2
};


class Settings {
public:

	static E_ANALYSIS_MODE AnalysisMode;				//the work to do
	//static std::string H5DataFromFilename;
	static unsigned int SimID;							//consistent identifier of the data analysis run

	static double ShearModulus;
	static double BurgersVector;
	static double SurfaceTension;
	static double MobGBm0;
	static double MobGBH0;
	static double MobDisloM0;
	static double MobDisloH0;

	static double MobGBMin;								//GB migration activation enthalpy attentuation factor in H0*(1+-value)
	static double MobGBIncr;
	static double MobGBMax;
	static double MobDisloMin;							//Dislocation glide activation enthalpy attentuation factor in H0*(1+-value)
	static double MobDisloIncr;
	static double MobDisloMax;
	static double MacroStrainRateMin;					//homogeneous macrostrain rate 10^value 1/s
	static double MacroStrainRateIncr;
	static double MacroStrainRateMax;
	static double MicroStrainRateAMin;					//inhomogeneous microstrain grain homogenized
	static double MicroStrainRateBMin;
	static double TemperatureMin;						//degree Celsius isothermal assumed deformation temperature
	static double TemperatureIncr;
	static double TemperatureMax;

	static double RhoIntegratorAccuracy;
	static unsigned int OutputIncrTotal;

	static void readXML(string filename = "");
	static bool checkUserInput( void );
};

#endif
