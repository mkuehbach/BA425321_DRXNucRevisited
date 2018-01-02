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


#include "DRXNUC_Settings.h"

using namespace rapidxml;


E_ANALYSIS_MODE Settings::AnalysisMode = E_DEBUG;

unsigned int Settings::SimID = 0;

double Settings::ShearModulus = 75.0e9;
double Settings::BurgersVector = 2.58e-10;
double Settings::SurfaceTension = 0.8;
double Settings::MobGBm0 = 1.0;
double Settings::MobGBH0 = 1.0;
double Settings::MobDisloM0 = 1.0;
double Settings::MobDisloH0 = 1.0;

double Settings::MobGBMin = -0.10;								//GB migration activation enthalpy attentuation factor in H0*(1+-value)
double Settings::MobGBIncr = 0.10;
double Settings::MobGBMax = +0.10;

double Settings::MobDisloMin = -0.10;							//Dislocation glide activation enthalpy attentuation factor in H0*(1+-value)
double Settings::MobDisloIncr = 0.10;
double Settings::MobDisloMax = +0.10;

double Settings::MacroStrainRateMin = -4.0;					//homogeneous macrostrain rate 10^value 1/s
double Settings::MacroStrainRateIncr = 1.0;
double Settings::MacroStrainRateMax = +2.0;
double Settings::MicroStrainRateAMin = -3.0;					//inhomogeneous microstrain grain homogenized
double Settings::MicroStrainRateBMin = -3.0;
double Settings::TemperatureMin = 800.0;						//degree Celsius isothermal assumed deformation temperature
double Settings::TemperatureIncr = 100.0;
double Settings::TemperatureMax = 1200.0;
//double Settings::RhoIntegratorAccuracy = 1.0e11;
double Settings::RhoIntegratorAccuracy = 1.0e-4;
unsigned int Settings::OutputIncrTotal = 1000;
//string Settings::TargetGIDFromFilename;


void Settings::readXML(string filename) {

	//find the desired .xml file
	if ( 0 == filename.compare("") )
		filename = string("DRXNUC.Input.xml");

	ifstream file( filename );
	if (file.fail()) {
		throw runtime_error(string("Unable to locate input file ") + filename);
	}

	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("Undefined parameters file!");
	}

	unsigned int mode = 0;
	if (0 != rootNode->first_node("AnalysisMode")) {
		mode = std::stoul( rootNode->first_node("AnalysisMode")->value());
	}
	if ( mode == E_DEBUG )
		AnalysisMode = E_DEBUG;

	if (0 != rootNode->first_node("ShearModulus"))
		ShearModulus = std::stod( rootNode->first_node("ShearModulus")->value() );
	if (0 != rootNode->first_node("BurgersVector"))
		BurgersVector = std::stod( rootNode->first_node("BurgersVector")->value() );
	if (0 != rootNode->first_node("SurfaceTension"))
		SurfaceTension = std::stod( rootNode->first_node("SurfaceTension")->value() );

	if (0 != rootNode->first_node("MobGBm0"))
		MobGBm0 = std::stod( rootNode->first_node("MobGBm0")->value() );
	if (0 != rootNode->first_node("MobGBH0"))
		MobGBH0 = std::stod( rootNode->first_node("MobGBH0")->value() );
	if (0 != rootNode->first_node("MobDisloM0"))
		MobDisloM0 = std::stod( rootNode->first_node("MobDisloM0")->value() );
	if (0 != rootNode->first_node("MobDisloH0"))
		MobDisloH0 = std::stod( rootNode->first_node("MobDisloH0")->value() );

	if (0 != rootNode->first_node("MobGBMin")) 
		MobGBMin = std::stod( rootNode->first_node("MobGBMin")->value() );
	if (0 != rootNode->first_node("MobGBIncr")) 
		MobGBIncr = std::stod( rootNode->first_node("MobGBIncr")->value() );
	if (0 != rootNode->first_node("MobGBMax")) 
		MobGBMax = std::stod( rootNode->first_node("MobGBMax")->value() );

	if (0 != rootNode->first_node("MobDisloMin")) 
		MobDisloMin = std::stod( rootNode->first_node("MobDisloMin")->value() );
	if (0 != rootNode->first_node("MobDisloIncr")) 
		MobDisloIncr = std::stod( rootNode->first_node("MobDisloIncr")->value() );
	if (0 != rootNode->first_node("MobDisloMax")) 
		MobDisloMax = std::stod( rootNode->first_node("MobDisloMax")->value() );

	if (0 != rootNode->first_node("MacroStrainRateMin")) 
		MacroStrainRateMin = std::stod( rootNode->first_node("MacroStrainRateMin")->value() );
	if (0 != rootNode->first_node("MacroStrainRateIncr")) 
		MacroStrainRateIncr = std::stod( rootNode->first_node("MacroStrainRateIncr")->value() );
	if (0 != rootNode->first_node("MacroStrainRateMax")) 
		MacroStrainRateMax = std::stod( rootNode->first_node("MacroStrainRateMax")->value() );

	if (0 != rootNode->first_node("MicroStrainRateAMin"))
		MicroStrainRateAMin = std::stod(rootNode->first_node("MicroStrainRateAMin")->value());
	if (0 != rootNode->first_node("MicroStrainRateBMin"))
		MicroStrainRateBMin = std::stod(rootNode->first_node("MicroStrainRateBMin")->value());
	
	if (0 != rootNode->first_node("TemperatureMin")) 
		TemperatureMin = CELSIUS2KELVIN(std::stod( rootNode->first_node("TemperatureMin")->value() ));
	if (0 != rootNode->first_node("TemperatureIncr")) 
		TemperatureIncr = std::stod( rootNode->first_node("TemperatureIncr")->value() );
	if (0 != rootNode->first_node("TemperatureMax")) 
		TemperatureMax = CELSIUS2KELVIN(std::stod( rootNode->first_node("TemperatureMax")->value() ));

	if (0 != rootNode->first_node("RhoIntegratorAccuracy"))
		RhoIntegratorAccuracy = std::stod( rootNode->first_node("RhoIntegratorAccuracy")->value() );
	if (0 != rootNode->first_node("OutputIncrementsTotal"))
			OutputIncrTotal = std::stol( rootNode->first_node("OutputIncrementsTotal")->value() );


	//convert units to SI
}


bool Settings::checkUserInput( void )
{
	//MK::check user input for validity and good sense
	//default mode is E_DEBUG
	//no other checks at the moment
	//if ( Settings::RhoIntegratorAccuracy < 1.0e11 ) {
	//	Settings::RhoIntegratorAccuracy = 1.0e11;
	//}
	//if ( Settings::RhoIntegratorAccuracy > 1.0e13 ) {
	//	Settings::RhoIntegratorAccuracy = 1.0e13;
	//}
	if ( Settings::OutputIncrTotal < 0 || Settings::OutputIncrTotal > 15000 ) {
		OutputIncrTotal = 14000;
	}
	cout << "SimulationID\t\t\t\t" << Settings::SimID << endl;
	cout << "ShearModulus\t\t\t\t" << Settings::ShearModulus << endl;
	cout << "BurgersVector\t\t\t\t" << Settings::BurgersVector << endl;

	cout << "MobilityGSm0\t\t\t\t" << Settings::MobGBm0 << endl;
	cout << "MobilityGSH0\t\t\t\t" << Settings::MobGBH0 << endl;
	cout << "MobilityDislom0\t\t\t\t" << Settings::MobDisloM0 << endl;
	cout << "MobilityDisloH0\t\t\t\t" << Settings::MobDisloH0 << endl;

	cout << "RhoIntegratorAccuracy\t\t\t\t" << Settings::RhoIntegratorAccuracy << endl;
	cout << "OutputIncrementsTotal\t\t\t\t" << Settings::OutputIncrTotal << endl;
	return true;
}

