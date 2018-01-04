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

#include "DRXNUC_Kernel.h"
//using namespace std;

//parameter handed over
#define SIMID				1
#define CONTROLFILE			2


void info( void ) {
	//MK::provide definitions and other pieces of information relevant to the user
	cout << "\t\tRight-handed x,y,z Cartesian coordinate system!" << endl;
}


int main(int argc, char** argv)
{
//SETUP
	cout << "Starting up DRXNUC v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << endl;
	info();
	if ( argc < 2 ) {
		cout << "\t\tERROR::We need at least a simulation id <unsigned int> and a XML control file <*.xml>!" << endl; return 0;
	}
	Settings::SimID = std::stoul( argv[SIMID] );
	try {
		Settings::readXML(argv[CONTROLFILE]);
	}
	catch (exception& e) {
		cout << "\t\tERROR::Unable to parse control file! Details:\n" << e.what() << endl; return 0;
	}
	if ( Settings::checkUserInput() == false ) { 
		cout << "\t\tERROR::Control file settings failed the validity check!" << endl; return 0;
	}
	cout << "\t\tInput valid and simulation ID = " << Settings::SimID << endl;

//go MPI parallel with OpenMP threading capability
	int supportlevel_desired = MPI_THREAD_FUNNELED;
	int supportlevel_provided;
	MPI_Init_thread( &argc, &argv, supportlevel_desired, &supportlevel_provided);

	int nr, r;
	if ( supportlevel_provided < supportlevel_desired ) { std::cout << "\t\tERROR::Insufficient threading capabilities of the MPI library!" << endl; return 0; }
	else {
		cout << "\t\tInitialize MPI_COMM_WORLD..." << endl;
		MPI_Comm_size(MPI_COMM_WORLD, &nr);
		MPI_Comm_rank(MPI_COMM_WORLD, &r);
		cout << "\t\tMPI firing up rank " << r << " within MPI_COMM_WORLD of " << nr << endl;
	}

	//generate MPI rank solverHdl instance managing which of the parameter configuration rank 0 does
	solverHdl* hdl = NULL;
	hdl = new solverHdl;
	hdl->myRank = r;
	hdl->nRanks = nr;

	/*
	//all processes generate all parameter sets and number them consecutively	
	//generate naive MPI style nested for loop for parameter studies to distribute round robin on available ranks
	for ( double epsr = Settings::MacroStrainRateMin; epsr <= Settings::MacroStrainRateMax; epsr += Settings::MacroStrainRateIncr ) { //epsrate
		for ( double mGB = Settings::MobGBMin; mGB <= Settings::MobGBMax; mGB += Settings::MobGBIncr ) { //GB mobility loop
			for ( double mDisl = Settings::MobDisloMin; mDisl <= Settings::MobDisloMax; mDisl += Settings::MobDisloIncr ) { // Disloc
				for ( unsigned int rhocase = CONSTANT_AHEAD; rhocase <= GRADIENT_AHEAD; rhocase++ ) { //CON
					for ( double temp = Settings::TemperatureMin; temp <= Settings::TemperatureMax; temp += Settings::TemperatureIncr ) {
						//unsigned int backstress = BACKSTRESS_ACCUMULATE;
						for ( unsigned int backstress = BACKSTRESS_NONE; backstress <= BACKSTRESS_ACCUMULATE; backstress++ ) {//backstress
							//double rho_da = 1.0;
							for ( double rho_da = 1.0; rho_da <= 100.0; rho_da *= 10.0 ) { //gradient strength, i.e. 1.0, 10.0, 100.0
								double rho_d0 = 1.0e13;

								//construct parameter configuration
								hdl->configuration.push_back ( runparm(mGB, mDisl, epsr, temp, rho_d0, rho_da, backstress, rhocase) );
								hdl->solver_instances.push_back ( NULL );
							}//gradient strength
						}//backstress model
					}//Temperature
	
				}//rhocase
			}//disl loop
		}//mobility loop
	}//epsrate
	*/


	//all processes generate all parameter sets and number them consecutively
	//generate naive MPI style nested for loop for parameter studies to distribute round robin on available ranks
	//for ( double mGB = Settings::MobGBMin; mGB <= Settings::MobGBMax; mGB += Settings::MobGBIncr ) { //GB mobility loop
	//	for ( double mDisl = Settings::MobDisloMin; mDisl <= Settings::MobDisloMax; mDisl += Settings::MobDisloIncr ) { // Disloc
			for ( double temp = Settings::TemperatureMin; temp <= Settings::TemperatureMax; temp += Settings::TemperatureIncr ) {
				for ( double eps = Settings::MacroStrainRateMin; eps <= Settings::MacroStrainRateMax;
						eps += Settings::MacroStrainRateIncr) {
					double epsmult[3] = {1.0, 2.0, 5.0};
					for ( unsigned int i = 0; i < 3; ++i ) {
						double epsa = pow(10.0, eps) * epsmult[i];
						double epsb = pow(10.0, eps);

	//					double epsa = pow(10.0, Settings::MacroStrainRateMin);
	//					double epsb = epsa;

						double lgrho[5] = {1.0e12, 1.0e13, 1.0e14, 1.0e15, 1.0e16};
						for ( unsigned int small = 0; small < 5; ++small) { //was < 5
							for ( unsigned int large = 0; large < 5; ++large) {
								if (small < large) {
									//hdl->configuration3.push_back(runparm3(
									//		mGB, mDisl, pow(10.0, epsa), pow(10.0, epsa), pow(10.0, epsa) * epsb,
									//		temp, lgrho[small], lgrho[large], 0.5e6, 0.5e6));


									hdl->configuration3.push_back(runparm3(
										epsa, epsb, temp, lgrho[small], lgrho[large], 0.5e6, 0.5e6 ));

									hdl->solver_instances.push_back(NULL);
								}
							}
						} //rholoop
					} //local strain variation loop
				}//epsa
			}//Temperature
	//	}//disl loop
	//}//mobility loop

/*

	hdl->configuration3.push_back(runparm3( pow(10.0, Settings::MacroStrainRateMin), pow(10.0, Settings::MacroStrainRateMin),
		Settings::TemperatureMin, 1.0e13, 1.0e14, 0.5e6, 0.5e6));
	hdl->solver_instances.push_back(NULL);
*/


	if ( hdl->myRank == MASTER ) 
		cout << "Parameter configurations generated, starting now to solve " << hdl->configuration3.size() << " configurations collectively!" << endl;

	MPI_Barrier( MPI_COMM_WORLD );

	//utilize adaptive integration
	//bool integrator = true; //atoi(argv[3]);if ( integrator == true ) cout << "Adaptive integrator!" << endl;else cout << "Fixed error integrator!" << endl;
	for ( size_t c = 0; c < hdl->configuration3.size(); c++ ) {
		if ( ((1+c) % hdl->nRanks) == hdl->myRank ) { //I do
			double tic = MPI_Wtime();
			solver* sr = NULL;
			sr = new solver;
			hdl->solver_instances.at(c) = sr;

			//runparm theparameter = hdl->configuration.at(c);
			//sr->run2(&theparameter, integrator );

			runparm3 theparameter = hdl->configuration3.at(c);
			sr->run3( &theparameter, false ); //true if rhofields output desired

			//particular instance of a solver pointed to by sr holds results...

			hdl->ascii_io3( &theparameter, sr );

			double toc = MPI_Wtime();
			string str = "MPIRank " + std::to_string(hdl->myRank) + " processed case " + std::to_string(c) + " in " + std::to_string((toc-tic)) + " seconds" + "\n";
			std::cout << str;
		}
		//else other ranks take care of it
	}

	//destructor solverHdl delete !NULL solver instances and configurations cleared automatically!

	cout << "\t\tStopping MPI_COMM_WORLD..." << endl;
	MPI_Finalize();
	return 0;
}
