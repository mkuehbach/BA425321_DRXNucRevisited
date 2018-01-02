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


#ifndef __DRXNUC_KERNEL_H__
#define __DRXNUC_KERNEL_H__

#include "DRXNUC_Datatypes.h"

//using namespace std;

class solver;

#define LOG_COORS			(0x00)
#define LOG_RHOS			(0x01)

class solverHdl //: public orimath
{
public:
	solverHdl();
	~solverHdl();

	vector<solver*> solver_instances;
	vector<runparm> configuration;
	vector<runparm3> configuration3;

	void ascii_io( const struct runparm* parms, const solver* thatsolver );
	void ascii_io3(const struct runparm3* parms, const solver* thatsolver);
	//MPI related
	int myRank;														//my MPI ID in the MPI_COMM_WORLD
	int nRanks;														//total number of MPI processes that work in the world

	//OpenMP related

	//performance counters
	profiler prof;

	//organization
	string simwrkdir;
};
typedef class solverHdl * solverHdlP;


class solver
{
public:
	solver();
	~solver();

	void run( const struct runparm * parms, const bool adaptive );
	void run2( const struct runparm * parms, const bool adaptive );
	void run3(const struct runparm3 * parms, const bool getfields );
	vector<physlog> results;
	vector<physlog3> results3;

	string fnprefix;

private:
	void setup_parameterization( const double mm, const double MM, const double edot, const double TT, const double rhod0, const double rhoda );
	void setup_fields( const unsigned int rho_ahead_model );
	double probe_strain( const double dt );
	double probe_recovery( const double dt );
	double probe_strain_and_recovery( void );
	void apply_strain_and_recovery(const double dt );
	double min_rho( const unsigned int istart, const unsigned int iend );
	inline unsigned int gbpos( const double x );
	inline unsigned long long gbpos2(const double x );

	world wo;
	gb f;
	rec r;
	unsigned int logpoint;

	vector<double> x;
	vector<double> rho; //split because operation on x and rho mostly independent

/////////////////////////////////////////////////////////
	void setup_parameterization3(const struct runparm3 * parameter); // const real mm, const real MM, const real edot, const real TT, const real rhod0, const real rhoda);
	void setup_fields3( const unsigned int gbx );
	real probe_strain_and_recovery3(const unsigned int gbpos);
	void apply_strain_and_recovery3(const unsigned int gbpos, const real dt);
	real min_rho3(const unsigned int istart, const unsigned int iend);
	inline unsigned int gbpos3(const real x);
	bool log_fields3( const char what );
	void write_binary_fields3( void );
	void build_fname3( const struct runparm3 * parms );



	world3 wo3;
	gb3 f3;
	rec3 r3;
	unsigned int logpoint3;
	vector<real> x3;
	vector<real> rho3; //split because operation on x and rho mostly independent
	vector<real*> rhofields3;

///////////////////////////////////////////////////////////
	//MPI related
	int myRank;														//my MPI ID in the MPI_COMM_WORLD
	int nRanks;														//total number of MPI processes that work in the world
	//solverHdlP owner;



	//OpenMP related

	/*//performance counters
	profiler prof;

	//organization
	string folder;
	string prefix;
	string suffix;*/
};
typedef class solver * solverP;

#endif
