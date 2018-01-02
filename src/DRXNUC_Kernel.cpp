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

//inline bool SortDistAsc(const ptmp &aa1 , const ptmp &aa2) { return aa1.d < aa2.d; }

solverHdl::solverHdl()
{
	myRank = MASTER;
	nRanks = 1;
	simwrkdir = "";
}

solverHdl::~solverHdl()
{
	for (unsigned int i = 0; i < this->solver_instances.size(); i++ ){
		if (solver_instances.at(i) != NULL ) {
			delete solver_instances.at(i);
		}
	}
	//configuration clears itself by vector template class destructor
}


void solverHdl::ascii_io(const struct runparm* parms, const solver* thatsolver )
{
	//results at all?
	
	//build speaking filename
	string fn;
	if ( parms->rhocase == CONSTANT_AHEAD )			fn = "CNST.";
	else											fn = "GRAD.";
	if ( parms->backstress == BACKSTRESS_NONE )		fn += "BKNONE."; //no backstresses, boundary sees only stored energy in dislocation field in front of it
	else											fn += "BKACCU."; //accumulating backstresses
	fn += "GB." + std::to_string((int)(parms->msc * 10.0)) + ".DI." + std::to_string((int)(parms->Msc*10.0));
	fn += ".EPS." + std::to_string((int)parms->epsrate) + ".T." + std::to_string((int)parms->temperature);
	fn += ".RDA." + std::to_string((int)parms->rho_da) + ".lgRD0.1E" + std::to_string(int(log(parms->rho_d0)/log(10.0)+DBL_EPSILON)) + ".csv";

	if ( thatsolver->results.size() < 1 ) { cout << "Lacking results for case study " << fn.c_str() << endl; return; }

	//single file per rank sequential I/O or MPI_COMM_SELF MPI I/O is implicit by control logic in DRXNUC_Main.cpp
	ofstream jpdf;
	jpdf.open( fn.c_str() );

	//jpdf << "Time;GBPosition;GBSpeed;RhoAhead;RhoBehind;RhoDiff;NetRhoChange;MaxNetRhoChangeAcc;rcritBailey;rcritBate;CapillaryFcap;Logic;dEsurf;dEdisl;dEnet\n";
	//jpdf << "s;m;m/s;m^-2;m^-2;m^-2;m^-2;m^-2;m;m;J/m^3;bool;J;J;J\n";
	jpdf << "Time;GBPosition;GBSpeed;RhoAhead;RhoBehind;RhoDiff;dEsurf;dEdisl;dEnet;Fcap+p;dx/LD\n";
	jpdf << "s;m;m/s;m^-2;m^-2;m^-2;J;J;J;fraction\n";

	for ( unsigned int r = 0; r < thatsolver->results.size(); r++ ) {
		jpdf << setprecision(8) << thatsolver->results[r].t << ";" << setprecision(8) << thatsolver->results[r].x << ";" << setprecision(8) << thatsolver->results[r].v << ";";
		jpdf << setprecision(12) << thatsolver->results[r].rhoD << ";" << setprecision(12) << thatsolver->results[r].rhoR << ";" << setprecision(12) << thatsolver->results[r].rhoDR << ";";
		jpdf << setprecision(12) << thatsolver->results[r].dEsurf << ";" << setprecision(12) << thatsolver->results[r].dEdisl << ";" << setprecision(12) << thatsolver->results[r].dEnet << ";";
		jpdf << setprecision(12) << (thatsolver->results[r].Fbal + thatsolver->results[r].p) << ";" << setprecision(6) << thatsolver->results[r].dr_max_acc << "\n";
	}

	jpdf.flush();
	jpdf.close();
}


void solverHdl::ascii_io3(const struct runparm3* parms, const solver* thatsolver)
{
	//results at all?

	//build speaking filename
    string fn = thatsolver->fnprefix + ".csv";
   	if (thatsolver->results3.size() < 1) { cout << "Lacking results for case study " << fn.c_str() << endl; return; }

	//single file per rank sequential I/O or MPI_COMM_SELF MPI I/O is implicit by control logic in DRXNUC_Main.cpp
	ofstream jpdf;
	jpdf.open(fn.c_str());

	jpdf << "Time;GBPosition;GBSpeed;DrivingForce;RhoBehind;RhoAhead;RhoDiff\n";
	jpdf << "s;m;m/s;J/m^3;m^-2;m^-2;m^-2\n";
	jpdf << setprecision(8);

	for (unsigned int r = 0; r < thatsolver->results3.size(); r++) {
		jpdf << thatsolver->results3[r].t << ";" << thatsolver->results3[r].x << ";" << thatsolver->results3[r].v << ";" << thatsolver->results3[r].p;
		jpdf << ";" << thatsolver->results3[r].rhoA << ";" << thatsolver->results3[r].rhoB << ";" << thatsolver->results3[r].rhoBA << "\n";
	}

	jpdf.flush();
	jpdf.close();
}





solver::solver()
{
	myRank = MASTER;
	nRanks = 1;
}

solver::~solver()
{
	for (size_t lp = 0; lp < rhofields3.size(); ++lp ) {
		if ( rhofields3.at(lp) != NULL ) {
			delete [] rhofields3.at(lp);
			rhofields3.at(lp) = NULL;
		}
	}
	rhofields3.clear();
}

void solver::setup_parameterization( const double mm, const double MM, const double edot, const double TT, const double rhod0, const double rhoda )
{
	wo.len = 15.1e-6;
	wo.dx = 2.0e-9;
	wo.n = (unsigned int) (wo.len/wo.dx);
	wo.i = 1;
	wo.imax = 1.0e8;
	wo.edot = edot;
	wo.T = TT;

	f.gamma = 0.8; //J/m^2
	f.x0 = NUCR0; //meter
	f.xmax = 15.0e-6; //meter
	f.x = f.x0;
	f.t = 0.0;
	f.m0 = 3.5e4; //m^4/Js
	f.Q = 360.0e3 * (1.0+mm);
	f.m = f.m0 * exp(-1.0*f.Q/(R*wo.T));

	r.rhoGB = 1.0e10;
	r.rhoReps = 1.0e10;
	r.rhoD0 = rhod0;
	r.rhoDa = rhoda;
	r.rhoDb = 0.5e6;
	r.G = 7.5 * 1.0e4 * 1.0e6;
	r.b = 0.25e-9;
	r.halfGb = 0.5*r.G*SQR(r.b);
	r.C = 10;
	r.M0= 5.3e7;
	r.Q = 400.0e3 * (1.0+MM);
	r.M = r.M0 * exp(-1.0*r.Q/(R*wo.T));
}

void solver::setup_fields( const unsigned int rho_ahead_model )
{
	//if constant == true, dislocation density in front of boundary constant as in Roberts/Ahlblom model, otherwise gradient
	x.reserve( wo.n );
	double binend = wo.dx;
	unsigned int i = 0;
	for ( i = 0; i < wo.n; ++i ) { //integration grid nanometer-spaced
		x.push_back( binend );
		binend += wo.dx;
	}

	rho.reserve( wo.n );
	for ( i = 0; i < wo.n; ++i ) {
		rho.push_back( 0.0 );
	}

	unsigned int gbx = gbpos( f.x );
	for ( i = 0; i <= gbx; ++i ) { //##MK::initial dislocation density behind the boundary
		rho[i] = r.rhoReps;
	}

	if ( rho_ahead_model == GRADIENT_AHEAD ) { //gradient
		for ( i = gbx+1; i < wo.n; ++i ) { //in front of the boundary
			rho[i] = r.rhoD0 * (1.0 + r.rhoDa*exp(-1.0*r.rhoDb*(x[i]-f.x)));
		}
	}
	else { //CONSTANT_AHEAD
		for ( i = gbx+1; i < wo.n; ++i ) {
			rho[i] = r.rhoD0 * (1.0 + r.rhoDa*1.0 ); //max of gradient when exp(-0)
		}
	}
}

inline unsigned int solver::gbpos( const double xpos )
{
	//numerical stable integer boundary position
	double tmp = floor( xpos/wo.dx );
	unsigned int pos = tmp; // - DBL_EPSILON; //1000.0 --> 999, 1000.5 --> 1000.0 --> 999, 1001.0, 1000
	return pos;
}

inline unsigned long long solver::gbpos2( const double xpos )
{
	//numerical stable integer boundary position
	double tmp = xpos/wo.dx;
	unsigned long long pos = tmp;
	return pos;
}

double solver::probe_strain( const double dt )
{
	double drho_max = 0.0; //strain does not reduce dislocation density so do not set to std::numeric_limits<double>:: lowest();
	double drho;
	for ( unsigned int i = 0; i < wo.n; i++ ) {
		drho = (wo.edot / (r.b * r.C)) * sqrt(rho[i]);
		if ( drho < drho_max ) 
			continue;
		//else 
		drho_max = drho;
	}
	return drho_max;
}

double solver::probe_recovery( const double dt )
{
	double drho_min = 0.0; //recovery contribution is negative
	double drho;
	for ( unsigned int i = 0; i < wo.n; i++ ) {
		drho = -2.0 * r.M * r.halfGb * SQR(rho[i]);
		if ( drho > drho_min ) 
			continue;
		//else 
		drho_min = drho;
	}
	return drho_min;
}

double solver::probe_strain_and_recovery( void )
{
	//what is currently the strongest change in dislocation density for a voxel?
	double drho_s, drho_r, drho_net;
	double drho_net_max = 0.0;
	double _ebC = wo.edot / (r.b*r.C);
	double _twoM = -2.0 * r.M * r.halfGb;
	for ( unsigned int i = 0; i < wo.n; ++i ) {
		drho_s = _ebC * sqrt(rho[i]); //strain increase multiplication with dt not required absolute change rate
		drho_r = _twoM * SQR(rho[i]); //recovery decrease
		drho_net = fabs(drho_s + drho_r);
		if ( drho_net < drho_net_max ) //most likely
			continue;
		//implicit else
		drho_net_max = drho_net;
	}
	return drho_net_max;
}

void solver::apply_strain_and_recovery( const double dt )
{
	//MK::mind order of integration!
	double drho_s, drho_r;
	double _ebC = dt * wo.edot / (r.b*r.C);
	double _twoM = dt * -2.0 * r.M * r.halfGb;
	for ( unsigned int i = 0; i < wo.n; ++i ) {
		drho_s = _ebC * sqrt(rho[i]);
		drho_r = _twoM * SQR(rho[i]);
		rho[i] = rho[i] + drho_s + drho_r;
	}
}

double solver::min_rho( const unsigned int istart, const unsigned int iend)
{
	//determine minimum of rho on interval [istart,iend]
	double minval = std::numeric_limits<double>::max();
	for (unsigned int i = istart; i <= iend; ++i) {
		if (rho[i] > minval) //most likely finding not the minimum
			continue;
		//implicit else
		minval = rho[i];
	}
	return minval;
}

void solver::run( const struct runparm * parms, const bool adaptive )
{
	//if adaptive integration error is 0.01*min(rho)
	/*folder = "case4/";
	if ( parms->rhocase == CONSTANT_AHEAD )			prefix = "CONST.";
	else if ( parms->rhocase == GRADIENT_AHEAD )	prefix = "GRADI.";
	else											return;
	suffix = "TEST";*/

	setup_parameterization( parms->msc, parms->Msc, (1.0*pow(10.0, parms->epsrate)), parms->temperature, parms->rho_d0, parms->rho_da );
	setup_fields( parms->rhocase );

	logpoint = 1;
	unsigned int lp = 0;
	bool stillgood = true;

	unsigned int subcycle_cnt = 0;
	unsigned int gbx_now = gbpos(f.x);
	unsigned int gbx_max = gbpos(f.xmax);
	if ( f.xmax >= x.at(x.size()-1) ) //maximum migration too close to right/max 1d boundary
		stillgood = false;

	//starting position of boundary to close to left/min 1d boundary
	while ( stillgood == true && gbx_now < gbx_max && wo.i < wo.imax ) { //simplified header conditions assuming f.xma
		double rhoR = rho[gbx_now]; //MK::fetch earlier entry first for caching issues
		double rhoD = rho[gbx_now+1];

		if ( (rhoD - rhoR) > 0.0 ) { //sufficient driving force?
			double p = r.halfGb * (rhoD - rhoR);
			double v = f.m * p;
			double dx = wo.dx; //definition boundary migrates an entire voxel
			double dtt = dx/v;

			//smooth integration of dislocation density field values via subcycling
			double dr_max_now = dtt * probe_strain_and_recovery();
			//INTEGRATOR ACCURACY either fixed error or adjustable error
			//double curr_minimum_rho = min_rho( 0, rho.size()-1 );
			//double dr_max_acc = (adaptive == true) ? (0.001 * curr_minimum_rho) : (DR_MAX_ACC);
			//double dr_max_acc = DR_MAX_ACC;
			//if ( adaptive == true && curr_minimum_rho > dr_max_acc ) 
			//	dr_max_acc = DR_MAX_FRACTION * curr_minimum_rho;
			//double dr_max_acc = 0.01 * curr_minimum_rho;

			double dr_max_acc = 0.001 * (rhoD - rhoR); //maximum allowed fraction of driving force variation!

			if ( dr_max_now <= dr_max_acc ) { //no subcycling apply total rho change for entire time dtt directly in which time the boundary advances by one cell
				apply_strain_and_recovery( dtt );
			}
			else { // too strong net change in dislocation density requires subcycling dr > DR_MAX_ACC
				//subcycling of recovery and straining
				double dt_subcycle_total = dtt;
				double dt_subcycle_incrm = dt_subcycle_total * (dr_max_acc / dr_max_now); //MK::mind that dr > DR_MAX_ACC !
				subcycle_cnt = 0;
				double dr_max_sub = 0.0;
				while ( (dt_subcycle_total - dt_subcycle_incrm) > DBL_EPSILON ) { //unable to integrate dislocation density, 
					apply_strain_and_recovery( dt_subcycle_incrm );
					dt_subcycle_total = dt_subcycle_total - dt_subcycle_incrm; //less time to subcycle remaining
					//probe response for next subcycling loop...
					dr_max_sub = dt_subcycle_total * probe_strain_and_recovery();
					//can we jump at once over the remaining interval?
					if ( dr_max_sub <= dr_max_acc ) { //yes, then finishing subcycling
						apply_strain_and_recovery( dt_subcycle_total );
						break;//now still _total > _incrm so we jump to the while head and next increment where the remaining time
					}
					//not broken out, so continue to subcycle
					dt_subcycle_incrm = dt_subcycle_total * (dr_max_acc / dr_max_sub);
					subcycle_cnt++;
				}
			} //advance definately boundary
			gbx_now++;
			f.x += dx;

			//std::cout << "SubcycleCnt/dr_max_acc = " << subcycle_cnt << "\t\t" << setprecision(32) << dr_max_acc << std::endl;

			//backstress model
			if ( parms->backstress == BACKSTRESS_NONE ) {
				rho[gbx_now] = r.rhoReps;
			}
			else if ( parms->backstress == BACKSTRESS_ACCUMULATE ) {
				rho[gbx_now] = min_rho( 0, (gbx_now-1) ); //MK::-1 because migration increment was already accounted for
			}
			else {
				stillgood = false;
				std::cout << "ERROR::Unknown backstress model!" << std::endl;
			}
			//real time advanced now for sure by increment dtt
			f.t += dtt;
			wo.i++;

			lp = (unsigned int) ((f.x-f.x0)/(f.xmax-f.x0)*750.0);
			if ( lp > logpoint-1) {
				double rcBAI = 2.0*f.gamma / ( r.halfGb * (rhoD - rhoR) );
				double rcBAT = 2.0*f.gamma / ( r.halfGb * rhoD * sqrt(1.0-(rhoR/rhoD)) );
				double Fcap = -1.0*f.gamma * 1.0 / (1.0 * (f.x-dx));//approximate capillary radius
				double dEsurf = +1.0*( 8.0*PI*(f.x-dx)*dx*f.gamma );
				double dEdisl = -1.0*( 4.0/3.0*PI*(3.0*SQR(f.x-dx))*p ); //..to first order ... + 3.0*(f.x-dx)*SQR(dx)

				results.push_back( physlog( f.t, dtt, f.x, v, p, rhoD, rhoR, (rhoD-rhoR), dx, dr_max_acc, rcBAI, rcBAT, Fcap, (Fcap < p), dEsurf, dEdisl, (dEsurf+dEdisl) ) );
				//std::cout << "GBPos/rhoD/rhoR/subcycle_cnt = " << f.x << "\t\t" << rhoD << "\t\t" << rhoR << "\t\t\t" << subcycle_cnt << std::endl;
				logpoint++;
			}
		}
		else {
			stillgood = false;
		}
	} //next loop

	/*
	unsigned int gbx_old = gbpos(f.x); //initial boundary position, for instance (unsigned int) 1000.0 becomes 999
	unsigned int gbx_new = gbpos(f.x);
	if ( f.xmax >= x.at(x.size()-1) ) //maximum migration too close to right/max 1d boundary
		stillgood = false;
	//starting position of boundary to close to left/min 1d boundary

	while ( stillgood == true && f.x < f.xmax && wo.i < wo.imax ) { //simplified header conditions assuming f.xmax < 
		gbx_old = gbpos(f.x);

		double rhoR = rho[gbx_old]; //fetch earlier entry first for caching issues
		double rhoD = rho[gbx_old+1];
		if ( (rhoD - rhoR) > 0.0 ) { //good
			double rcBAI = 2.0*f.gamma / (r.halfGb * (rhoD - rhoR));
			double rcBAT = 2.0*f.gamma / (r.halfGb * rhoD * sqrt(1.0-(rhoR/rhoD)));
			double Fcap = -1.0*f.gamma * 1.0 / (1.0 * f.x);

			double p = r.halfGb * (rhoD - rhoR);
			double v = f.m * p;
			double dx = CHI * wo.dx;
			double dtt = dx/v;

			double dEsurf =+1.0*(8.0*PI*f.x*dx*f.gamma);
			double dEdisl = -1.0*(4.0/3.0*PI* (3.0*SQR(f.x) + 3.0*f.x*SQR(dx)) * p);

			double drdt_def_max = probe_strain( dtt );
			double drdt_rec_max = probe_recovery( dtt );
			double dr = fabs(drdt_def_max + drdt_rec_max) * dtt;
			double dtt_reduced = dtt;
			if ( dr > DR_MAX_ACC ) {
				dtt_reduced = dtt * (DR_MAX_ACC / dr);
				dx = v*dtt_reduced;
				dtt = dtt_reduced;
			}
			f.x += dx;

			//new integer boundary position, did boundary swept the cell?
			gbx_new = gbpos(f.x);

			apply_strain_and_recovery( dtt );

			//backstress model
			if ( gbx_new > gbx_old ) { //boundary migrated into next cell
				if ( parms->backstress == BACKSTRESS_NONE ) {
					rho[gbx_old] = r.rhoReps;
				}
				else if ( parms->backstress == BACKSTRESS_ACCUMULATE ) {
					rho[gbx_old] = min_rho( 0, (gbx_old-1) );
				}
				else {
					stillgood = false;
					cout << "ERROR::Unknown backstress model!" << endl;
				}
			}

			f.t += dtt;
			wo.i++;

			double tmp = (f.x-f.x0)/(f.xmax-f.x0);
			unsigned int lp = tmp; lp *= 750;
			if ( lp > logpoint-1) {
				//##MK::log
				logpoint++;
			}
		}
		else {
			stillgood = false;
		}
	} //next loop
	*/

	//cout << "Boundary at " << (unsigned int) (f.x/wo.dx) << " after " << f.t << " seconds" << endl;
	//I/O
}


#define MINUS_FOUR_PI			((-4.0)*(PI))
#define PLUS_EIGHT_PI			((+8.0)*(PI))
#define INTEGRATOR_ACCURACY		(0.01) //the accepted dislocation density change per migration step unless the boundary is quicker then a fraction of this increment is applied
//this fixed limit error for the strong strain hardening regions but chokes integration speed when stronlgy strain harden 
#define ACCURACY_BKACCU(RHOR,RHOD)		(1.0e12) //MK:: was 1.0e12 (INTEGRATOR_ACCURACY)*(RHOR)) //mentality assure that the cells with lowest density, who strain hardenen in this simplified law strongest do not increase rho in a migration step by more than
#define ACCURACY_BKNONE(RHOR,RHOD)		(1.0e12)

#define ACCURACY_RELAXATION				(1.0e12)

void solver::run2( const struct runparm * parms, const bool adaptive )
{
	setup_parameterization( parms->msc, parms->Msc, (1.0*pow(10.0, parms->epsrate)), parms->temperature, parms->rho_d0, parms->rho_da );
	setup_fields( parms->rhocase );

	if ( f.xmax >= x.at(x.size()-1) ) { return; } //maximum migration too close to right/max 1d boundary
	logpoint = 1;
	unsigned int lp = 0;

	unsigned long long gbx_now = gbpos2(f.x);
	unsigned long long gbx_next = gbx_now;
	unsigned long long gbx_max = gbpos2(f.xmax); //starting position of boundary to close to left/min 1d boundary
	double rhoR, rhoD;
	//get first approximation of a suitable time step increment to advance boundary not farther than it takes the dislocation density behnd to change net by at most 0.01
	rhoR = rho[gbx_now];
	rhoD = rho[gbx_now+1];
	if ( (rhoD-rhoR) < DBL_EPSILON ) { cout << "Lacking sufficient driving force!" << endl; return; }
	double dt;
	if ( parms->backstress == BACKSTRESS_NONE )
		dt = ACCURACY_BKNONE(rhoR, rhoD) / probe_strain_and_recovery(); //INTEGRATOR_ACCURACY*rhoR
	else
		dt = ACCURACY_BKACCU(rhoR, rhoD) / probe_strain_and_recovery();
	
	double p = r.halfGb*(rhoD-rhoR);
	double v = f.m * p;
	double dx = v*dt;
	if ( dx > wo.dx ) {
		dt = wo.dx/v;
		dx = wo.dx;
	} //limit migration jump to +1 cell at most

	//get first position
	results.push_back( physlog( f.t, dt, f.x, v, p, rhoD, rhoR, (rhoD-rhoR), dx, (dx/wo.dx), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ) );
			

	while ( gbx_now < gbx_max && wo.i < wo.imax ) { //forward integration 1d grain boundary migration with simultaneous recovery and straining ##MK::simplified header conditions assuming f.xma
		apply_strain_and_recovery( dt );

		f.x += dx;
		f.t += dt;

		//log output?
		lp = (unsigned int) ((f.x-f.x0)/(f.xmax-f.x0)*750.0);
		if ( lp > logpoint-1) {
			double rcBAI = +2.0*f.gamma / ( r.halfGb*(rhoD-rhoR) );
			double rcBAT = +2.0*f.gamma / ( r.halfGb*rhoD*sqrt(1.0-(rhoR/rhoD)) );
			double Fcap = -1.0*f.gamma * 1.0 / (1.0*f.x); //approximate capillary radius
			double dEsurf = PLUS_EIGHT_PI * f.x * dx * f.gamma;
			double dEdisl = MINUS_FOUR_PI * SQR(f.x) * dx * p; //..to first order ... + 3.0*(f.x-dx)*SQR(dx)

			//##MK::dr_max_acc currently overwritten with 0.0
			results.push_back( physlog( f.t, dt, f.x, v, p, rhoD, rhoR, (rhoD-rhoR), dx, (dx/wo.dx), rcBAI, rcBAT, Fcap, (Fcap < p), dEsurf, dEdisl, (dEsurf+dEdisl) ) );
			//std::cout << "GBPos/rhoD/rhoR/subcycle_cnt = " << f.x << "\t\t" << rhoD << "\t\t" << rhoR << "\t\t\t" << subcycle_cnt << std::endl;
			logpoint++;
		}

		//cout << setprecision(12) << dx/wo.dx << "\t\t" << setprecision(12) << f.x <<  endl; //"dx/wo.dx = " 
		//prepare new timestep, potentially boundary migrating over domain boundary...
		gbx_next = gbpos2(f.x);

		if ( gbx_next > gbx_now ) { //entering next finite voxel..., apply backstress model!
			if ( parms->backstress == BACKSTRESS_NONE ) { //dislocation density in swept volume immediately behind is lowest
				rho[gbx_now] = r.rhoReps;
			}
			else { //parms->backstress == BACKSTRESS_ACCUMULATE{
				rho[gbx_now] = min_rho( 0, gbx_now ); //MK::-1 because migration increment was already accounted for
			}
			rhoR = rho[gbx_now];
			rhoD = rho[gbx_now+1];
		}
		else { //remaining in the same voxel, update working dislocation densities the boundary sees
			//evaluate rho
			rhoR = rho[gbx_now-1];
			rhoD = rho[gbx_now];
		}

		if ( (rhoD-rhoR) > DBL_EPSILON ) { 
			gbx_now = gbx_next;

			if ( parms->backstress == BACKSTRESS_NONE )
				dt = ACCURACY_BKNONE(rhoR, rhoD) / probe_strain_and_recovery(); //INTEGRATOR_ACCURACY*rhoR
			else
				dt = ACCURACY_BKACCU(rhoR, rhoD) / probe_strain_and_recovery();

			p = r.halfGb*(rhoD-rhoR);
			v = f.m*p;
			dx = v*dt;
			if ( dx > wo.dx ) {
				dt = wo.dx/v;
				dx = wo.dx;
			} //limit migration jump to +1 cell at most

			wo.i++; //prepare for next iteration
			//std::cout << wo.i << "\t\t" << f.x << std::endl;
			continue;
		}
		else {
			cout << "Lacking sufficient driving force!" << endl; return;
		}
	} //next loop
}





void solver::setup_parameterization3( const struct runparm3 * parameter )
{
	wo3.len = 15.1e-6;
	wo3.dx = 2.0e-9;
	wo3.edot = 0.0; //parameter->epsrate;
	wo3.edotA = parameter->aepsrateA;
	wo3.edotB = parameter->aepsrateB;
	wo3.T = parameter->temperature;
	wo3.n = static_cast<unsigned int>(wo3.len / wo3.dx);
	wo3.i = 1;
	wo3.imax = static_cast<unsigned int>(1.0e8);
	

	f3.x0 = NUCR0; //meter
	f3.xmax = 15.0e-6; //meter
	f3.x = f3.x0;
	f3.t = 0.0;
	//because of significant uincertainty in mobility values, we scan order of magnitude differences
	//f3.m0 = 3.5e4; //m^4/Js
	//f3.Q = 360.0e3 * (1.0 + parameter->msc);
	//f3.m = f3.m0 * exp(-1.0*f3.Q / (R*wo3.T));

	f3.gamma = Settings::SurfaceTension; //J/m^2
	f3.m0 = Settings::MobGBm0; //m^4/Js
	f3.Q = Settings::MobGBH0; // in eV
	f3.m = f3.m0 * exp(-1.0*f3.Q / (JOULE2EV(kboltzmann)*wo3.T));


	r3.rhoGB = 1.0e10;
	r3.rhoA0 = parameter->rho_A0;
	r3.rhoB0 = parameter->rho_B0;
	r3.cA0 = parameter->rho_cA;
	r3.cB0 = parameter->rho_cB;
	r3.CHoltA = 10;
	r3.CHoltB = 10;
	//r3.G = 7.5 * 1.0e4 * 1.0e6;
	//r3.b = 0.25e-9;
	//r3.halfGb = 0.5* r3.G * SQR(r3.b);
	r3.G = Settings::ShearModulus;
	r3.b = Settings::BurgersVector;
	r3.halfGb = 0.5 * r3.G * SQR(r3.b);

	//r3.M0 = 5.3e7;
	//r3.Q = 400.0e3 * (1.0 + parameter->Msc);
	//r3.M = r3.M0 * exp(-1.0*r3.Q / (R*wo3.T));

	r3.M0 = Settings::MobDisloM0; //m^3/Js
	r3.Q = Settings::MobDisloH0; //in eV
	r3.M = r3.M0 * exp(-1.0*r3.Q / (JOULE2EV(kboltzmann)*wo3.T) );

	cout << "LineEnergy\t\t\t\t" << r3.halfGb << endl;
	cout << "BoundaryMob\t\t\t\t" << f3.m << endl;
	cout << "BoundaryGamma\t\t\t\t" << f3.gamma << endl;
	cout << "DisloMobility\t\t\t\t" << r3.M << endl;

}

void solver::setup_fields3( const unsigned int gbx )
{
	//MK::implements the GRADIENT_AHEAD_AND_BEHIND and BACKPRESSURE_RELAXING model

	//init one-dimensional world
	x3.reserve(wo3.n);
	real binend = wo3.dx;
	for (unsigned int i = 0; i < wo3.n; ++i) { //integration grid nanometer-spaced
		x3.push_back(binend);
		binend += wo3.dx;
	}

	//init dislocation density field
	rho3.reserve(wo3.n);
	for (unsigned int i = 0; i < wo3.n; ++i) {
		rho3.push_back(0.0);
	}

	for (unsigned int i = 0; i <= gbx; ++i) { //##MK::initial dislocation density behind the boundary
		unsigned int distance = gbx - i;
		real d = static_cast<real>(distance);
		d *= wo3.dx;
		rho3[i] = r3.rhoA0 * (1+exp(-1.0*r3.cA0 * d));
	}
	for (unsigned int i = gbx+1; i < wo3.n; ++i ) { //in front of the boundary
		unsigned int distance = i - gbx;
		real d = static_cast<real>(distance);
		d *= wo3.dx;
		rho3[i] = r3.rhoB0 * (1+exp(-1.0*r3.cB0 * d));
	}

}

real solver::probe_strain_and_recovery3( const unsigned int gbpos )
{
	//what is currently the strongest change in dislocation density for a voxel?
	real drho_s, drho_r, drho_net; //s - strain, r - recovery, net - ... guess what ...
	real drho_net_max = 0.0;

	real _ebCA = wo3.edotA / (r3.b * r3.CHoltA);
	real _ebCB = wo3.edotB / (r3.b * r3.CHoltB);

	real _twoM = -2.0 * r3.M * r3.halfGb;

	//recovery part A of bicrystal
	for (unsigned int i = 0; i <= gbpos; ++i) {
		drho_s = _ebCA * sqrt(rho3[i]); //strain increase multiplication with dt not required absolute change rate
		drho_r = _twoM * SQR(rho3[i]); //recovery decrease
		drho_net = fabs(drho_s + drho_r);
		if (drho_net < drho_net_max) //most likely
			continue;
		//implicit else
		drho_net_max = drho_net;
	}

	//recovery part B of bicrystal
	for (unsigned int i = gbpos + 1; i < wo3.n; ++i) {
		drho_s = _ebCB * sqrt(rho3[i]); //strain increase multiplication with dt not required absolute change rate
		drho_r = _twoM * SQR(rho3[i]); //recovery decrease
		drho_net = fabs(drho_s + drho_r);
		if (drho_net < drho_net_max) //most likely
			continue;
		//implicit else
		drho_net_max = drho_net;
	}

	//conservative forward integration, maximum increment change determines time integration
	return drho_net_max;
}

void solver::apply_strain_and_recovery3(const unsigned int gbpos, const real dt)
{
	//MK::mind order of integration!
	real drho_s, drho_r;
	real _ebCA = dt * wo3.edotA / (r3.b*r3.CHoltA);
	real _ebCB = dt * wo3.edotB / (r3.b*r3.CHoltB);
	real _twoM = dt * -2.0 * r3.M * r3.halfGb;

	for (unsigned int i = 0; i <= gbpos; ++i) {
		drho_s = _ebCA * sqrt(rho3[i]);
		drho_r = _twoM * SQR(rho3[i]);
		rho3[i] = rho3[i] + drho_s + drho_r;
	}

	for (unsigned int i = gbpos+1; i < wo3.n; ++i) {
		drho_s = _ebCB * sqrt(rho3[i]);
		drho_r = _twoM * SQR(rho3[i]);
		rho3[i] = rho3[i] + drho_s + drho_r;
	}
}

real solver::min_rho3(const unsigned int istart, const unsigned int iend)
{
	//determine minimum of rho on interval [istart,iend]
	real minval = std::numeric_limits<real>::max();
	for (unsigned int i = istart; i <= iend; ++i) {
		if (rho3[i] > minval) //most likely finding not the minimum
			continue;
		//implicit else
		minval = rho3[i];
	}
	return minval;
}

inline unsigned int solver::gbpos3(const real xpos)
{
	//numerical stable integer boundary position
	real tmp = floor(xpos / wo3.dx);
	unsigned int pos = static_cast<unsigned int>(tmp); // - DBL_EPSILON; //1000.0 --> 999, 1000.5 --> 1000.0 --> 999, 1001.0, 1000
	return pos;
}

void solver::run3(const struct runparm3 * parms, const bool getfields )
{
	bool status = true;

	build_fname3( parms );

	setup_parameterization3( parms );
	//cout << "Parameter setup!" << endl;

	logpoint = 1;
	unsigned int lp = 0;

	unsigned int gbx_now = gbpos3(f3.x);
	unsigned int gbx_max = gbpos3(f3.xmax);

	setup_fields3( gbx_now );
	
	if ( x3.size() != rho3.size() ) {
		cout << "Field discretization is inconsistent!" << endl;
		return;
	}

	//cout << "Fields setup!" << endl;

	if (f3.xmax >= x3.back() ) { //at(x3.size() - 1)) { 
		cout << "Maximum migration too close to right/max 1d boundary!" << endl;
		return;
	}

	real rhoA = rho3[gbx_now]; //current rho in part A of bicrystal aka "behind the boundary"
	real rhoB = rho3[gbx_now+1]; //current rho in part B of bicrystal aka "in front of the boundary"

	if ((rhoB - rhoA) < MYEPSILON) { 
		cout << "Lacking sufficient driving force!" << endl;
		cout << "Time/GBPosition/rhoB/rhoA\t\t" << f3.t << "\t\t" << f3.x << "\t\t" << rhoB << "\t\t" << rhoA << endl;

		if ( getfields == true ) {
			write_binary_fields3();
		}

		return;
	}

	if ( getfields == true )
		status = log_fields3( LOG_COORS );

	real acc = Settings::RhoIntegratorAccuracy * min_rho3( 0, gbx_max );
	real dt = acc / probe_strain_and_recovery3( gbx_now );
	//real dt = Settings::RhoIntegratorAccuracy / probe_strain_and_recovery3( gbx_now );
	
	real p = r3.halfGb*(rhoB - rhoA); //BACKSTRESS_RELAXATION MODEL 
	real v = f3.m * p;
	real dx = v * dt;

	if (dx > wo3.dx) { //limit migration jump to +1 cell at most
		dt = wo3.dx / v;
		dx = wo3.dx;
	} 

	//log first position
	results3.push_back(physlog3(f3.t, f3.x, v, p, rhoB, rhoA, (rhoB-rhoA)) );
	
	//forward integration 1d grain boundary migration with simultaneous recovery and straining ##MK::simplified header conditions assuming f.xma
	while (gbx_now < gbx_max && wo3.i < wo3.imax) { 
		
		apply_strain_and_recovery3(gbx_now, dt);

		f3.x += dx;
		f3.t += dt;

		//log output? partitioning total migration distance into 1000.0 equal increments
		real tmp = (f3.x - f3.x0) / (f3.xmax - f3.x0) * static_cast<real>(Settings::OutputIncrTotal);
		lp = static_cast<unsigned int>(tmp);

		if ( lp > (logpoint-1) ) {

			/*real rcBAI = +2.0*f3.gamma / (r3.halfGb*(rhoB - rhoA));
			real rcBAT = +2.0*f3.gamma / (r3.halfGb*rhoB*sqrt(1.0-(rhoA/rhoB)));
			real Fcap = -1.0*f3.gamma * 1.0 / (1.0*f3.x); //approximate capillary radius
			real dEsurf = PLUS_EIGHT_PI * f3.x * dx * f3.gamma;
			real dEdisl = MINUS_FOUR_PI * SQR(f3.x) * dx * p; //..to first order ... + 3.0*(f.x-dx)*SQR(dx)*/

			//##MK::dr_max_acc currently overwritten with 0.0
			results3.push_back( physlog3(f3.t, f3.x, v, p, rhoB, rhoA, (rhoB - rhoA)) ); /*rcBAI, rcBAT, Fcap, (Fcap < p), dEsurf, dEdisl, (dEsurf + dEdisl)));*/
			//cout << "Time/GBPosition/rhoB/rhoA\t\t" << f3.t << "\t\t" << f3.x << "\t\t" << rhoB << "\t\t" << rhoA << endl;

			if ( getfields == true )
				status = log_fields3( LOG_RHOS );

			logpoint++;
		}

		//cout << setprecision(8) << f3.t << "\t\t" << f3.x << "\t\t" << rhoB << "\t\t" << rhoA << "\t\t" << acc << endl;
		//prepare new timestep, potentially boundary migrating over domain boundary...
		unsigned int gbx_next = gbpos3(f3.x);
		if ( gbx_next > gbx_now) {
			//let dislocation density remain as we assume the pipe up to relax immediately back into the boundary
			//rho[gbx_now] = rho[gbx_now - 1];

			real rho3_relaxback = rho3[gbx_now];

			//relax dislocation content
			rho3[gbx_now+1] = rho3_relaxback;

			//migration increment accounting
			gbx_now++;

			//reevaluate dislocation state entering next finite voxel..., apply dislocation density relaxation model
			rhoA = rho3[gbx_now];
			rhoB = rho3[gbx_now + 1]; //because jumping not more than length of one finite voxel
		}
		else { //remaining in the same voxel, update working dislocation densities the boundary sees
			//evaluate rho
			//rhoA = rho3[gbx_now - 1];
			//rhoB = rho3[gbx_now];

			rhoA = rho3[gbx_now];
			rhoB = rho3[gbx_now + 1];
		}

		if ((rhoB - rhoA) >= MYEPSILON) {
			gbx_now = gbpos3(f3.x);

			acc = Settings::RhoIntegratorAccuracy * min_rho3( 0, gbx_max );
			dt = acc / probe_strain_and_recovery3( gbx_now );
			//dt = Settings::RhoIntegratorAccuracy / probe_strain_and_recovery3(gbx_now);

			p = r3.halfGb*(rhoB - rhoA);
			v = f3.m * p;
			dx = v * dt;
			if (dx > wo3.dx) {
				dt = wo3.dx / v;
				dx = wo3.dx;
			} //limit migration jump to +1 cell at most
			
			//prepare for next iteration
			wo3.i++; 
		}
		else {
			cout << "Lacking sufficient driving force!" << endl;
			cout << "Time/GBPosition/rhoB/rhoA\t\t" << f3.t << "\t\t" << f3.x << "\t\t" << rhoB << "\t\t" << rhoA << endl;

			if ( getfields == true ) {
				write_binary_fields3();
			}

			return;
		}
	} //next loop


	if ( getfields == true ) {
		write_binary_fields3();
	}
}

bool solver::log_fields3( const char what )
{
	//rows along +x, columns along time steps
	unsigned int nrows = 1 + static_cast<unsigned int>(x3.size()); //first is time, rest are field values beginning with 0
	real* bucket = NULL;
	try { bucket = new real[nrows]; }
	catch (std::bad_alloc &memexc) {
		return false;
	}

	if ( what == LOG_RHOS ) {
		bucket[0] = f3.t;
		for ( unsigned int i = 0; i < rho3.size(); ++i ) {
			bucket[1+i] = rho3[i];
		}
		rhofields3.push_back( bucket ); //copyconstructor
	}
	else if ( what == LOG_COORS ) {
		bucket[0] = f3.t;
		for ( unsigned int i = 0; i < x3.size(); ++i ) {
			bucket[1+i] = x3[i];
		}
		rhofields3.push_back( bucket ); //copyconstructor
	}
	else {
		//nothing to log
	}
	return true;
}


void solver::write_binary_fields3( void )
{
	//store matrix of dislocation density fields
	size_t ncols = rhofields3.size();
	size_t nrows = 1 + x3.size();

    string fn = fnprefix + ".NR." + to_string(nrows) + ".NC." + to_string(ncols) + ".RhoFields.bin";
	ofstream myfile ( fn.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	if ( myfile.is_open() ) {
		real* rbucket = NULL;
		for ( size_t cols = 0; cols < ncols; ++cols ) {
			rbucket = rhofields3.at(cols);
			//cout << "\t\tWriting column " << cols << " of ncols = " << ncols << " nrows = " << nrows << endl;
			myfile.write( reinterpret_cast<const char*>(rbucket), std::streamsize(nrows*sizeof(real)) );
		}
		myfile.close();
	}
	else {
		cout << "Unable to write binary field data to file!" << endl;
	}
}

void solver::build_fname3( const struct runparm3 * parms )
{
	fnprefix = "RELAX.SimID." + to_string(Settings::SimID);
	//fnprefix += ".GB." + to_string((int)(parms->msc*10.0)) + ".DI." + to_string((int)(parms->Msc*10.0));
	//fnprefix += ".GB." + to_string((int)(log10(f3.m))) + ".DI." + to_string((int)(log10(r3.M)));
	fnprefix += ".EPSA." + to_string((int)(log10(parms->aepsrateA)*100.0)) + ".EPSB." + to_string((int)(log10(parms->aepsrateB)*100.0));
	fnprefix += ".T." + to_string((int)parms->temperature);
	fnprefix += ".RA." + to_string((int)log10(parms->rho_A0)) + ".RB." + to_string((int)log10(parms->rho_B0));
}

