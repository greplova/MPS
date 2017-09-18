#include "core.h"
#include "sites/spinhalf.h"
#include "autompo.h"
#include "hams/Heisenberg.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;
using namespace itensor;
using std::vector;

int
main(int argc, char* argv[])
    {
		
	/* initialize random seed: */
    //srand (time(NULL));	

	int N = 19; // Number of sites.
	auto sites = SpinHalf(N); // Define lattice sites.	
	auto psi = MPS(sites); // Create a random MPS on the sites.
	auto psiStart = psi;

// Create the measurement operator Omega(mu+-1) = 1/2 * (expm(-1i*phi*A)+-1i*expm(1i*phi*A)) which measures on the site indexSite;
	int indexSite = 5;
	auto ampo1 = AutoMPO(sites); // Mu+1 Measurement Operator on indexSite
	auto ampo2 = AutoMPO(sites); // Mu-1 Measurement Operator on indexSite
	auto ampo3 = AutoMPO(sites); // Omega^\dagger Omega for Mu+1 Operator on indexSite
	auto ampo4 = AutoMPO(sites); // Omega^\dagger Omega for Mu-1 Operator on indexSite
	auto ampo5 = AutoMPO(sites);

	ampo1 += 0.3608,"Sz",indexSite; 			// Assign real value to entry (1,1)
	ampo1 += 0.3608/2,"Id",indexSite;			// Assign real value to entry (1,1)
	ampo1 += 0.3608*Cplx_i,"Sz",indexSite;		// Assign complex value to entry (1,1)
	ampo1 += 0.3608/2*Cplx_i,"Id",indexSite;	// Assign complex value to entry (1,1)

	ampo1 += -0.6082,"Sz",indexSite;			// Assign real value to entry (2,2)
	ampo1 += 0.6082/2,"Id",indexSite;			// Assign real value to entry (2,2)
	ampo1 += -0.6082*Cplx_i,"Sz",indexSite;		// Assign complex value to entry (2,2)
	ampo1 += 0.6082/2*Cplx_i,"Id",indexSite;	// Assign complex value to entry (2,2)

	ampo2 += 0.6082,"Sz",indexSite;				// Real (1,1)
	ampo2 += 0.6082/2,"Id",indexSite;			// Real (1,1)
	ampo2 += -0.6082*Cplx_i,"Sz",indexSite;		// Complex (1,1)
	ampo2 += -0.6082/2*Cplx_i,"Id",indexSite;	// Complex (1,1)

	ampo2 += -0.3608,"Sz",indexSite;
	ampo2 += 0.3608/2,"Id",indexSite;
	ampo2 += 0.3608*Cplx_i,"Sz",indexSite;
	ampo2 += -0.3608/2*Cplx_i,"Id",indexSite;
	
	ampo3 += 0.2603,"Sz",indexSite;
	ampo3 += 0.2603/2,"Id",indexSite;
	ampo3 += -0.7397,"Sz",indexSite;
	ampo3 += 0.7397/2,"Id",indexSite;

	ampo4 += 0.7397,"Sz",indexSite;
	ampo4 += 0.7397/2,"Id",indexSite;
	ampo4 += -0.2603,"Sz",indexSite;
	ampo4 += 0.2603/2,"Id",indexSite;

	ampo5 += 0.3608,"Sz",1;
	ampo5 += 0.3608/2,"Id",1;
	ampo5 += 0.3608*Cplx_i,"Sz",1;
	ampo5 += 0.3608/2*Cplx_i,"Id",1;

	ampo5 += -0.6082,"Sz",1;
	ampo5 += 0.6082/2,"Id",1;
	ampo5 += -0.6082*Cplx_i,"Sz",1;
	ampo5 += 0.6082/2*Cplx_i,"Id",1;
	
	auto omegaPlus = MPO(ampo1);
	auto omegaMinus = MPO(ampo2);
	auto omegaPlusDagOmegaPlus = MPO(ampo3);
	auto omegaMinusDagOmegaMinus = MPO(ampo4);
	auto startingMeasurement = MPO(ampo5);


// Create the Heisenberg XYZ-chain Hamiltonian
	auto ampo = AutoMPO(sites);
	Real jInt = 1;
	for(int b = 1; b < N; ++b)
	{
		ampo += jInt/2,"S+",b,"S-",b+1;
		ampo += jInt/2,"S-",b,"S+",b+1;
		ampo += jInt,"Sz",b,"Sz",b+1;
	}
	auto H = MPO(ampo);
	
	
// Create the Correlators S1 = Sz(x-2)Sz(x+2), S2 = Sz(x-1)Sz(x+1), S3 = Sz(x)Sz(x+1), S4 = Sz(x)Sz(x+2)
/*	auto ampoS1 = AutoMPO(sites);
	ampoS1 += 1,"Sz",indexSite-2,"Sz",indexSite+2;
	auto S1 = MPO(ampoS1);
	
	auto ampoS2 = AutoMPO(sites);
	ampoS2 += 1,"Sz",indexSite-1,"Sz",indexSite+1;
	auto S2 = MPO(ampoS2);
	
	auto ampoS3 = AutoMPO(sites);
	ampoS3 += 1,"Sz",indexSite,"Sz",indexSite+1;
	auto S3 = MPO(ampoS3);
	
	auto ampoS4 = AutoMPO(sites);
	ampoS4 += 1,"Sz",indexSite,"Sz",indexSite+2;
	auto S4 = MPO(ampoS4);
	*/
	
// Ground State Search via DRMG algorithm on psi2
  // Define DMRG sweeps
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
  // Do DMRG on psi
    println("GS Calculation: Started");
    dmrg(psi,H,sweeps,"Quiet");
    println("GS Calculation: Completed");
	
	
// Apply Measurement operator (in this case Sz) on site measureSite to create a local measurement
	auto argsMeasure = Args("Cutoff",1E-9,"Maxm",3000);
//	fitApplyMPO(psi,omegaPlus,psi,argsMeasure);
//	psi /= psi.norm();

	
// Time Evolution for psi.
 	auto tau = 0.01; // Time step size
	auto ttotal = 5; // Total time
	auto expH = toExpH<ITensor>(ampo,tau*Cplx_i); //auto expH = toExpH<ITensor>(ampo,tau*Cplx_i) if a complex exponential function is desired

  // Arguments for the fitApplyMPO algorithm
	auto args = Args("Cutoff",1E-9,"Maxm",3000);
	auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
	
	auto epsilon = ((double) rand() / (RAND_MAX));
	Real probability;
	Real probabilityTest;
	auto phi = psi;
	
	println("Time Evolution: Started");
	
// Pre-loop stuff

	// Selection of measurement
	psi.position(1);
	epsilon = ((double) rand() / (RAND_MAX));
	probability = real(psiHphi(psi,omegaPlusDagOmegaPlus,psi));
	probabilityTest = real(psiHphi(psi,omegaMinusDagOmegaMinus,psi));
	
	if(epsilon <= probability)
	{
		fitApplyMPO(psi,omegaPlus,psi,args);
		psi /= psi.norm();
		printfln("eps %.10f, prob %.10f, prob2 %.10f, mu = +1",epsilon,probability,probabilityTest);
	}
	else
	{
		fitApplyMPO(psi,omegaMinus,psi,args);
		psi /= psi.norm();
		printfln("eps %.10f, prob %.10f, prob2 %.10f, mu = -1",epsilon,probability,probabilityTest);
	}
	psi /= psi.norm();
	psi.position(1);
	//fitApplyMPO(psi,startingMeasurement,psi,args);
	//psi /= psi.norm();	
	
	for(int n = 1; n <= nt; ++n)
    {
		
		//*// Calculate Spin Projection
		std::string filename = "./StrongMeasurementsOnSite5/Projection";
		filename += to_string(n);
		filename += ".txt";
		
		ofstream dataFile; // See http://www.cplusplus.com/doc/tutorial/files/ for a tutorial on writing to files :)
		dataFile.open(filename);
		dataFile << "site Sx Sy Sz\n";
		dataFile.close();

		for(int k = 1; k <= N; ++k)
		{
			psi.position(k); // Gauges the MPS to the specified site.
			
			Real SxProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sx",k)*psi.A(k)).toComplex());
			Real SyProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sy",k)*psi.A(k)).toComplex());
			Real SzProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sz",k)*psi.A(k)).toComplex());
			
			dataFile.open(filename,ios::app); // ios:app tells the datafile that it should append.
			dataFile << k << " " << SxProjection << " " << SyProjection << " " << SzProjection << "\n";
			dataFile.close();
		}
		//*// Done Calculating Spin Projection

		//*// Calculate Correlators
		/*std::string filename2 = "./DataCorrelatorPhi500Site7/Correlator";
		filename2 += to_string(n);
		filename2 += ".txt";
		
		ofstream dataFile2; // See http://www.cplusplus.com/doc/tutorial/files/ for a tutorial on writing to files :)
		dataFile2.open(filename2);
		dataFile2 << "Sz(x-2)Sz(x+2) Sz(x-1)Sz(x+1) Sz(x)Sz(x+1) Sz(x)Sz(x+2)\n";
		dataFile2.close();
		
		Real CorrelatorXM2XP2 = real(psiHphi(psi,S1,psi));
		Real CorrelatorXM1XP1 = real(psiHphi(psi,S2,psi));
		Real CorrelatorXXP1 = real(psiHphi(psi,S3,psi));
		Real CorrelatorXXP2 = real(psiHphi(psi,S4,psi));
		
		dataFile.open(filename2,ios::app);
		dataFile << CorrelatorXM2XP2 << " " << CorrelatorXM1XP1 << " " << CorrelatorXXP1 << " " << CorrelatorXXP2 << "\n";*/
		//*// Done Calculating Correlators
		
		//*// Calculate Measurement
		if(n % 5 == 0) // conduct measurement
		{
			// Selection of measurement
			psi.position(1);
			epsilon = ((double) rand() / (RAND_MAX));
			probability = real(psiHphi(psi,omegaPlusDagOmegaPlus,psi));
			probabilityTest = real(psiHphi(psi,omegaMinusDagOmegaMinus,psi));
			
			if(epsilon <= probability)
			{
				fitApplyMPO(psi,omegaPlus,psi,args);
				psi /= psi.norm();
				printfln("eps %.10f, prob %.10f, prob2 %.10f, mu = +1",epsilon,probability,probabilityTest);
			}
			else
			{
				fitApplyMPO(psi,omegaMinus,psi,args);
				psi /= psi.norm();
				printfln("eps %.10f, prob %.10f, prob2 %.10f, mu = -1",epsilon,probability,probabilityTest);
			}
		}
		//*// Done Calculating Measurement
		
		//*// Time Evolution
		psi.position(1);
		fitApplyMPO(psi,expH,psi,args);
		psi /= psi.norm();
		printfln("Step %.10f of %.10f",n,nt);
		//*// Time Evolution Done

    }
	println("Time Evolution: Completed");

return 0;
    }
