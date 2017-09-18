#include "core.h"
#include "sites/spinhalf.h"
#include "autompo.h"
#include "hams/Heisenberg.h"

#include <iostream>
#include <fstream>

using namespace std;
using namespace itensor;
using std::vector;

int
main(int argc, char* argv[])
    {

	int N = 20; // Number of sites.
	auto sites = SpinHalf(N); // Define lattice sites.	
	auto psi = MPS(sites); // Create a random MPS on the sites.
	auto psiStart = psi;

// Create the Heisenberg XYZ-chain Hamiltonians
	auto ampo = AutoMPO(sites);
	Real jInt = 1;
	for(int b = 1; b < N; ++b)
	{
		ampo += jInt/2,"S+",b,"S-",b+1;
		ampo += jInt/2,"S-",b,"S+",b+1;
		ampo += jInt,"Sz",b,"Sz",b+1;
	}	
	auto H = MPO(ampo);
	auto ampoTemp = ampo;

// Ground State Search via DRMG algorithm on psi2
  // Define DMRG sweeps
    Sweeps sweeps(5);
    sweeps.maxm() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
  // Do DMRG on psi
    println("GS Calculation: Started");
    dmrg(psi,H,sweeps,"Quiet");
    println("GS Calculation: Completed");
	
// Time Evolution for psi.
 	auto tau = 0.01; // Time step size
	auto ttotal = 5; // Total time
	auto expH = toExpH<ITensor>(ampo,tau*Cplx_i); //auto expH = toExpH<ITensor>(ampo,tau*Cplx_i) if a complex exponential function is desired

  // Arguments for the fitApplyMPO algorithm
	auto args = Args("Cutoff",1E-9,"Maxm",3000);
	auto nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
	Real bFieldZero = 20*jInt;
	Real pi = 3.14159265;
	Real period = 25/jInt;
	Real omega = 2*pi/period;
	auto expHTemp = expH;
	
	println("Time Evolution: Started");
	for(int n = 1; n <= nt; ++n)
    {
		std::string filename = "./DataPeriod25/Projection";
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
		
		/*std::string filename2 = "./DataCorrelator/Correlator";
		filename2 += to_string(n);
		filename2 += ".txt";
		
		ofstream dataFile2; // See http://www.cplusplus.com/doc/tutorial/files/ for a tutorial on writing to files :)
		dataFile2.open(filename2);
		for(int k = 2; k<= N; ++k)
		{
			std::string correlatorName = "S1S" + to_string(k);
			dataFile << correlatorName << " ";
		}
		dataFile2 << "\n";
		dataFile2.close();
		
		for(int k = 1; k <= N; ++k)
		{
			psi.position(k); // Gauges the MPS to the specified site.
			
			Real SxProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sx",k)*psi.A(k)).toComplex());
			Real SyProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sy",k)*psi.A(k)).toComplex());
			Real SzProjection = real((dag(prime(psi.A(k),Site))*sites.op("Sz",k)*psi.A(k)).toComplex());
			
			dataFile2.open(filename2,ios::app); // ios:app tells the datafile that it should append.
			dataFile2 << k << " " << SxProjection << " " << SyProjection << " " << SzProjection << "\n";
			dataFile2.close();
		}*/
		
		printfln("Step %.10f of %.10f",n,nt);
		
		// The following code creates the time evolution hamilatonian operator
		Real bFieldTemp = bFieldZero * cos(omega * n);
		PrintData(bFieldTemp);
		auto ampoTemp2 = ampoTemp;
		ampoTemp2 += bFieldTemp,"Sz",1;
		auto expHTemp =toExpH<ITensor>(ampoTemp2,tau*Cplx_i);
		
		psi /= psi.norm();
		psi.position(1);
		fitApplyMPO(psi,expHTemp,psi,args);

    }
	println("Time Evolution: Completed");



return 0;
    }
