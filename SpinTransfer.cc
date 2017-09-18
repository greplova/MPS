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
	//auto psi = MPS(sites); // Create a random MPS on the sites.
	//auto psiStart = psi;

// Initial state
	InitState init(sites,"Up");
	for(int j = 1; j <= N; j += 2)
	{
    	init.set(j,"Dn");
	}
	MPS psi(init);
	
	
	
// Create the spin transfer operators which transfers all spin one to the right.
	auto ampoST = AutoMPO(sites);
	for(int l = 1; l<N;l++)
	{
		ampoST += 1,"S+",l,"S-",l+1;
		ampoST += 1,"S-",l,"S+",l+1;
		//ampoST += 1,"S+",l,"S-",l,"S+",l+1,"S-",l+1;
		//ampoST += 1,"S-",l,"S-",l+1,"S+",l,"S+",l+1;
	}

	auto spinTOp = MPO(ampoST);
	
	PrintData(spinTOp);

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
	
	std::string filename = "./Projection";
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
	
	
	auto args = Args("Cutoff",1E-9,"Maxm",3000);
	fitApplyMPO(psi,spinTOp,psi,args);

	std::string filename2 = "./Projection2";
	filename2 += ".txt";
	
	ofstream dataFile2; // See http://www.cplusplus.com/doc/tutorial/files/ for a tutorial on writing to files :)
	dataFile2.open(filename2);
	dataFile2 << "site Sx Sy Sz\n";
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
	}
return 0;
    }
