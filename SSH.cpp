/**
 *	TIGHT-BINDING MODEL FOR SSH MODEL (TBG)
 *  (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de
 * 	
 * This code allows a selfconsistent calculation of the temperature dependent equilibrium displacement u aswell as the Peierls-undiced electronic dynamics + nuclear Ehrenfest synamcis  
 *  
 * Output:
 * -equilibrium u, rho
 * -t.-d. rho(t), u(t), p(t) --> E_electronic(t), E_nuclear(t) 
 * 
 * Necessary input:
 * -None
 *
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>


// PARAMETERS ##########################################################

// intrinsic parameters
// electronic
#define NATOM     2                            						    // Number of atoms (orbitals) in unit cell
#define tt        0.4                                                  // bare hopping
#define Alpha     0.28                                                  // electron phonon coupling
#define BETA      10000.                       					     	// inverse temperature
#define MU        0.0                                                   // Chemical potential

// lattice parameters
#define K         0.55                                                  // spring constant
#define w0        0.00347  											    // Free Oscillation Frequqency of ions: 2*np.sqrt(K/M) (in eV)
#define NN        1024                                                  // # of nuclei -> # of electrons (k-points) = 0.5*NN half filling

// numerical paramters
#define dev       1e-10                    					        	// exit treshold for while loop in groundstate()
#define DELTA     1e-7												    // correction prefactor for chemical potential in groundstate()
#define u_init    0.01										     	    // initially guessed displacement

// propagations parameters
#define starttime 0.0 
#define endtime   10000.0
#define timesteps 1e5   
#define fac       10                                                     // Every fac-th value is stored on disc		

#define SWITCH_OFF 1.0                                                  // switch nuclear dynamics on (1.0) / off (0.0)

// Electronic driving
#define w_peierls      0.188                                            // Frequency of Applied Field (in eV)
#define AA_peierls     0.10                                              // Amplitude of Applied Field (E**2 in eV**2*A)
#define tm_peierls     1500.                                             // mean expectation of Gauss pulse defiel1()
#define sig_peierls    274.                                              // standard deviation of Gauss pulse defiel1() -> 300 fs FWHM on intensity orfile



// CALCULATION OPTIONS #################################################

#ifndef NO_MPI
    #include <mpi.h>
#endif

#define PI        3.14159265359

using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);



// DEFINITION OF FUNCTIONS #############################################

//LAPACK (Fortran 90) functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//routine to find eigensystem of matrix
extern "C" {
/** 
 *  Computes the eigenvalues and, optionally, the eigenvectors for a Hermitian matrices H
 */
    void zheev_(char* jobz, char* uplo, int* N, cdouble* H, int* LDA, double* W, cdouble* work, int* lwork, double* rwork, int *info);
}
//'N','V':  Compute eigenvalues only, and eigenvectors
char    jobz = 'V';       
//'U','L':  Upper, Lower triangle of H is stored 
char    uplo = 'U';  
// The order of the matrix H.  NATOM >= 0
int     matsize = NATOM;    
// The leading dimension of the array H.  lda >= max(1, NATOM)
int     lda = NATOM;             
// The length of the array work.  lwork  >= max(1,2* NATOM-1)
int     lwork = 2*NATOM-1;    
// dimension (max(1, 3* NATOM-2))
double  rwork[3*NATOM-2];  
// dimension (MAX(1,LWORK))
cdouble work[2*NATOM-1];  
// Info
int	    info;


void diagonalize(cvec &Hk, dvec &evals)
{
/**
 *  Diagonalization of matrix Hk. Stores eigenvalues in real vector evals and eigenvectors in complex vector Hk
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */
    zheev_(&jobz, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}


//INLINE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline int fq(int i, int j, int N)
/**
 *  MAT[i,j] = Vec[fq(i,j,N)] with row index i and column index j
 */
{
    return i*N+j;
}


inline double delta(int a, int b)
/**
 *  Delta function
 */
{
	if (a==b)
		return 1.;
	else
		return 0.;
}


template <class Vec>
inline void print(Vec vec)
/**
 * Print out vector
 */
{
	for(int i=0; i<vec.size(); i++)
		{
	    	cout << vec[i] << " ";
	    }	
	cout << endl;
}


template <class Vec>
inline double get_deviation(Vec &M1, Vec &M2)
/**
 * Calculates deviation of two different vectors M1 und M2
 */
{
    double deviations = 0.;
	for(int i=0; i<M1.size(); i++)
	{
		deviations += abs(M1[i]-M2[i]);  	
    }
    return deviations;
}


inline double fermi(double energy, double mu)
{
/**
 *	Fermi distribution:
 *	-energy: Energy eigenvalue
 *	-mu: Chemical potential
 */
    return 1./(exp((energy-mu)*BETA) + 1.);
}


inline double gauss(double time, double delay, double sigma)
/**
 *	Normalized Gauss distribution
 *	-time: time coordinate
 *	-delay: mean expectation value
 *	-sigma: standard deviation 
 **/
{
	return 1./(sigma*sqrt(2.*PI))*exp(-0.5*pow((time-delay)/sigma,2.));
}


inline double A_peierls(double time)
{
/**
 *	Peierls field for electrons in x-direction:
 *  -time: Real time coordinate
 */
    return AA_peierls*exp(-0.5*pow((time-tm_peierls)/sig_peierls,2.))*sin(w_peierls*time);
}


// VOID FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ReadIn(vector<dvec> &MAT, const string& filename)
{
/**
 *	Read in real valued matrix
 */
	ifstream in(filename);
	string record;
	if(in.fail()){
		cout << "file" << filename << "could not be found!" << endl;
	}
	while (getline(in, record))
	{
		istringstream is( record );
		dvec row((istream_iterator<double>(is)),
		istream_iterator<double>());
		MAT.push_back(row);
	}
	in.close();
}

template <class Vec>
void times(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product of quadratic matrices: $C = A \cdot B$
 */
{
    int dim = sqrt(A.size());
	Vec TEMP(dim*dim);
    // Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
	    for(int j=0; j<dim; j++) {
		    TEMP[fq(j,i,dim)] = B[fq(i,j,dim)];
		   }
    }
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += A[fq(i,k,dim)]*TEMP[fq(j,k,dim)]; 
			}
		}
	}	
}


template <class Vec>
void times_dn(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of first factor: $C = A^\dagger \cdot B$
 */
{
	int dim = sqrt(A.size());
	Vec TEMP1(dim*dim);
	Vec TEMP2(dim*dim);
	// Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++) {
			TEMP1[fq(j,i,dim)] = A[fq(i,j,dim)];
			TEMP2[fq(j,i,dim)] = B[fq(i,j,dim)];
		}
	}		
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += conj(TEMP1[fq(i,k,dim)])*TEMP2[fq(j,k,dim)];
			}
		}
	}		
}


template <class Vec>
void times_nd(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of second factor: $C = A \cdot B^\dagger$
 */
{
	int dim = sqrt(A.size());	
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
					C[fq(i,j,dim)] += A[fq(i,k,dim)]*conj(B[fq(j,k,dim)]);
			}
		}
	}	
}


void set_Hk0(double k, double &u, cvec &Hk)
/**
 *	Set electronic eq. Hamiltonian
 *  -k: k-pint of 1d reduced Brilluoin zone
 *  -u: displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 */
{
	Hk[0] = -2.*tt*cos(k);
	Hk[1] = 4.*Alpha*u*sin(k);
	Hk[2] = 4.*Alpha*u*sin(k);
	Hk[3] = +2.*tt*cos(k); 				
}


void set_Hk(double k, double &u, cvec &Hk, double time)
/**
 *	Set t.d. electronic Hamiltonian
 *  -k: k-pint of 1d reduced Brilluoin zone
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -time: Real time coordinate
 */
{
	Hk[0] = -2.*tt*cos(k-A_peierls(time));
	Hk[1] = 4.*Alpha*u*sin(k-A_peierls(time));
	Hk[2] = 4.*Alpha*u*sin(k-A_peierls(time));
	Hk[3] = +2.*tt*cos(k-A_peierls(time)); 				
}


void set_dHkdu0(double k, double &u, cvec &Hk)
/**
 *	Set derivative of e.q. electronic Hamiltonian by displacement u
 *  -k: k-pint of 1d reduced Brilluoin zone
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 */
{
	Hk[0] = 0.;
	Hk[1] = 4.*Alpha*sin(k);
	Hk[2] = 4.*Alpha*sin(k);
	Hk[3] = 0.; 				
}


void set_dHkdu(double k, double &u, cvec &Hk, double time)
/**
 *	Set derivative of electronic Hamiltonian by displacement u
 *  -k: k-pint of 1d reduced Brilluoin zone
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -time: Real time coordinate
 */
{
	Hk[0] = 0.;
	Hk[1] = 4.*Alpha*sin(k-A_peierls(time));
	Hk[2] = 4.*Alpha*sin(k-A_peierls(time));
	Hk[3] = 0.; 				
}


void SET_u0(double &u, cvec &Hk, dvec &BZ, vector<cvec> &RHO_0, int &numprocs, int &myrank)
/**
 *	Self-consistent calculation of u by density matrix
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{	
	cvec TEMP(4);
	u = 0.;
	for(int k=myrank; k<NN/2; k+=numprocs)
	{		
		set_dHkdu0(BZ[k], u, Hk);
		times(RHO_0[k], Hk, TEMP);
		u += -2./double(4.*K*NN)*real(TEMP[fq(0,0,2)]+TEMP[fq(1,1,2)]);      
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &u, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
}


void set_Rhok(double &u, dvec &evals, cvec &Hk, dvec &BZ, vector<cvec> &RHO_0, int &numprocs, int &myrank)
{	
/**
 *	Set density matrix in k-orbital basis
 *  -evals: Real vector[2] of eigenvalues
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
	cvec TEMP(4,0.);
	cvec RHO_BAND(4,0.);
	
	for(int k=myrank; k<NN/2; k+=numprocs)
	{
		set_Hk0(BZ[k], u, Hk);
		diagonalize(Hk, evals);                                             
		
		for(int i=0; i<2; i++)
			RHO_BAND[fq(i,i,2)] = fermi(evals[i], MU);         	

		times(RHO_BAND, Hk, TEMP);                                             
		times_dn(Hk, TEMP, RHO_0[k]);
	}
}


void groundstate(double &u, double &p, dvec &evals, cvec &Hk, dvec &BZ, vector<cvec> &RHO_0, vector<dvec> &E_TOT, int &numprocs, int &myrank)
/**
 *	Calculate initial displacement u0
 *  -u: nuclear displacement 
 *  -p: nuclear momentum
 *  -evals: Real vector[2] of eigenvalues
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{	
	int count = 0;                                                      // count # of loops of self-consistency
	u = u_init;
	
	double N_tot;
	double u_old;	    
	double deviation = 1.0;
	
	double E_tot; 
	double E_nuc;                                               
	
	cvec TEMP(4);
		
	while(deviation > dev)
	{
		count++;				
		u_old = u;	    
		N_tot = 0.;
		E_tot = 0.;

		set_Rhok(u, evals, Hk, BZ, RHO_0, numprocs, myrank);
		SET_u0(u, Hk, BZ, RHO_0, numprocs, myrank);

		deviation = abs(u-u_old);
	
		for(int k=myrank; k<NN/2; k+=numprocs)
		{
			// Calculation of total electronic energy per unit cell (1 electron per unit cell)
			set_Hk0(BZ[k], u, Hk);
			times(RHO_0[k], Hk, TEMP);
			E_tot += real(TEMP[fq(0,0,2)]+TEMP[fq(1,1,2)])/(double(NN)/2.);	
			N_tot += real(RHO_0[k][fq(0,0,2)]+RHO_0[k][fq(1,1,2)])/(double(NN)/2.);	
		}	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &E_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &N_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
		// Energy of nuclear system per unit cell (2 atoms per unit cell)
		E_nuc = 2.*2.*K*u*u;
		
		if(myrank==0){
			cout << "loop #" << count << ": deviation = " << deviation << endl;
			cout << "u: " << u << endl;
			cout << "gap: " << 8.*Alpha*u  << endl;
			cout << "total particle number N = " << N_tot << endl;
			cout << "E_el_tot: " << E_tot << endl;
		    cout << "E_nu_tot: " << E_nuc << endl;
			cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
		}	 
	}
	E_TOT[0][0] = E_tot;
	E_TOT[0][1] = E_nuc;
	p=0.0;
	if(myrank==0)
	{
		ofstream myfile ("PARAMETERS.txt");
		if (myfile.is_open())
		{
			myfile <<  u  << " " <<  8.*Alpha*u << endl;
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{		
		ofstream myfile ("RHO0.txt");
		if (myfile.is_open())
		{
			for(int k=0; k<NN/2; k++)
			{
				for(int i=0; i<4; i++)
				{
					myfile << RHO_0[k][i] << " " ;
				}
			myfile  << endl;
			}
		myfile.close();
		}	
		else cout << "Unable to open file" << endl;		
	}
}


void Hk_bands(vector<dvec> BANDS, double &u, cvec &Hk, dvec &evals, dvec &BZ, double time, const string& filename)
/**
 *	Calculate equlibrium bandstructure
 *  -evals: Real vector[2] of eigenvalues
 *  -u: Displacement 
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{
	for(int k=0; k<NN/2; k++)
	{
		set_Hk0(BZ[k], u, Hk);
		diagonalize(Hk, evals);
		for(int m=0; m<2; m++)
			BANDS[k][m] = evals[m];
	}
	ofstream myfile (filename);
	if (myfile.is_open())
	{
		for(int k=0; k<NN/2; k++)
		{
			for(int m=0; m<2; m++)
			{
				myfile << BANDS[k][m] << " " ;
			}
		myfile  << endl;
		}
	myfile.close();
	}
    else cout << "Unable to open file" << endl;
}

				
void set_dRHOdt(cvec &TEMP1, cvec &TEMP2, cvec &RHO_t_tk, cvec &dRHO_dt, cvec &Hk)
/**
 *  Calculation of the time-derivative of the density matrix
 *  -TEMP1, TEMP2: Complex helper matrix 
 *  -RHO_t_tk: Complex vector[4] of k- and time-dependent density matrix
 *  -dRHO_dt: Complex vector[4] of temporal change of density matrix
 *  -Hk: Complex vector[4] to store Hamiltonian matrix
 */
{	
	// COHERENT PART
	times(Hk, RHO_t_tk, TEMP1);										
	times(RHO_t_tk, Hk, TEMP2);
	for(int i=0; i<4; i++)
	{
		dRHO_dt[i] = -II*(TEMP1[i]-TEMP2[i]);	
	}
}	


void set_dudp(double &u, double &p, double &du, double &dp, cvec &TEMP, cvec &RHO_t_tk, cvec &dHkdu)
/**
 *  Calculate t.d. change of nuclear displacement u and momentum p
 *  -TEMP1, TEMP2: Complex helper matrix 
 *  -RHO_t_tk: Complex vector[4] of k- and time-dependent density matrix
 *  -dRHO_dt: Complex vector[4] of temporal change of density matrix
 *  -Hk: Complex vector[4] to store Hamiltonian matrix
 */
{	
	
	du = p/(4.*K/pow(w0,2.));
	times(RHO_t_tk, dHkdu, TEMP);
	dp += -2./double(NN)*real(TEMP[fq(0,0,2)]+TEMP[fq(1,1,2)]) - 4.*K*u/(double(NN)/2.);
}	


void AB2_propatation(vector<dvec> &E_TOT, dvec &evals, vector<cvec> &RHO_0, vector<cvec> &dRHO_dt0, vector<cvec> &dRHO_dt1, double &u, double &p, cvec &Hk, dvec &BZ, vector<cvec*> RHO_t, vector<dvec> &ORDER_t, int &numprocs, int &myrank)
/**
 *	Two-step Adams-Bashforth prdictor corrector method for propagation of full (el+nuc) system
 *  -E_TOT: Vector[timesteps] of real vectors[2] to store enrgies
 *  -evals: Real vector[2] of eigenvalues
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -dRHO_dt0, dRHO_dt1: Vector of complex vectors[64] to store change of density matrix
 *  -u: nuclear displacement 
 *  -p: nuclear momentum
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_t_tk: Vector[4] ofcomplex vector pointers to store t.-d. density matrix
 *  -ORDER_t: Vector[timesteps] of real vectors[2] to store t.-d. u(t) and p(t)
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{
	double du0, du1, dp0, dp1, uu, E_tot, E_nuc, time, h;
	cvec *temp0, *temp1, *temp2;
    h = (endtime-starttime)/timesteps;
    
    cvec TEMP1(4,0.);												    // Helper arrays for set_dRhodt()
	cvec TEMP2(4,0.);
    
    //p=10.;
    
	// initial Order
	ORDER_t[0][0] = u;                                                  // intial nuclear displacement
	ORDER_t[0][1] = p;                                                  // intial nuclear momentum
	
	for(int k=myrank; k<NN/2; k+=numprocs)
	{
		for(int i=0; i<4; i++)
		{
			(*RHO_t[fq(0, k, NN/2)])[i] = RHO_0[k][i];
		}	
	}
	
	// Propagation	
	for(int t=0; t<timesteps-1; t++)
	{   
		// 1st Euler step 
		if(t==0)
		{    
			dp0 = 0.;  
			for(int k=myrank; k<NN/2; k+=numprocs)
			{ 
				set_dHkdu(BZ[k], ORDER_t[0][0], Hk, h*double(t));
				set_dudp(ORDER_t[0][0], ORDER_t[0][1], du0, dp0, TEMP1, RHO_t[fq(0, k, NN/2)][0], Hk); 
				set_Hk(BZ[k], ORDER_t[0][0], Hk, h*double(t));
				set_dRHOdt(TEMP1, TEMP2, RHO_t[fq(0,k,NN/2)][0], dRHO_dt0[k], Hk);
				for(int i=0; i<4; i++)
				{
					(*RHO_t[fq(1,k,NN/2)])[i] = (*RHO_t[fq(0,k,NN/2)])[i] + h*dRHO_dt0[k][i]; 
					RHO_0[k][i] = (*RHO_t[fq(1,k,NN/2)])[i]; 
				}
			}		
			
#ifndef 	NO_MPI		
			MPI_Allreduce(MPI_IN_PLACE, &dp0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
			//cout << "dp0: " << dp0 <<  endl;
			ORDER_t[t+1][0] = ORDER_t[t][0] + h*du0*SWITCH_OFF;
			ORDER_t[t+1][1] = ORDER_t[t][1] + h*dp0*SWITCH_OFF; 
			
			E_tot = 0.0;
			for(int k=myrank; k<NN/2; k+=numprocs)
			{
				// calculation of total electronic energy per unit cell
				set_Hk(BZ[k], ORDER_t[t+1][0], Hk, h*double(t+1));
				times(RHO_0[k], Hk, TEMP1);
				for(int i=0; i<2; i++)
				{
					E_tot += real(TEMP1[fq(i,i,2)])/(double(NN)/2.);	
				}
			}
			MPI_Allreduce(MPI_IN_PLACE, &E_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			E_nuc = 2.*(2.*K*pow(ORDER_t[t][0],2.)+ pow(ORDER_t[t][2],2)/(2.*4.*K/pow(w0,2.)));
			E_TOT[t+1][0] = E_tot;
			E_TOT[t+1][1] = E_nuc; 
		}
		// Two step Adams–Bashforth method
		else
		{	// 2-step Adams predictor	
			dp0 = 0.;
			dp1 = 0.;
			for(int k=myrank; k<NN/2; k+=numprocs)
			{
				set_dHkdu(BZ[k], ORDER_t[0][0], Hk, h*double(t-1));
				set_dudp(ORDER_t[t-1][0], ORDER_t[t-1][1], du0, dp0, TEMP1, RHO_t[fq(0, k, NN/2)][0], Hk);
				set_dHkdu(BZ[k], ORDER_t[0][0], Hk, h*double(t));
				set_dudp(ORDER_t[t][0], ORDER_t[t][1], du1, dp1, TEMP1, RHO_t[fq(1, k, NN/2)][0], Hk);
				set_Hk(BZ[k], ORDER_t[t-1][0], Hk, h*double(t-1));					
				set_dRHOdt(TEMP1, TEMP2, RHO_t[fq(0,k,NN/2)][0], dRHO_dt0[k], Hk);
				set_Hk(BZ[k], ORDER_t[t][0], Hk, h*double(t));					
				set_dRHOdt(TEMP1, TEMP2, RHO_t[fq(1,k,NN/2)][0], dRHO_dt1[k], Hk);
				
				for(int i=0; i<4; i++)
				{
					// P_{n+1} = y_{n} + 3/2*h*f(t_{n},y_{n}) - 0.5*h*f(t_{n-1},y_{n-1})	
					(*RHO_t[fq(2,k,NN/2)])[i] = (*RHO_t[fq(1,k,NN/2)])[i] + h*(3./2.*dRHO_dt1[k][i] - 0.5*dRHO_dt0[k][i]); 		
					RHO_0[k][i] = (*RHO_t[fq(2,k,NN/2)])[i]; 
				}
			}			
							
#ifndef 	NO_MPI		
			MPI_Allreduce(MPI_IN_PLACE, &dp0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, &dp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
			ORDER_t[t+1][0] = ORDER_t[t][0] + h*(3./2.*du1 - 0.5*du0)*SWITCH_OFF;
			ORDER_t[t+1][1] = ORDER_t[t][1] + h*(3./2.*dp1 - 0.5*dp0)*SWITCH_OFF;  
				
			// 2-step Moulton corrector
			dp0 = 0.;
			for(int k=myrank; k<NN/2; k+=numprocs)
			{
				set_dHkdu(BZ[k], ORDER_t[0][0], Hk, h*double(t));
				set_dudp(ORDER_t[t+1][0], ORDER_t[t+1][1], du0, dp0, TEMP1, RHO_t[fq(2, k, NN/2)][0], Hk);
				set_Hk(BZ[k], ORDER_t[t+1][0], Hk, h*double(t+1));				
				set_dRHOdt(TEMP1, TEMP2, RHO_t[fq(2,k,NN/2)][0], dRHO_dt0[k], Hk);   // comment out -> works!!!
				
				for(int i=0; i<4; i++)
				{
					// y_{n+1} = y_{n} + 1/2*h*(f(t_{n+1},P_{n+1}) + f(t_{n},y_{n}))
					(*RHO_t[fq(2,k,NN/2)])[i] = (*RHO_t[fq(1,k,NN/2)])[i] + 0.5*h*(dRHO_dt0[k][i] + dRHO_dt1[k][i]); 		   
					RHO_0[k][i] = (*RHO_t[fq(2,k,NN/2)])[i]; 
				}
			
				// Cyclic exchange of pointers
				temp0 = RHO_t[fq(0,k,NN/2)];
				temp1 = RHO_t[fq(1,k,NN/2)];
				temp2 = RHO_t[fq(2,k,NN/2)];
				
				RHO_t[fq(0,k,NN/2)] = temp1;
				RHO_t[fq(1,k,NN/2)] = temp2;
				RHO_t[fq(2,k,NN/2)] = temp0;
			}				

#ifndef 	NO_MPI		
			MPI_Allreduce(MPI_IN_PLACE, &dp0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
			ORDER_t[t+1][0] = ORDER_t[t][0] + 0.5*h*(du0 + du1)*SWITCH_OFF;
			ORDER_t[t+1][1] = ORDER_t[t][1] + 0.5*h*(dp0 + dp1)*SWITCH_OFF;
			
			E_tot = 0.0;
			for(int k=myrank; k<NN/2; k+=numprocs)
			{
				// calculation of total electronic energy per unit cell
				set_Hk(BZ[k], ORDER_t[t+1][0], Hk, h*double(t+1));
				times(RHO_0[k], Hk, TEMP1);
				for(int i=0; i<2; i++)
				{
					E_tot += real(TEMP1[fq(i,i,2)])/(double(NN)/2.);	
				}
			}
			MPI_Allreduce(MPI_IN_PLACE, &E_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			E_nuc = 2.*(2.*K*pow(ORDER_t[t][0],2.)+ pow(ORDER_t[t][2],2)/(2.*4.*K/pow(w0,2.)));
			E_TOT[t+1][0] = E_tot;
			E_TOT[t+1][1] = E_nuc;
			
			if(myrank==0) 
			{	
				cout << " time step: " << t+1 <<  endl;
				cout << " Peierls Driving: " << A_peierls(h*double(t+1)) << endl;
				cout << " u: "  << ORDER_t[t+1][0] << " p: "  << ORDER_t[t+1][1] << endl;
				cout << " E_tot electronic: " << E_TOT[t+1][0] << " E_tot nuclear: " << E_TOT[t+1][1] << endl;
				cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
			}		
		}
	}
	if(myrank==0)
	{
		ofstream myfile ("ORDER_t.txt");
		if (myfile.is_open())
		{
			for(int t=0; t<timesteps; t+=fac)
			{
				for(int i=0; i<2; i++)
					myfile << ORDER_t[t][i]  << " ";
				myfile  << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile ("ENERGY_t.txt");
		if (myfile.is_open())
		{
			for(int t=0; t<timesteps-1; t+=fac)
			{
				for(int i=0; i<2; i++)
					myfile << E_TOT[t][i]  << " ";
				myfile  << endl;	
			}
		}
		else cout << "Unable to open file" << endl;	
	}	
}


void ExpU_propatation(double &u, double &p, dvec &evals, cvec &Hk, dvec &BZ, vector<cvec> &RHO_0, int &numprocs, int &myrank)
/**
 *	Unitary propagation of electronic density operator for p,u = const (use SWITCH_OFF 0.0 for comparison with AB2_propatation()!)
 *  -u: nuclear displacement 
 *  -p: nuclear momentum
 *  -evals: Real vector[2] of eigenvalues
 *  -Hk: Complex vector[4] to store Hamiltonian
 *  -BZ: Real vector[NN/2] of 1d reduced Brilluoin zone
 *  -RHO_0: Vector[num_k] of real vectors[4] to store density matrix 
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{
	const int dimH = NATOM*NATOM;
	const double h = (endtime-starttime)/timesteps;	
	
	cvec *TEMP0 = new cvec(dimH);									
	cvec *TEMP1 = new cvec(dimH);									
	cvec *TEMP2 = new cvec(dimH);
	cvec *TEMP3 = new cvec(dimH);
	dvec E_t(timesteps, 0.);

	int count = 0;
	
	// Propagation
	for(int k=myrank; k<NN/2; k+=numprocs)	
	{   
		if(myrank==0) 
		{	
			cout << " k: " << k <<  endl;
		}	
		for(int i=0; i<NATOM*NATOM; i++)
			(*TEMP0)[i] = RHO_0[k][i];
		
		for(int t=0; t<timesteps-1; t++)
		{ 
			set_Hk(BZ[k], u, Hk, h*(double(t)+0.5));;  // Euler Mid-point rule
			diagonalize(Hk, evals);                                      
			
			for(int i=0; i<NATOM; i++)
			{
				for(int j=0; j<NATOM; j++)
				{
					(*TEMP1)[fq(i,j,NATOM)] = exp(-II*evals[i]*h)*delta(i,j);
				}
			}		
			
			times(TEMP1[0], Hk, TEMP2[0]);                              
			times_dn(Hk, TEMP2[0], TEMP1[0]);
			
			times(TEMP0[0], TEMP1[0], TEMP2[0]);                        // propagation: U(t+h)rho(t)U(t+h)^+
			times_dn(TEMP1[0], TEMP2[0], TEMP0[0]);
				
		// calculation of total energy	
			set_Hk(BZ[k], u, Hk, h*double(t+1));
			times(TEMP0[0], Hk, TEMP1[0]);
			for(int i=0; i<NATOM; i++)
			{
				E_t[t+1] += real((*TEMP1)[fq(i,i,NATOM)])/(double(NN)/2.);
			}								
		}
		count++;
	}
#ifndef NO_MPI		
	MPI_Allreduce(MPI_IN_PLACE, &E_t[0], timesteps, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
	// write td mean field parameters to file
	if(myrank==0)
	{
		ofstream myfile ("ENERGY_t_EXP.dat");
		if (myfile.is_open())
		{
			for(int t=0; t<timesteps; t+=fac)
			{
				myfile << E_t[t] << endl;
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
delete TEMP0, TEMP1, TEMP2, TEMP3;
}

// main() function #####################################################

int main(int argc, char * argv[])
{
    //************** MPI INIT ***************************
  	int numprocs=1, myrank=0, namelen;
    
#ifndef NO_MPI
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Get_processor_name(processor_name, &namelen);
    
	cout << "Process " << myrank << " on " << processor_name << " out of " << numprocs << " says hello." << endl;
	MPI_Barrier(MPI_COMM_WORLD);
    
#endif
	if(myrank==0) cout << "\n\tProgram running on " << numprocs << " processors." << endl;

   
	// DECLARATION AND INTITALIZATIO
	
	//vector of high-symmetry path vectors
	dvec BZ(NN/2);
	for(int ii=0; ii<NN/2; ii++)
	{
		BZ[ii] = PI*double(ii)/double(NN/2.)-PI/(2.);                        
	}
	
	// vector for eigenvalues
	dvec evals(2);

	// vector for Hamiltonian Hk
	cvec Hk(4);
	
	// bands
	vector<dvec> BANDS(NN/2,dvec(2));
 	
	// allocation of matrices RHO[k]
	vector<cvec> RHO_0(NN/2, cvec(4,0.0));                                       		
	
    // initial EOP
    double u;
    double p;
	
	// calculate driving field
	double h = (endtime-starttime)/timesteps;
	double Ap;
	if(myrank==0)
	{
		ofstream myfile ("DRIVING_t.txt");
		if (myfile.is_open())
		{
			for(int t=0; t<timesteps;  t+=fac)
			{			
				myfile << A_peierls(h*double(t)) << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}

	// allocation of matrix dRHO_dt[k,t]
	vector<cvec> dRHO_dt0(NN/2, cvec(4));  
	vector<cvec> dRHO_dt1(NN/2, cvec(4));    
	
	// dynamic allocation of matrix RHO[k,t]
	vector<cvec*> RHO_t(NN/2*3);                                        //<-- 2 step A-B need 
	for(int kt=0; kt<NN/2*3; kt++)
		RHO_t[kt] = new cvec(4);	
			
	// NUclear coordinates
	vector<dvec> ORDER_t(timesteps, dvec(2));                                       

	// Total electronic[0] and nuclear energy[1]                          
	vector<dvec> E_TOT(timesteps, dvec(2, 0.0));	

	
	// CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	const clock_t begin_time = clock();
	
	// caclulation of groundstate
	groundstate(u, p, evals, Hk, BZ, RHO_0, E_TOT, numprocs, myrank);
	
	if(myrank==0) 
	{
		cout << u;
		//caclulation of groundstate bands
		Hk_bands(BANDS, u, Hk, evals, BZ, 0.0, "bands0.txt");
	}
	//ExpU_propatation(u, p, evals, Hk, BZ, RHO_0, numprocs, myrank);
	AB2_propatation(E_TOT, evals, RHO_0, dRHO_dt0, dRHO_dt1, u, p, Hk, BZ, RHO_t, ORDER_t, numprocs, myrank);

	
	if(myrank==0) cout << "Calculations lasted: " << float(clock() - begin_time)/CLOCKS_PER_SEC << " seconds" << endl;
	
#ifndef NO_MPI
	MPI_Finalize();
#endif	

	// free memory	
	for(int kt=0; kt<NN/2*3; kt++)
	{                            
		delete RHO_t[kt];
	}		
}
