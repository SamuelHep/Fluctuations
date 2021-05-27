
#include <iostream>
#include "FlucContainer.h" 

using std::cout;
using std::cerr;
using std::endl;
using std::vector;


//__________________________________________________________________________________
FlucContainer::FlucContainer(){
  
  maxOrder = 6;
  ncm.resize(2*maxOrder+1);
  C.resize(2*maxOrder+1);
  
}

//__________________________________________________________________________________
void FlucContainer::SetMaxOrder( const int order ){
	maxOrder = order;
}

//__________________________________________________________________________________
void FlucContainer::SetMoments( const vector<double> &vtmpNCM ){
	for(int i=1; i<=maxOrder; i++){
		ncm[i] = vtmpNCM[i];
	}
}

//__________________________________________________________________________________
void FlucContainer::SetMoments( const int n, const double val ){
	ncm[n] = val;
}

//__________________________________________________________________________________
void FlucContainer::SetCumulants( const int n, const double val ){
	C[n] = val;
}

//_________________________________________________________________________________________
double FlucContainer::GetCumulant( int n ){
	if(n>maxOrder){
		cerr << n <<"-th cumulant has not been calculated."<< endl;
		return -9999;
	}
	if(n==1) return ncm[1]; 
	C[n] = ncm[n];
	for(int m=1; m<=n-1; m++){
		C[n] -= TMath::Binomial(n-1,m-1)*GetCumulant(m)*ncm[n-m];
	}
	return C[n];
}

//_________________________________________________________________________________________
double FlucContainer::GetMoment( const int n ) const {
	return ncm[n];
}

//__________________________________________________________________________________
double FlucContainer::GetMomentsFromCumulants( int n ){
	if(n==0) return 1;
	if(n==1) return C[1];
	double val = 0;
	for(int i=0; i<=n-1; i++){
		val += TMath::Binomial(n-1,i)*C[i+1]*GetMomentsFromCumulants(n-1-i);
	}
	return val;
}

double FlucContainer::GetMomentsFromCumulants_v2( int n )
{
  ncm[n] = GetMomentsFromCumulants(n);
  return ncm[n];

}
