#ifndef __Fluc_h__
#define __Fluc_h__

#include <vector>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>

//__________________________________________________________________________________
class FlucContainer {

	private :
		double nevent;
		std::vector<double> ncm;
		std::vector<double> C;

	protected :
		int maxOrder;

	public :
		FlucContainer();
		virtual ~FlucContainer(){};
		void SetMaxOrder( const int order);
		void Init();
		
		// Setter
		void SetMoments( const std::vector<double> &vtmpNCM );
		void SetMoments( const int n, const double val );
		void SetCumulants( const int n, const double val );
		// Getter
		double GetCumulant( int n );
		double GetMoment( const int n ) const;
		double GetMomentsFromCumulants( int n );
		double GetMomentsFromCumulants_v2( int n );
				

};
#endif
