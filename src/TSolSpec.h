#ifndef __TSOLSPEC_H
#define __TSOLSPEC_H

#include "TObject.h"
#include "THaSpectrometer.h"
#include "types.h"

class TSolGEMPlane;

class TSolSpec : public THaSpectrometer {
    public:
	TSolSpec( const char *name, const char *desc );
	virtual ~TSolSpec() {;}

	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();

	Int_t TrackCalc() { return 0; }

	Int_t FindVertices(TClonesArray &);

	void MakePrefix(){ return; }
    private:

    public:
	ClassDef(TSolSpec,1)

};

#endif//__TSOLSPEC_H
