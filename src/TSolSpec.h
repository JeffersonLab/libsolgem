#ifndef __TSOLSPEC_H
#define __TSOLSPEC_H

#include "TList.h"

#include "THaSpectrometer.h"

#include "TSolGEMChamber.h"
#include "types.h"

class TSolSpec : public THaSpectrometer {
    public:
	TSolSpec( const char *name, const char *desc );
	virtual ~TSolSpec() {;}

	Int_t AddGEM (TSolGEMChamber& pdet);

	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();

	Int_t TrackCalc() { return 0; }

	Int_t FindVertices(TClonesArray &);
//	void MakePrefix(){ return; }

	UInt_t GetNChambers() const { return fChambers.GetEntries(); }
	TSolGEMChamber &GetChamber(Int_t i) const 
	  { return *(static_cast<TSolGEMChamber*>(fChambers.At(i))); }

	void Print();

    private:
	TList  fChambers;

    public:
	ClassDef(TSolSpec,1)

};

#endif//__TSOLSPEC_H
