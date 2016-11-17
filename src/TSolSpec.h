#ifndef __TSOLSPEC_H
#define __TSOLSPEC_H

#include "THaSpectrometer.h"

#include "types.h"
#include <vector>

class TSolGEMChamber;

class TSolSpec : public THaSpectrometer {
    public:
	TSolSpec( const char *name, const char *desc );
        virtual ~TSolSpec();

	Int_t AddGEM (TSolGEMChamber* pdet);

	// Useless: the actual job is done by TreeSearch.
	// However, those methods seem to have to be declared, 
	// perhaps for reasons of inheritence from THaSpectrometer
	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();

	Int_t TrackCalc() { return 0; }

	Int_t FindVertices(TClonesArray &);
	void MakePrefix(){ return; }

	UInt_t GetNChambers() const { return fChambers.size(); }
	TSolGEMChamber &GetChamber(Int_t i) const { return *(fChambers.at(i)); }

	void Print() const;

    private:
        std::vector<TSolGEMChamber*>  fChambers;

    public:
	ClassDef(TSolSpec,0)
};

#endif//__TSOLSPEC_H
