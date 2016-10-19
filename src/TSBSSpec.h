#ifndef __TSBSSPEC_H
#define __TSBSSPEC_H

#include "THaSpectrometer.h"

#include "types.h"
#include <vector>

class TSBSGEMChamber;

class TSBSSpec : public THaSpectrometer {
    public:
	TSBSSpec( const char *name, const char *desc );
        virtual ~TSBSSpec();

	Int_t AddGEM (TSBSGEMChamber* pdet);

	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();

	Int_t TrackCalc() { return 0; }

	Int_t FindVertices(TClonesArray &);
//	void MakePrefix(){ return; }

	UInt_t GetNChambers() const { return fChambers.size(); }
	TSBSGEMChamber &GetChamber(Int_t i) const { return *(fChambers.at(i)); }

	void Print() const;

    private:
        std::vector<TSBSGEMChamber*>  fChambers;

    public:
	ClassDef(TSBSSpec,0)
};

#endif//__TSBSSPEC_H
