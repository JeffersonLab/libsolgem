#ifndef __TSBSSPEC_H
#define __TSBSSPEC_H

#include "THaSpectrometer.h"

#include "types.h"
#include <vector>

class TSBSGEMChamber;

// Class TSBS Spec is more or less a "container" to store the information of all GEM chambers.
// Ultimately, it will also contain information on reconstrcuted tracks, 
// but that needs more work.

class TSBSSpec : public THaSpectrometer {
    public:
        //Constructor and destructor
	TSBSSpec( const char *name, const char *desc );
        virtual ~TSBSSpec();

	Int_t AddGEM (TSBSGEMChamber* pdet);

	// Functions to reconstruct tracks, etc. Still void (return 0) 
	// TO-DO (???): develop these
	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();

	Int_t TrackCalc() { return 0; }

	Int_t FindVertices(TClonesArray &);
//	void MakePrefix(){ return; }
	
	//Access to GEM chambers info 
	UInt_t GetNChambers() const { return fChambers.size(); }
	TSBSGEMChamber &GetChamber(Int_t i) const { return *(fChambers.at(i)); }
	
	//Print spectrometer info, with each individual GEM chamber
	void Print() const;

    private:
        std::vector<TSBSGEMChamber*>  fChambers;

    public:
	ClassDef(TSBSSpec,0)
};

#endif//__TSBSSPEC_H
