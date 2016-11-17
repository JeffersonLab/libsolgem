#ifndef __TSBSSPEC_H
#define __TSBSSPEC_H

#include "THaSpectrometer.h"

#include "types.h"
#include <vector>

class TSBSGEMChamber;

// Class TSBSSpec is more or less a "container" to store the information of all GEM chambers.
// Ultimately, it will also contain information on reconstrcuted tracks, 
// but that needs more work.
// This class inherits from class THaSpectrometer, 
// which grants it all the functions from its class
// (see http://hallaweb.jlab.org/podd/doc/html_v16/ClassIndex.html for more info).

class TSBSSpec : public THaSpectrometer {
    public:
        //Constructor and destructor
	TSBSSpec( const char *name, const char *desc );
        virtual ~TSBSSpec();

	Int_t AddGEM (TSBSGEMChamber* pdet);

	// Useless: the actual job is done by TreeSearch.
	// However, those methods seem to have to be declared, 
	// perhaps for reasons of inheritence from THaSpectrometer
	Int_t CoarseTrack();
	Int_t CoarseReconstruct();
	Int_t Track();
	Int_t Reconstruct();
	
	Int_t TrackCalc() { return 0; }
	
	Int_t FindVertices(TClonesArray &);
	/* void MakePrefix(){ return; } */
	
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
