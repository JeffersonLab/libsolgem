#ifndef __TSOLGEMPLANE_H
#define __TSOLGEMPLANE_H

#include "THaDetector.h"
#include "THaEvData.h"

#include "types.h"

class TSolGEMCluster;
class TClonesArray;

class TSolGEMPlane : public THaDetector {
    public:
	TSolGEMPlane(const char *name, const char *desc);
	virtual ~TSolGEMPlane() {;}

	TClonesArray *GetClusters() { return fClusters; }

	Int_t Decode( const THaEvData &);
	GEMDir_t GetDirection(){ return fDir; }
	TSolGEMPlane *GetPairedPlane() { return fPairPlane; }


    private:
	TClonesArray  *fClusters; // Clusters
	GEMDir_t fDir;		 // Plane orientation
	TSolGEMPlane *fPairPlane; // Paired plane

    public:
	ClassDef(TSolGEMPlane,1)

};

#endif//__TSOLGEMPLANE_H
