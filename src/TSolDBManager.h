#ifndef __TSOLDBMANAGER_H
#define __TSOLDBMANAGER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include "Rtypes.h"
#include "VarDef.h"
#include "TMath.h"
#include "types.h"

using namespace std;

class TSolDBManager {
public:
    ~TSolDBManager();
    static TSolDBManager* GetInstance() {
        if (fManager == NULL) fManager = new TSolDBManager();
        return fManager;
    }
    
    void LoadGeneralInfo(const string& fileName);
    void LoadGeoInfo(const string& fileName);
    
    const int    &   DoMapSector()         { return fDoMapSector; }
    const int    &   DoSelfDefineSector()  { return fDoSelfDefinedSector; }
    const int    &   GetSectorMapped()     { return fMappedSector; }
    const int    &   GetNTracker()         { return fNTracker; }
    const int    &   GetNSector()          { return fNSector; }
    const int    &   GetNReadOut()         { return fNReadOut; }
    const int    &   GetGEMDriftID()       { return fGEMDriftID; }
    const int    &   GetGEMCopperFrontID() { return fGEMCopperFrontID; }
    const int    &   GetGEMCopperBackID()  { return fGEMCopperBackID; }
    const int    &   GetGEMStripID()       { return fGEMStripID; }
    const int    &   GetNSigParticle()     { return fNSigParticle; }
    const int    &   GetFAECID()           { return fFAECID; }
    const int    &   GetLAECID()           { return fLAECID; }
    const int    &   GetChanPerSlot()      { return fChanPerSlot; }
    const int    &   GetModulesPerReadOut(){ return fModulesPerReadOut; }
    const int    &   GetModulesPerChamber(){ return fModulesPerChamber; }
    const int    &   GetChambersPerCrate() { return fChambersPerCrate; }
    
    const int    &   GetSigPID(unsigned int i);
    const int    &   GetSigTID(unsigned int i);
    
    const double &   GetSectorZ(int i, int j);
    const double &   GetSectorRMin(int i, int j);
    const double &   GetSectorRMax(int i, int j);
    const double &   GetSectorPhiStart(int i, int j);
    const double &   GetSectorPhiCover(int i, int j);
    const double &   GetSectorStripAngle(int i, int j, int k);
    const double &   GetSectorPitch(int i, int j, int k);
    
    int              GetSectorIDFromPos(double& x, double& y, int& itracker);
    
protected:
    TSolDBManager();
    int    LoadDB(ifstream& inp, DBRequest* request, const string& prefix);
    string FindKey( ifstream& inp, const string& key );
    bool   CheckIndex(int i, int j=0, int k=0);
    
    static TSolDBManager* fManager;

    //variable for data base information
    int fDoMapSector;
    int fMappedSector;
    int fDoSelfDefinedSector;
    
    int    fNTracker;
    int    fNSector;
    int    fNReadOut;
    int    fNSigParticle;
    int    fGEMDriftID;
    int    fGEMCopperFrontID;
    int    fGEMCopperBackID;
    int    fGEMStripID;
    int    fFAECID;
    int    fLAECID;
    
    int    fChanPerSlot;
    int    fModulesPerReadOut;
    int    fModulesPerChamber;
    int    fChambersPerCrate;
    
    int    fErrID;
    double fErrVal;
    
    vector<int>    fSigPID;
    vector<int>    fSigTID;
    
    vector<double> fTrackerZ;
    map< int, vector<GeoInfo> > fGeoInfo;
    
};

#endif
