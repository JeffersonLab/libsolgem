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
    void LoadGeoInfo(const string& prefix);
    
    const int    &   DoMapSector()          { return fDoMapSector;         }
    const int    &   DoSelfDefineSector()   { return fDoSelfDefinedSector; }
    const int    &   GetSectorMapped()      { return fMappedSector;        }
    const int    &   GetNChamber()          { return fNChamber;            }
    const int    &   GetNSector()           { return fNSector;             }
    /* // see comment l. 93-95. */
    /* const int    &   GetNChamber2()          { return fNChamber2;            } */
    /* const int    &   GetNSector2()           { return fNSector2;             } */
    const int    &   GetNReadOut()          { return fNReadOut;            }
    const int    &   GetGEMDriftID()        { return fGEMDriftID;          }
    const int    &   GetGEMCopperFrontID()  { return fGEMCopperFrontID;    }
    const int    &   GetGEMCopperBackID()   { return fGEMCopperBackID;     }
    const int    &   GetGEMStripID()        { return fGEMStripID;          }
    const int    &   GetNSigParticle()      { return fNSigParticle;        }
    const int    &   GetFAECID()            { return fFAECID;              }
    const int    &   GetLAECID()            { return fLAECID;              }
    const int    &   GetChanPerSlot()       { return fChanPerSlot;         }
    const int    &   GetModulesPerReadOut() { return fModulesPerReadOut;   }
    const int    &   GetModulesPerChamber() { return fModulesPerChamber;   }
    const int    &   GetChambersPerCrate()  { return fChambersPerCrate;    }
    
    const int    &   GetSigPID(unsigned int i);
    const int    &   GetSigTID(unsigned int i);
    
    const int    &   Getg4sbsDetectorType() { return fg4sbsDetectorType;   }
    const double &   Getg4sbsZSpecOffset()  { return fg4sbsZSpecOffset;    }
    
    const double &   GetZ0()                { return fgZ0;                 }
    const double &   GetCaloZ()             { return fgCaloZ;              }
    const double &   GetCaloRes()           { return fgCaloRes;            }
    const int    &   DoCalo()               { return fgDoCalo;             }

    void     SetZ0( Double_t z0 ) { fgZ0 = z0; }
    // Support for calorimeter emulation. Static functions to allow script access
    void     EmulateCalorimeter( Bool_t f = true ) { fgDoCalo = f; }
    void     SetCaloZ( Double_t z )     { fgCaloZ   = z; }
    void     SetCaloRes( Double_t res ) { fgCaloRes = res; }
    
    const double &   GetDMag(int i, int j);
    const double &   GetD0(int i, int j);
    const double &   GetXOffset(int i, int j);
    const double &   GetDX(int i, int j);
    const double &   GetDY(int i, int j);
    const double &   GetThetaH(int i, int j);
    const double &   GetThetaV(int i, int j);
    const double &   GetStripAngle(int i, int j, int k);
    const double &   GetPitch(int i, int j, int k);
    
    int GetSectorIDFromPos(int ichamber, double x, double y = 0);
    
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
    
    int    fNChamber;
    int    fNSector;
    /* // the two following lines have to be added for BB GEMs, */
    /* // since there are two types of trackers which it is not relevant to treat independently */
    /* // (in contrast with FT/FPP). */
    /* int    fNChamber2; */
    /* int    fNSector2; */
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
    
    // Parameters for TSBSGeant4File
    int fg4sbsDetectorType;// flag to determine which type of GEM should be read.
    //Choices are: 1 - BB GEMs
    //Choices are: 2 - SIDIS SBS GEMs
    //Choices are: 3 - GEP SBS GEMs: FT
    //Choices are: 3 - GEP SBS GEMs: FPP
    double fg4sbsZSpecOffset;
    //Offset between the local z value recorded in g4sbs and the actual distance 
    //of the GEM from the midplane pivot
    
    // Parameters for TSBSSimDecoder
    // Calorimeter emulation
    double fgCaloZ;     // z position of emulated calorimeter
    double fgCaloRes;   // Resolution (sigma) of emulated calorimeter (m)
    int    fgDoCalo;    // Enable calorimeter emulation
    double fgZ0;        // z position of first chamber plane
    
    int    fErrID;
    double fErrVal;
    
    vector<int>    fSigPID;
    vector<int>    fSigTID;
    
    /* vector<double> fChamberZ; */
    map< int, vector<GeoInfo> > fGeoInfo;
    
};

#endif
