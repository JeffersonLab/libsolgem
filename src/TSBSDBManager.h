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

class TSBSDBManager {
public:
    ~TSBSDBManager();
    static TSBSDBManager* GetInstance() {
        if (fManager == NULL) fManager = new TSBSDBManager();
        return fManager;
    }
    
    void LoadGeneralInfo(const string& fileName);
    void LoadGeoInfo(const string& prefix);
   
    int       DoMapSector() const          { return fDoMapSector;         }
    int       DoSelfDefineSector() const   { return fDoSelfDefinedSector; }
    int       GetSectorMapped() const      { return fMappedSector;        }
    int       GetNChamber() const          { return fNChamber;            }
    int       GetNSector() const           { return fNSector;             }
    int       GetNGEMPlane() const         { return fNGEMPlane;           }
    int       GetNModule(int plane) const  { if(plane<0) return 0; else return fNModule[plane]; }
    int       GetPlaneID(int igem)         { return fmIgemtoPlane[igem];  }
    int       GetModuleID(int igem)        { return fmIgemtoModule[igem]; }
    int       GetGEMID(int ip,int im)      { return fmPMtoIgem[ip][im];   }
    /* // see comment l. 93-95. */
    /* int      GetNChamber2() const          { return fNChamber2;            } */
    /* int      GetNSector2() const           { return fNSector2;             } */
    int       GetNReadOut() const          { return fNReadOut;            }
    int       GetGEMDriftID() const        { return fGEMDriftID;          }
    int       GetGEMCopperFrontID() const  { return fGEMCopperFrontID;    }
    int       GetGEMCopperBackID() const   { return fGEMCopperBackID;     }
    int       GetGEMStripID() const        { return fGEMStripID;          }
    int       GetNSigParticle() const      { return fNSigParticle;        }
    int       GetFAECID() const            { return fFAECID;              }
    int       GetLAECID() const            { return fLAECID;              }
    int       GetChanPerSlot() const       { return fChanPerSlot;         }
    int       GetModulesPerReadOut() const { return fModulesPerReadOut;   }
    int       GetModulesPerChamber() const { return fModulesPerChamber;   }
    int       GetChambersPerCrate() const  { return fChambersPerCrate;    }
    
    int       GetSigPID(unsigned int i) const;
    int       GetSigTID(unsigned int i) const;
    
    int       Getg4sbsDetectorType() const { return fg4sbsDetectorType;   }
    double    Getg4sbsZSpecOffset() const  { return fg4sbsZSpecOffset;    }
    
    double    GetCaloThreshold() const     { return fCaloThr;             }
    double    GetZ0() const                { return fgZ0;                 }
    double    GetCaloZ() const             { return fgCaloZ;              }
    double    GetCaloRes() const           { return fgCaloRes;            }
    int       DoCalo() const               { return fgDoCalo;             }

    void     SetZ0( Double_t z0 ) { fgZ0 = z0; }
    // Support for calorimeter emulation. Static functions to allow script access
    void     EmulateCalorimeter( Bool_t f = true ) { fgDoCalo = f; }
    void     SetCaloZ( Double_t z )     { fgCaloZ   = z; }
    void     SetCaloRes( Double_t res ) { fgCaloRes = res; }
    
    double    GetDMag(int i, int j);
    double    GetD0(int i, int j);
    double    GetXOffset(int i, int j);
    double    GetDX(int i, int j);
    double    GetDY(int i, int j);
    //double    GetThetaH(int i, int j);
    double    GetThetaV(int i, int j);
    double    GetStripAngle(int i, int j, int k);
    //double    GetPitch(int i, int j, int k);
    
    int GetModuleIDFromPos(int iplane, double x, double y = 0);
    double GetPosFromModuleStrip(int iproj, int iplane, int isector, int istrip);

protected:
    TSBSDBManager();
    int    LoadDB(ifstream& inp, DBRequest* request, const string& prefix);
    string FindKey( ifstream& inp, const string& key ) const;
    bool   CheckIndex(int i, int j=0, int k=0) const;
    
    static TSBSDBManager* fManager;

    //variable for data base information
    int fDoMapSector;
    int fMappedSector;
    int fDoSelfDefinedSector;
    
    int    fNChamber;
    int    fNSector;
    int    fNGEMPlane;
    std::vector<Int_t> fNModule;
    std::map<Int_t,std::map<Int_t, Int_t> > fmPMtoIgem;
    std::map<Int_t, Int_t> fmIgemtoPlane;
    std::map<Int_t, Int_t> fmIgemtoModule;

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

    // Parameters for simple calorimeter analysis
    double fCaloThr;
    
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
    map< int, vector<GeoInfo> > fPMGeoInfo; //plane module format geo info
    
};

#endif
