/////////////////////////////////////////////////////////////////////
//
//   TSolSimEvent
//
//   Common class definitions (event, MC track, etc.) for SoLID
//   simulation decoder
//
//   Ole Hansen (ole@jlab.org)  December 2011
//
/////////////////////////////////////////////////////////////////////

#ifndef __TSolSim_h
#define __TSolSim_h

#include "TVector3.h"
#include <vector>

class TClonesArray;

//-----------------------------------------------------------------------------
class TSolSimTrack : public TObject {
 public: 
  TSolSimTrack() : fType(0), fLayer(0), fPID(0), fTimeOffset(0.0), fHit1(-1) {}
  TSolSimTrack( Int_t type, Int_t layer, Int_t pid, Double_t toff )
    : fType(type), fLayer(layer), fPID(pid), fTimeOffset(toff), fHit1(-1) {}

  Int_t    fType;        // Track type
  Int_t    fLayer;       // Material type where track originated (unused)
  Int_t    fPID;         // Track particle ID (PDG)
  Double_t fTimeOffset;  // Time offset of the track (s)
  TVector3 fOrigin;      // Origin (m)
  TVector3 fMomentum;    // Momentum (GeV)
  Int_t    fHit1;        // Chamber where this track makes its first hit

  Double_t TX()     { return fOrigin.X(); }
  Double_t TY()     { return fOrigin.Y(); }
  Double_t TTheta() { return fMomentum.Px()/fMomentum.Pz(); }
  Double_t TPhi()   { return fMomentum.Py()/fMomentum.Pz(); }
  Double_t P()      { return fMomentum.Mag(); }

  virtual void Print( const Option_t* opt="" ) const;
  
  ClassDef(TSolSimTrack, 1) // A true MC track
};

//-----------------------------------------------------------------------------
// Kludgy hardcoded parameters necessary because I can't get ROOT to allocate
// arrays dynamically via TTree::GetEntry
//#define MC_MAXC  10000
//#define MC_MAXS  20000
#define MC_MAXSAMP 10

#define treeName "digtree"
#define eventBranchName "event"

class TSolSimEvent : public TObject {
public:
  TSolSimEvent();                 // Default constructor, for ROOT I/O
  TSolSimEvent( UInt_t ntracks ); // Construct and initialize fMCTracks
  virtual ~TSolSimEvent();

  UInt_t NGEMClust () const {return fGEMClust.size();};
  Short_t GetCChamber (UInt_t ic) const {return fGEMClust[ic].fChamber;};
  Float_t GetCCharge (UInt_t ic) const {return fGEMClust[ic].fCharge;};
  Int_t GetCRefEntry (UInt_t ic) const {return fGEMClust[ic].fRefEntry;};
  Int_t GetCRefFile (UInt_t ic) const {return fGEMClust[ic].fRefFile;};
  Float_t GetCTime (UInt_t ic) const {return fGEMClust[ic].fTime;};
  TVector3 GetCP (UInt_t ic) const {return fGEMClust[ic].fP;};
  Int_t GetCPID (UInt_t ic) const {return fGEMClust[ic].fPID;};
  Int_t fGetCSize (UInt_t ic, UInt_t i) const {return fGEMClust[ic].fSize[i];};
  Int_t fGetCStart (UInt_t ic, UInt_t i) const {return fGEMClust[ic].fStart[i];};
  TVector3 GetCXEntry (UInt_t ic) const {return fGEMClust[ic].fXEntry;};

  UInt_t NGEMStrips () const {return fGEMStrips.size();};
  Short_t GetSGEM (UInt_t i) const {return fGEMStrips[i].fGEM;};
  Short_t GetSPlane (UInt_t i) const {return fGEMStrips[i].fPlane;};
  Short_t GetSNum (UInt_t i) const {return fGEMStrips[i].fNum;};
  Short_t GetSSigType (UInt_t i) const {return fGEMStrips[i].fSigType;};
  Short_t GetSTrack (UInt_t i) const {return fGEMStrips[i].fTrack;};
  Float_t GetSCharge (UInt_t i) const {return fGEMStrips[i].fCharge;};
  Float_t GetSTime1 (UInt_t i) const {return fGEMStrips[i].fTime1;};
  Short_t GetSADC (UInt_t i, UInt_t isamp) const {return fGEMStrips[i].fADC[isamp];};

  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  // MC tracks
  TClonesArray*   fMCTracks;      //-> MC tracks in this event
  
  // Cluster variables (MC generated)
  // Note: fNSignal <= fClust.size() (equal if no background)
  // The first fNSignal elements in the array of clusters are from the signal
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)

  struct GEMCluster {
    Short_t   fChamber;   // Chamber number of cluster
    Float_t   fCharge;    // Charge of avalanche
    Int_t     fRefEntry;  // Reference entry in MC
    Int_t     fRefFile;   // Particle type associated with hit
    Float_t   fTime;      // Arrival time at electronics (w/o ToF)
    TVector3  fP;         // Momentum of particle generating the cluster
    Int_t     fPID;       // PDG ID of particle generating the cluster
    Int_t     fSize[2];   // Number of strips in cluster per axis
    Int_t     fStart[2];  // Number of first strip in cluster per axis
    TVector3  fXEntry;    // Truth position of hit
  };

  std::vector<GEMCluster> fGEMClust;  // All MC-generated clusters in the GEMs
  
  // Digitized strip amplitude data
  struct DigiGEMStrip {
    Short_t   fGEM;       // Location index of strip
    Short_t   fPlane;     // Plane number
    Short_t   fNum;       // Strip number
    Short_t   fSigType;   // Signal type (0=signal(+bck), 1=bck)
    Short_t   fTrack;     // Track index (-1=none)
    Float_t   fCharge;    // Total charge in strip
    Float_t   fTime1;     // Time of first sample
                             //   relative to event start in target (TBC)
    Short_t   fADC[MC_MAXSAMP]; // ADC samples
  };

  std::vector<DigiGEMStrip> fGEMStrips; // Digitized strip amplitudes of the GEMs

  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ClassDef(TSolSimEvent, 1) // Simulated data for one event
};


#endif
