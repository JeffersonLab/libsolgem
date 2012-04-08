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
#define MC_MAXC  10000
#define MC_MAXS  20000
#define MC_MAXSAMP 10

class TSolSimEvent : public TObject {
public:
  TSolSimEvent();
  virtual ~TSolSimEvent();

  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number
  Int_t     fRefFile;             // Reference file in MC

  // MC tracks
  TClonesArray*   fMCTracks;      //-> MC tracks in this event

  // Cluster variables (MC generated)
  Int_t     fNClust;              // Number of clusters (sum over all planes)
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)
  // Note: fNsignal <= fNClust (equal if no background)
  // The first fNsignal elements in these arrays are from the signal
  Short_t   fClsChamber[MC_MAXC]; // [fNClust] Chamber number of cluster
  Float_t   fClsCharge[MC_MAXC];  // [fNClust] Charge of avalanche
  Int_t     fClsRefEntry[MC_MAXC];// [fNClust] Reference entry in MC
  Float_t   fClsTime[MC_MAXC];    // [fNClust] Arrival time at electronics (w/o ToF)
  TVector3  fClsP[MC_MAXC];       // [fNClust] Momentum of particle generating the cluster
  Int_t     fClsPID[MC_MAXC];     // [fNClust] PDG ID of particle generating the cluster
  Int_t     fClsSizeX[MC_MAXC];   // [fNClust] Number of strips in cluster on x-axis
  Int_t     fClsSizeY[MC_MAXC];   // [fNClust] Number of strips in cluster on y-axis
  Int_t     fClsStartX[MC_MAXC];  // [fNClust] Number of first strip in cluster on x-axis
  Int_t     fClsStartY[MC_MAXC];  // [fNClust] Number of first strip in cluster on y-axis

  // Digitized strip amplitude data
  Int_t     fNStrips;             // Number of strips firing
  Short_t   fStpGEM[MC_MAXS];     // [fNStrips] Location index of strip
  Short_t   fStpPlane[MC_MAXS];   // [fNStrips] Plane number
  Short_t   fStpNum[MC_MAXS];     // [fNStrips] Strip number
  Short_t   fStpSigType[MC_MAXS]; // [fNStrips] Signal type (0=signal(+bck), 1=bck)
  Short_t   fStpTrack[MC_MAXS];   // [fNStrips] Track index (-1=none)
  Float_t   fStpCharge[MC_MAXS];  // [fNStrips] Total charge in strip
  Float_t   fStpTime1[MC_MAXS];   // [fNStrips] Time of first sample
                                  //   relative to event start in target (TBC)
  Short_t   fStpADC[MC_MAXSAMP][MC_MAXS]; // [fNStrips] ADC samples

  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;

  ClassDef(TSolSimEvent, 1) // Simulated data for one event
};


#endif
