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

#include "SimDecoder.h"
#include "TVector3.h"
#include "TArrayS.h"
#include <vector>

class TClonesArray;

//-----------------------------------------------------------------------------
class TSolSimTrack : public Podd::MCTrack {
public:
  TSolSimTrack( Int_t number, Int_t pid,
		const TVector3& vertex, const TVector3& momentum );
  TSolSimTrack();

  ClassDef(TSolSimTrack,3)  // A MC physics track in SoLID
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

  virtual void Clear( const Option_t* opt="" );
  virtual void Print( const Option_t* opt="" ) const;
  TSolSimTrack* AddTrack( Int_t number, Int_t pid,
			  const TVector3& vertex, const TVector3& momentum );

  Int_t GetNclust()  const { return fGEMClust.size(); }
  Int_t GetNstrips() const { return fGEMStrips.size(); }
  Int_t GetNtracks() const;

  // Event identification
  Int_t     fRunID;               // Run number
  Int_t     fEvtID;               // Event number

  Double_t  fWeight;              // Event weight

  // MC tracks
  TClonesArray*   fMCTracks;      //-> Physics tracks

  // Cluster variables (MC generated)
  // Each "cluster" is a MC hit along with the digitized charge avalance data
  // Note: fNSignal <= fClust.size() (equal if no background)
  // The first fNSignal elements in the array of clusters are from the signal
  Int_t     fNSignal;             // Number of clusters from trigger track (signal)

  Bool_t    fSectorsMapped;       // Sectors are mapped to 0
  Int_t     fSignalSector;        // Sector of primary signal particle

  struct GEMCluster {
    Short_t   fID;        // Cluster number, cross-ref to GEMStrip
    // MC hit data
    Short_t   fSector;    // Sector number
    Short_t   fPlane;     // Plane number
    Short_t   fRealSector;// Real sector number (may be !=fSector if mapped)
    Int_t     fSource;    // MC data set source (0 = signal, >0 background)
    Int_t     fType;      // GEANT particle counter (1 = primary)
    Int_t     fPID;       // PDG ID of particle generating the cluster
    TVector3  fP;         // Momentum of particle generating the cluster [GeV]
    TVector3  fXEntry;    // Track at chamber entrance in lab coords [m]
    TVector3  fMCpos;     // Approx. truth position of hit in lab [m]
    TVector3  fHitpos;    // fMCpos in Tracker frame [m]
    // Digitization results for this hit
    Float_t   fCharge;    // Charge of avalanche
    Float_t   fTime;      // Arrival time at electronics
    Int_t     fSize[2];   // Number of strips in cluster per axis
    Int_t     fStart[2];  // Number of first strip in cluster per axis
    Float_t   fXProj[2];  // fMCpos along projection axis [m]
    TVector3  fVertex;    // Vertex
  };

  std::vector<GEMCluster> fGEMClust;  // All MC-generated clusters in the GEMs

  // Digitized strip amplitude data
  struct DigiGEMStrip {
    Short_t   fSector;    // Sector number
    Short_t   fPlane;     // Plane number
    Short_t   fProj;      // Readout coordinate ("x" = 0, "y" = 1)
    Short_t   fChan;      // Channel number
    Short_t   fSigType;   // Accumulated signal types (BIT(0) = signal)
    Float_t   fCharge;    // Total charge in strip
    Float_t   fTime1;     // Time of first sample
                          //   relative to event start in target (TBC)
    UShort_t  fNsamp;     // Number of ADC samples
    Int_t     fADC[MC_MAXSAMP]; // [fNsamp] ADC samples
    TArrayS   fClusters;  // Clusters contributing to signal on this strip
  };

  std::vector<DigiGEMStrip> fGEMStrips; // Digitized strip amplitudes of the GEMs

  ClassDef(TSolSimEvent, 5) // Simulated data for one event
};

#endif
