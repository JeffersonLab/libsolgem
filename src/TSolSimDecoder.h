#ifndef __TSolSimDecoder_h
#define __TSolSimDecoder_

/////////////////////////////////////////////////////////////////////
//
//   TSolSimDecoder
//
/////////////////////////////////////////////////////////////////////

#include "THaEvData.h"
#include "THaAnalysisObject.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TSolSimEvent.h"
#include <cassert>

class THaCrateMap;

//-----------------------------------------------------------------------------
// Helper classes for making decoded event data available via global variables

class TSolSimGEMHit : public TObject {
public:
  TSolSimGEMHit() {}
  TSolSimGEMHit( const TSolSimEvent::GEMCluster& cl );

  virtual void Print( const Option_t* opt="" ) const;

  Double_t P() const { return fP.Mag(); }
  Double_t X() const { return fMCpos.X(); }
  Double_t Y() const { return fMCpos.Y(); }
  Double_t Z() const { return fMCpos.Z(); }

  Int_t     fID;        // Cluster number, cross-ref to GEMStrip
  // MC hit data
  Int_t     fSector;    // Sector number
  Int_t     fPlane;     // Plane number
  Int_t     fType;      // GEANT particle counter (1 = primary)
  Int_t     fPID;       // PDG ID of particle generating the cluster
  TVector3  fP;         // Momentum of particle generating the cluster
  TVector3  fXEntry;    // Track at chamber entrance in lab coords [m]
  TVector3  fMCpos;     // Approx. truth position of hit in lab [m]
  TVector3  fHitpos;    // fMCpos in Tracker frame [m]
  // Digitization results for this hit
  Float_t   fCharge;    // Charge of avalanche
  Float_t   fTime;      // Arrival time at electronics (w/o ToF)
  Int_t     fUSize;     // Number of strips in u-cluster
  Int_t     fUStart;    // Number of first strip in u-cluster
  Float_t   fUPos;      // fMCpos along u-projection axis [m]
  Int_t     fVSize;     // Number of strips in v-cluster
  Int_t     fVStart;    // Number of first strip in v-cluster
  Float_t   fVPos;      // fMCpos along v-projection axis [m]

  ClassDef(TSolSimGEMHit, 0) // A Monte Carlo hit at a GEM tracking chamber
};

//-----------------------------------------------------------------------------
class TSolSimBackTrack : public TObject {
public:
  TSolSimBackTrack() {}
  TSolSimBackTrack( const TSolSimEvent::GEMCluster& cl );
  
  Double_t X()         const { return fOrigin.X(); }
  Double_t Y()         const { return fOrigin.Y(); }
  Double_t P()         const { return fMomentum.Mag(); }
  Double_t ThetaT()    const { return fMomentum.Px()/fMomentum.Pz(); }
  Double_t PhiT()      const { return fMomentum.Py()/fMomentum.Pz(); }
  Double_t R()         const { return fOrigin.Perp(); }
  Double_t Theta()     const { return fOrigin.Theta(); }
  Double_t Phi()       const { return fOrigin.Phi(); }
  Double_t ThetaDir()  const { return fMomentum.Theta(); }
  Double_t PhiDir()    const { return fMomentum.Phi(); }
  Double_t HX()        const { return fHitpos.X(); }
  Double_t HY()        const { return fHitpos.Y(); }

  virtual void Print( const Option_t* opt="" ) const;

  Int_t    GetType()    const { return fType; }
  Int_t    GetHitBits() const { return fHitBits; }
  void     SetHitBit( UInt_t i )  { SETBIT(fHitBits,i); }
  Bool_t   TestHitBit( UInt_t i ) { return TESTBIT(fHitBits,i); }

  Int_t    Update( const TSolSimEvent::GEMCluster& cl );

private:

  Int_t    fType;        // GEANT particle number
  Int_t    fPID;         // Track particle ID (PDG)
  Int_t    fSector;      // Sector where this track detected
  TVector3 fOrigin;      // Position at first plane in lab frame (m)
  TVector3 fHitpos;      // Position at first plane in Tracker frame [m]
  TVector3 fMomentum;    // Momentum (GeV)
  Int_t    fHitBits;     // Bitpattern of plane numbers hit by this back track

  ClassDef(TSolSimBackTrack, 0) // MC track registered at first tracking chamber
};

//-----------------------------------------------------------------------------
// SoLID simulation decode class
class TSolSimDecoder : public THaEvData {
 public:
  TSolSimDecoder();
  virtual ~TSolSimDecoder();

  virtual Int_t LoadEvent( const int* evbuffer, THaCrateMap* usermap );
  virtual void  Clear( Option_t* opt="" );

  Int_t  GetNBackTracks() const { return fBackTracks->GetLast()+1; }
  Int_t  GetNHits()       const { return fHits->GetLast()+1; }
  Int_t  GetNTracks()     const { return fTracks.GetSize(); }

  Int_t  DefineVariables( THaAnalysisObject::EMode mode =
			  THaAnalysisObject::kDefine );

  TSolSimBackTrack* GetBackTrack( Int_t i ) const {
    TObject* obj = fBackTracks->UncheckedAt(i);
    assert( dynamic_cast<TSolSimBackTrack*>(obj) );
    return static_cast<TSolSimBackTrack*>(obj);
  }
  TSolSimGEMHit* GetGEMHit( Int_t i ) const {
    TObject* obj = fHits->UncheckedAt(i);
    assert( dynamic_cast<TSolSimGEMHit*>(obj) );
    return static_cast<TSolSimGEMHit*>(obj);
  }

 protected:

  TList          fTracks;     // MC physics tracks
  TClonesArray*  fHits;       //-> MC hits (clusters)
  TClonesArray*  fBackTracks; //-> Primary particle tracks at first chamber

  Bool_t  fIsSetup;

  Int_t DoLoadEvent( const int* evbuffer, THaCrateMap* usermap );

  ClassDef(TSolSimDecoder,0) // Decoder for simulated SoLID spectrometer data
};


#endif
