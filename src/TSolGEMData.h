// manage the input data (output of the GEANT4 simulation)
// for the gem part

#ifndef __TSolGEMData__
#define __TSolGEMData__

#include <TRandom.h>
#include <TVector3.h>

class TSolGEMData 
{
 public:

  TSolGEMData (UInt_t h = 0);
  virtual ~TSolGEMData();

  void ClearHit();
  void InitHit (UInt_t h);
  
  void SetNHit (UInt_t h) {fNHit = h;};
  void SetEvent (UInt_t id) {fEvtID = id;};
  void SetRun (UInt_t r) {fRunID = r;};
  void SetMomentum (UInt_t k, const TVector3 p);
  void SetHitEntrance (UInt_t k, const TVector3 xEnt);
  void SetHitExit (UInt_t k, const TVector3 xExit);
  void SetHitReadout (UInt_t k, const TVector3 xRead);
  void SetHitEnergy (UInt_t k, Double_t e) {fEdep[k] = e;};
  void SetHitChamber (UInt_t k, UInt_t c) {fGem[k] = c;};
  void SetEntryNumber (UInt_t k, UInt_t eno) {fIdxV[k] = eno;};
  void SetParticleID (UInt_t k, UInt_t pid) {fPID[k] = pid;};
  void SetParticleType (UInt_t k, UInt_t type) {fType[k] = type;};

  UInt_t GetNHit() const {return fNHit;};
  UInt_t GetEvent() const {return fEvtID;};
  UInt_t GetRun() const {return fRunID;};
  TVector3 *GetMomentum (UInt_t k) const {return fMom[k];};
  TVector3 *GetHitEntrance (UInt_t k) const {return fXi[k];};
  TVector3 *GetHitExit (UInt_t k) const {return fXo[k];};
  TVector3 *GetHitReadout (UInt_t k) const {return fXr[k];};
  Double_t GetHitEnergy (UInt_t k) const {return fEdep[k];};
  UInt_t GetHitChamber (UInt_t k) const {return fGem[k];};
  UInt_t GetEntryNumber (UInt_t k) const {return fIdxV[k];};
  UInt_t GetParticleID (UInt_t k) const {return fPID[k];};
  UInt_t GetParticleType (UInt_t k) const {return fType[k];};

 private:

  UInt_t fNHit; // number of hits in event

  UInt_t fRunID, fEvtID;

  // Hits

  UInt_t *fGem; // chamber with hit
  Double_t *fEdep; // energy lost in drift
  UInt_t *fIdxV; // entry number
  UInt_t *fPID; // particle ID
  UInt_t *fType;  // 0 for signal, >0 for backgrounds
  TVector3 **fXi; // entrance point in drift
  TVector3 **fXo; // exit point in drift
  TVector3 **fXr; // entrance point in readout
  TVector3 **fMom; // momentum of the particle

};
#endif
