// manage the input data (output of the GEANT4 simulation)
// for the gem part

#ifndef __TSolGEMData__
#define __TSolGEMData__

#include <TRandom.h>
#include <TVector3.h>

class TSolGEMData 
{
 public:

  TSolGEMData (Int_t h = 0);
  virtual ~TSolGEMData();

  void ClearHit();
  void InitHit (Int_t h);
  
  void SetNHit (Int_t h) {fNHit = h;};
  void SetEvent (Int_t id) {fEvtID = id;};
  void SetRun (Int_t r) {fRunID = r;};
  void SetMomentum (Int_t k, const TVector3 p);
  void SetHitEntrance (Int_t k, const TVector3 xEnt);
  void SetHitExit (Int_t k, const TVector3 xExit);
  void SetHitReadout (Int_t k, const TVector3 xRead);
  void SetHitEnergy (Int_t k, Double_t e) {fEdep[k] = e;};
  void SetHitChamber (Int_t k, Int_t c) {fGem[k] = c;};
  void SetEntryNumber (Int_t k, Int_t eno) {fIdxV[k] = eno;};
  void SetParticleID (Int_t k, Int_t pid) {fPID[k] = pid;};

  Int_t GetNHit() const {return fNHit;};
  Int_t GetEvent() const {return fEvtID;};
  Int_t GetRun() const {return fRunID;};
  TVector3 *GetMomentum (Int_t k) const {return fMom[k];};
  TVector3 *GetHitEntrance (Int_t k) const {return fXi[k];};
  TVector3 *GetHitExit (Int_t k) const {return fXo[k];};
  TVector3 *GetHitReadout (Int_t k) const {return fXr[k];};
  Double_t GetHitEnergy (Int_t k) const {return fEdep[k];};
  Int_t GetHitChamber (Int_t k) const {return fGem[k];};
  Int_t GetEntryNumber (Int_t k) const {return fIdxV[k];};
  Int_t GetParticleID (Int_t k) const {return fPID[k];};

 private:

  Int_t fNHit; // number of hits in event

  Int_t fRunID, fEvtID;

  // Hits

  Int_t *fGem; // chamber with hit
  Double_t *fEdep; // energy lost in drift
  Int_t *fIdxV; // entry number
  Int_t *fPID; // particle ID
  TVector3 **fXi; // entrance point in drift
  TVector3 **fXo; // exit point in drift
  TVector3 **fXr; // entrance point in readout
  TVector3 **fMom; // momentum of the particle

};
#endif
