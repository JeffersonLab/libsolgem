// manage the input data (output of the GEANT4 simulation)
// for the gem part

#ifndef __TSolGEMData__
#define __TSolGEMData__

#include <TRandom.h>
#include <TVector3.h>
#include <vector>

class TSolGEMData 
{
 public:

  TSolGEMData (UInt_t h = 0);
  virtual ~TSolGEMData();

  void ClearEvent();
  void InitEvent (UInt_t h);
  
  void SetNHit (UInt_t h)   { fHitData.resize(h); }
  void SetEvent (UInt_t id) { fEvtID = id; }
  void SetRun (UInt_t r)    { fRunID = r; }
  void SetSource (Int_t s)  { fSource = s; }

  void SetMomentum (UInt_t k, const TVector3& p) { fHitData[k].fMom = p; }
  // Positions are in mm
  void SetHitEntrance (UInt_t k, const TVector3& xEnt) { fHitData[k].fXi = xEnt; }
  void SetHitExit (UInt_t k, const TVector3& xExit)    { fHitData[k].fXo = xExit; }
  void SetHitReadout (UInt_t k, const TVector3& xRead) { fHitData[k].fXr = xRead; }
  void SetVertex (UInt_t k, const TVector3& v)         { fHitData[k].fVert = v; }
  // Energy lost is in eV
  void SetHitEnergy (UInt_t k, Double_t e)     { fHitData[k].fEdep = e; }
  void SetHitTime(UInt_t k, Double_t t)        { fHitData[k].fTime = t; }
  void SetHitChamber (UInt_t k, UInt_t c)      { fHitData[k].fGem  = c; }
  void SetParticleID (UInt_t k, Int_t pid)     { fHitData[k].fPID  = pid; }
  void SetParticleType (UInt_t k, UInt_t type) { fHitData[k].fType = type; }

  UInt_t GetNHit()   const { return fHitData.size(); }
  UInt_t GetEvent()  const { return fEvtID; }
  UInt_t GetRun()    const { return fRunID; }
  Int_t  GetSource() const { return fSource; }

  const TVector3& GetMomentum (UInt_t k)    const { return fHitData[k].fMom; }
  const TVector3& GetHitEntrance (UInt_t k) const { return fHitData[k].fXi; }
  const TVector3& GetHitExit (UInt_t k)     const { return fHitData[k].fXo; }
  const TVector3& GetHitReadout (UInt_t k)  const { return fHitData[k].fXr; }
  const TVector3& GetVertex (UInt_t k)      const { return fHitData[k].fVert; }
  Double_t GetHitEnergy (UInt_t k)    const { return fHitData[k].fEdep; }
  Double_t GetHitTime(UInt_t k)       const { return fHitData[k].fTime; }
  UInt_t GetHitChamber (UInt_t k)     const { return fHitData[k].fGem; }
  Int_t  GetParticleID (Int_t k)      const { return fHitData[k].fPID; }
  UInt_t GetParticleType (UInt_t k)   const { return fHitData[k].fType; }

  void Print() const;
  void PrintHit (UInt_t k) const;

 private:

  UInt_t fRunID, fEvtID;
  Int_t  fSource; // MC source file ID (0 = signal, >0 background)

  // Hit data
  struct GEMHitData_t {
    UInt_t    fGem;
    Double_t  fEdep;  // energy lost in drift
    Double_t  fTime;  // hit time
    Int_t     fPID;   // particle ID
    UInt_t    fType;  // particle counter: 1 = primary, >1 secondaries
    TVector3  fXi;    // entrance point in drift
    TVector3  fXo;    // exit point in drift
    TVector3  fXr;    // entrance point in readout
    TVector3  fMom;   // momentum of the particle
    TVector3  fVert;  // vertex position
  };
  std::vector<GEMHitData_t> fHitData;

};
#endif
