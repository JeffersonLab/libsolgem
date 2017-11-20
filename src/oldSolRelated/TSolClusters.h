/////////////////////////////////////////////////////////////////////
//
//   TSolClusters
//
//   Clusters from digitized GEM hits
//
/////////////////////////////////////////////////////////////////////

#ifndef __TSolClusters_h
#define __TSolClusters_h

#include <TObject.h>
#include <vector>
#include <map>

typedef std::map < Int_t, Float_t > HitMap_t;
class TSolGEMPlane;
class TSolSimGEMDigitization;

struct Cluster_t
{
  Double_t fPos;          // coordinate in chamber
  Double_t fSize;         // width
  Double_t fCharge;       // total charge in ADC
  UInt_t fType;           // see comments in clustering code
  Double_t fResolution;   // estimated resolution
};

//-----------------------------------------------------------------------------
class TSolClusters: public TObject
{
 public: 
  TSolClusters() {Init();};
  virtual ~TSolClusters() {};

  void Init() {fRawHits.clear(); fClusters.clear();}
  void AddRawHit (UInt_t i, Float_t adc) {fRawHits[i] = adc;}

  void SetPos (UInt_t ic, Double_t vPos) {fClusters[ic].fPos = vPos;};
  void SetSize (UInt_t ic, Double_t vSize) {fClusters[ic].fSize = vSize;};
  void SetCharge (UInt_t ic, Double_t vCharge) {fClusters[ic].fCharge = vCharge; };
  void SetType (UInt_t ic, UInt_t vType) {fClusters[ic].fType = vType;};
  void SetResolution (UInt_t ic, Double_t vResolution) {fClusters[ic].fResolution = vResolution;};

  Double_t GetPos (UInt_t ic) const {return fClusters[ic].fPos;};
  Double_t GetSize (UInt_t ic) const {return fClusters[ic].fSize;};
  Double_t GetCharge (UInt_t ic) const {return fClusters[ic].fCharge;};
  UInt_t GetType (UInt_t ic) const {return fClusters[ic].fType;};
  Double_t GetResolution (UInt_t ic) const {return fClusters[ic].fResolution;};

  Int_t MakeClusters (const Double_t stripstart, 
		      const Double_t strippitch);
  Int_t ClusterPlane (TSolGEMPlane& gp, 
		      const Int_t ich,
		      const Int_t ip,
		      TSolSimGEMDigitization& ddd,
		      const Double_t cut);

  //  virtual void Print( const Option_t* opt="" ) const;

 private:

  static Double_t gSplitFrac;
  static UInt_t gMaxClusterSize;
  static UInt_t gMaxHits;
  static Double_t gBig;

  HitMap_t fRawHits;
  std::vector < Cluster_t > fClusters;

  ClassDef(TSolClusters, 1) // A true MC track
};

#endif
