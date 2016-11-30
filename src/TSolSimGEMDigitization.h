
// Class handling digitization of GEMs

#ifndef __TSolSimGEMDigitization__
#define __TSolSimGEMDigitization__

#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayI.h"

#include "THaAnalysisObject.h"

#include <vector>

class TFile;
class TTree;

class TSolGEMData;
class TSolGEMVStrip;
class TSolSpec;
class TSolSimEvent;
class TSolEVIOFile;

// First an auxiliary class

// The whole strip plane; used to cumulate virtual strips charges
// and generate real strips

class TSolDigitizedPlane {
private:
  // ADC sampled value of strip array of each axis

  //TODO: make this a struct inside an STL vector or similar
  TArrayI fStripADC;
  Short_t *fType;
  Int_t   *fTotADC;

  Float_t *fCharge;
  Float_t *fTime;

  UShort_t fNSamples;
  UShort_t fNStrips;
  Int_t    fThreshold;

  UShort_t  fNOT;   // # strips over threshold
  Short_t*  fOverThr;  // # list of strips over threshold

  std::vector< std::vector<Short_t> > fStripClusters; // Clusters associated with each strip

public:
  // init and reset physics strips arrays
  TSolDigitizedPlane (UShort_t nstrip,
		      UShort_t nsample = 10,
		      Int_t    threshold = 0 );
  ~TSolDigitizedPlane();
  void Clear();

  // cumulate hits (strips signals)
  void Cumulate (const TSolGEMVStrip *vv, Short_t type, Short_t clusterID );

  Short_t  GetType (Int_t n) const {return fType[n];}
  Int_t    GetTotADC (Int_t n) const {return fTotADC[n];}
  Float_t  GetTime (Int_t n) const {return fTime[n];}
  Float_t  GetCharge (Int_t n) const {return fCharge[n];}
  Int_t    GetADC (Int_t n, Int_t ks) const {return fStripADC[n*fNSamples+ks];}
  UShort_t GetNSamples() const {return fNSamples;}
  UShort_t GetNStrips() const {return fNStrips;}

  UShort_t Threshold (Int_t thr);

  UShort_t GetNOverThr() const {return fNOT;}
  Short_t  GetIdxOverThr (Int_t n) const {return fOverThr[n];}

  const std::vector<Short_t>& GetStripClusters(UInt_t n) const { return fStripClusters[n]; }
};



class TSolSimGEMDigitization: public THaAnalysisObject
{
 public:
  TSolSimGEMDigitization( const TSolSpec& spect,
			  const char* name = "testdigitization");
  virtual ~TSolSimGEMDigitization();

  void Initialize(const TSolSpec& spect);
  Int_t ReadDatabase (const TDatime& date);

  Int_t Digitize (const TSolGEMData& gdata, const TSolSpec& spect); // digitize event
  Int_t AdditiveDigitize (const TSolGEMData& gdata, const TSolSpec& spect); // add more digitized data, e.g. background
  void NoDigitize (const TSolGEMData& gdata, const TSolSpec& spect); // do not digitize event, just fill tree
  const TSolDigitizedPlane& GetDigitizedPlane (UInt_t ich, UInt_t ip) const {return *(fDP[ich][ip]);};
  void Print() const;
  void PrintCharges() const;
  void PrintSamples() const;

  Double_t GetGateWidth(){ return fGateWidth; }

  // Tree methods
  // To write a tree with digitization results:
  //   Call InitTree before main loop;
  //   Call SetTreeEvent in main loop (before or after Digitize)
  //   Call FillTree in main loop (after Digitize and SetTreeEvent)
  // Call WriteTree and CloseTree after main loop

  void InitTree (const TSolSpec& spect, const TString& ofile);
  void SetTreeEvent (const TSolGEMData& tsgd,
		     const TSolEVIOFile& f,
		     Int_t evnum = -1);
  Short_t SetTreeHit (const UInt_t ih,
		      const TSolSpec& spect,
		      TSolGEMVStrip* const *dh,
		      const TSolGEMData& tsgd,
		      Double_t t0 ); // called from Digitization
  void SetTreeStrips(); // called from Digitization
  void FillTree ();
  void WriteTree () const;
  void CloseTree () const;

  // Access to results
  Short_t GetType (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetType (n);}
  Int_t   GetTotADC (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetTotADC (n);}
  Float_t GetTime (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetTime (n);}
  Float_t GetCharge (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetCharge (n);}
  Int_t   GetADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks) const {return fDP[ich][ip]->GetADC (n, ks);}
  UInt_t   GetNChambers() const {return fNChambers;};
  UInt_t   GetNPlanes (const UInt_t i) const {return fNPlanes[i];}
  UShort_t GetNSamples (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNSamples();}
  UShort_t GetNStrips (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNStrips();}
  UShort_t Threshold (UInt_t ich, UInt_t ip, Int_t thr) {return fDP[ich][ip]->Threshold (thr);}
  UShort_t GetNOverThr (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNOverThr();}
  Short_t  GetIdxOverThr (UInt_t ich, UInt_t ip, Int_t n) const
  { return fDP[ich][ip]->GetIdxOverThr(n); }

  const std::vector<Short_t>& GetStripClusters(UInt_t ich, UInt_t ip, UInt_t n) const
  { return fDP[ich][ip]->GetStripClusters(n); }

  TSolSimEvent* GetEvent() const { return fEvent; }

  Bool_t IsMapSector() const { return fDoMapSector; }
  void SetMapSector( Bool_t b = true ) { fDoMapSector = b; }

 private:

  void IonModel (const TVector3& xi,
		 const TVector3& xo,
		 const Double_t elost );
  
  TSolGEMVStrip ** AvaModel (const Int_t ic,
			     const TSolSpec& spect,
			     const TVector3& xi,
			     const TVector3& xo,
			     const Double_t time_off);
  
  // Gas parameters
  Double_t fGasWion;               // eV
  Double_t fGasDiffusion;          // mm2/s
  Double_t fGasDriftVelocity;      // mm/s
  Double_t fAvalancheFiducialBand; // number of sigma defining the band around the avalanche in readout plane
  Int_t    fAvalancheChargeStatistics;  // 0 Furry, 1 Gaussian
  Double_t fGainMean;
  Double_t fGain0;

  // Electronics parameters
  Double_t fTriggerOffset;       // trigger offset (ns), incl latency & readout offset
  Double_t fTriggerJitter;       // trigger sigma jitter (ns)
  Int_t    fEleSamplingPoints;
  Double_t fEleSamplingPeriod;   // ns
  Double_t fPulseNoiseSigma; // sigma of the amplitude noise distribution on each sample
  Double_t fADCoffset;       // ADC offset
  Double_t fADCgain;         // ADC gain
  Int_t    fADCbits;         // ADC resolutions in bits
  Double_t fGateWidth;       // to be changed , ns - pulse shape width at ~1/10 max

  // Pulse shaping parameters
  Double_t fPulseShapeTau0;   // [ns] GEM model; = 50. in SiD model
  Double_t fPulseShapeTau1;   // [ns] GEM model only; if negative assume SiD model

  // Geometry
  Double_t fRoutZ;            // z-distance hit entrance to readout plane [mm]
  
  // Sector mapping
  Bool_t   fDoMapSector;
  Int_t    fSignalSector;

  TSolDigitizedPlane*** fDP; // 2D array of plane pointers indexed by chamber, plane #
  TSolGEMVStrip** fdh;// array of U & V GEM strips
  
  UInt_t fNChambers;  // # chambers
  UInt_t* fNPlanes;   // # planes in each chamber
  TRandom3 fTrnd;     // time randomizer
  UInt_t   fRNIon;    // number of ions
  struct IonPar_t {
    Double_t X;       // position of the point on the projection
    Double_t Y;
    Double_t Charge;  // Charge deposited by this ion
    Double_t SNorm;   // 3 x radius of ion diffusion area at readout
    Double_t R2;      // = SNorm^2 : radius of numerical integration area
    Double_t ggnorm;  // = Charge/R2/pi : charge per unit area
  };
  std::vector<IonPar_t> fRIon;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;
  
  std::vector<Double_t> fSumA;
  std::vector<Short_t>  fDADC;
  
  // Tree

  TFile* fOFile;          // Output ROOT file
  TTree* fOTree;          // Output tree
  TSolSimEvent* fEvent;   // Output event structure, written to tree

  Bool_t fFilledStrips;   // True if no data changed since last SetTreeStrips

  void MakePrefix() { THaAnalysisObject::MakePrefix(0); }
  void DeleteObjects();

  ClassDef (TSolSimGEMDigitization, 0)
};

#endif

