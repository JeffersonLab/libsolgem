// Class handling digitization of GEMs

#ifndef __TSolSimGEMDigitization__
#define __TSolSimGEMDigitization__

#include "Rtypes.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayS.h"

#include "THaAnalysisObject.h"

class TFile;
class TTree;

class TSolGEMData;
class TSolGEMVStrip;
class TSolSpec;

// First an auxiliary class

// The whole strip plane; used to cumulate virtual strips charges
// and generate real strips

class TSolDigitizedPlane 
{
 private:
  // ADC sampled value of strip array of each axis

  TArrayS *fPStripADC;
  Short_t *fType;
  Short_t *fTotADC;

  Float_t *fCharge;
  Float_t *fTime;

  Short_t fNSamples;
  Short_t fNStrips;

  UInt_t fNOT;   // # strips over threshold
  UInt_t* fOverThr;  // # list of strips over threshold

 public:
  // init and reset physics strips arrays
  TSolDigitizedPlane (Short_t nstrip,
		      Short_t nsample = 10);
  ~TSolDigitizedPlane();
  void Init();

  // cumulate hits (strips signals)
  void Cumulate (TSolGEMVStrip *vv, Int_t type) const;

  Short_t GetType (Int_t n) const {return fType[n];}
  Short_t GetTotADC (Int_t n) const {return fTotADC[n];}
  Float_t GetTime (Int_t n) const {return fTime[n];};
  Float_t GetCharge (Int_t n) const {return fCharge[n];};
  Short_t GetADC (Int_t n, Int_t ks) const {return fPStripADC->At(n*fNSamples+ks);};
  Short_t GetNSamples() const {return fNSamples;}
  Short_t GetNStrips() const {return fNStrips;}

  UInt_t Threshold (Short_t thr);

  UInt_t GetNOverThr() const {return fNOT;};
  UInt_t GetIdxOverThr (Int_t n) const {return fOverThr[n];}; 
};



class TSolSimGEMDigitization: public THaAnalysisObject
{
 public:
  TSolSimGEMDigitization(const TSolSpec& spect);
  virtual ~TSolSimGEMDigitization();

  void Initialize(const TSolSpec& spect);
  void SetGasParams (Double_t wion = 26., // eV 
		     Double_t diff = 250.,  // cm2/s
		     Double_t vDrift = 5.5*1e6, // cm/s
		     Double_t avFidBand = 10.0, // n sigma
		     Int_t    avChargeStats = 0, // 0 Furry, 1 Gaussian
		     Double_t gainMean = 8000.0,
		     Double_t gain0 = 20.0
		     );
  void SetEleParams (Double_t off=10.,  // trigger offset [ns] 
		     Double_t jit=5.,   // trigger sigma jitter [ns]
		     Double_t sampleTime = 25., // 25 ns
		     Int_t samplePoints = 10,
		     Double_t pulseNoiseSigma = 100.0, // sigma of the gaussian noise distribution on each sample
		     Int_t thr=1,       // electronics sparse readout threshold (adc unit)
		     Double_t adc_off=0., 
		     Double_t adc_gain=1., 
		     Int_t adc_bits=12,
		     Double_t gateWidth = 100.
		     );
  void SetPulseShaping (Double_t tau0=38.3, // [ns] GEM model; = 50. in SiD model 
			Double_t tau1=129. // [ns] GEM model only; if negative assume SiD model
			);
  Int_t ReadDatabase (const TDatime& date);

  void Digitize (const TSolGEMData& gdata, const TSolSpec& spect); // digitize event  
  const TSolDigitizedPlane& GetDigitizedPlane (UInt_t ich, UInt_t ip) const {return *(fDP[ich][ip]);};
  void Print() const;
  void PrintCharges() const;
  void PrintSamples() const;

  // Tree methodss
  void InitTree (const TSolSpec& spect, const TString& ofile);
  void FillTreeHit (const UInt_t ih, 
		    const UInt_t igem, 
		    TSolGEMVStrip** dh,
		    const TSolGEMData& tsgd);
  void FillTreeHit (const TSolGEMData& gdata);
  void FillTreeEvent (const TSolGEMData& gdata);
  void WriteTree () const;
  void CloseTree () const;

  // Access to results
  Short_t GetType (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetType (n);};
  Short_t GetTotADC (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetTotADC (n);};
  Float_t GetTime (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetTime (n);};
  Float_t GetCharge (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetCharge (n);};
  Short_t GetADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks) const {return fDP[ich][ip]->GetADC (n, ks);};
  UInt_t GetNChambers() const {return fNChambers;};
  UInt_t GetNPlanes (const UInt_t i) const {return fNPlanes[i];};
  Short_t GetNSamples (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNSamples();}
  Short_t GetNStrips (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNStrips();}
  UInt_t Threshold (UInt_t ich, UInt_t ip, Short_t thr) {return fDP[ich][ip]->Threshold (thr);};
  UInt_t GetNOverThr (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNOverThr();};
  UInt_t GetIdxOverThr (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetIdxOverThr (n);}; 


 private:

  void MakePrefix() {THaAnalysisObject::MakePrefix (NULL);}
  void IonModel (const TVector3& xi, 
		 const TVector3& xo, 
		 const Double_t elost, 
		 const TVector3& xrout);   // used only to calculate distance drift-readout (to be removed in future version)
  
  TSolGEMVStrip ** AvaModel (const Int_t ic,
			     const TSolSpec& spect,
			     const TVector3& xi, 
			     const TVector3& xo,
			     const Double_t time_off);
  

  // Gas parameters
  Double_t fGasWion;
  Double_t fGasDiffusion;
  Double_t fGasDriftVelocity;
  Double_t fAvalancheFiducialBand; // number of sigma defining the band around the avalanche in readout plane
  Int_t    fAvalancheChargeStatistics;  // 0 Furry, 1 Gaussian
  Double_t fGainMean;
  Double_t fGain0;

  // Electronics parameters
  Double_t fTriggerOffset;    // trigger offset (ns)
  Double_t fTriggerJitter;       // trigger sigma jitter
  Int_t    fEleSamplingPoints;   
  Double_t fEleSamplingPeriod;    // ns    
  Double_t fPulseNoiseSigma;  // sigma of the gaussian noise distribution on each sample
  Double_t fADCoffset;       // ADC offset
  Double_t fADCgain;         // ADC gain
  Int_t    fADCbits;         // ADC resolutions in bits
  Double_t fGateWidth;       // to be changed , ns - pulse shape width at ~1/10 max
   
  // Pulse shaping parameters
  Double_t fPulseShapeTau0;   // [ns] GEM model; = 50. in SiD model 
  Double_t fPulseShapeTau1;   // [ns] GEM model only; if negative assume SiD model

  //
  TSolDigitizedPlane*** fDP; // 2D array of plane pointers indexed by chamber, plane #

  UInt_t fNChambers;  // # chambers
  UInt_t* fNPlanes;   // # planes in each chamber
  TRandom3 fTrnd;   // time randomizer
  UInt_t   fRNIon;  // number of ions
  Double_t *fRX;
  Double_t *fRY;
  Double_t *fRSNorm;
  Double_t *fRCharge;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;

  // Tree

  TString fOFileName;  // Name of output ROOT file
  TFile* fOFile;       // Output ROOT file
  TTree* fOTree;

  Int_t fRunID;  // Run number
  Int_t fEvtID;  // Event number

  // cluster (hit) variables  
  Int_t fNTreeHits;   // total number of hits for tree
  Int_t fNSignal;    // number of hits of signal 

  static const Int_t gMAX_CHANNELS = 10000;
  
  Short_t fClsChamber[gMAX_CHANNELS]; // chamber/module index from 0 to 5 (or total number of chambers)
  Float_t fClsCharge[gMAX_CHANNELS]; // charge of the hits avalanche
  Int_t fClsRefEntry[gMAX_CHANNELS]; // reference entry in mc 
  Int_t fClsRefFile[gMAX_CHANNELS]; // reference file in mc (combined with the above can be used to lok at the Montecarlo raw data)
  Float_t fClsTime[gMAX_CHANNELS]; // time of arrival of the signal into the electronics (time of flight not included yet)
  Float_t fClsMx[gMAX_CHANNELS]; // momentum particle generating the hit (x component)
  Float_t fClsMy[gMAX_CHANNELS]; // momentum particle generating the hit (y component)
  Float_t fClsMz[gMAX_CHANNELS]; // momentum particle generating the hit (z component)
  Int_t fClsPID[gMAX_CHANNELS]; // particle ID (PDG code) generating the hit

  Int_t fClsSize[2][gMAX_CHANNELS];   // number of strips in cluster/hit on each axis
  Int_t fClsFirstStrip[2][gMAX_CHANNELS]; // address of first strip of the cluster on each axis
  
  Int_t fDNCh; // number of strips (channels) firing in event
  Short_t fDGEM[gMAX_CHANNELS]; // chamber[7,15] | module[3-6]] | axis[0-2] 
  Short_t fDPlane[gMAX_CHANNELS]; // axis or more generally a readout plane
  Short_t fDStrip[gMAX_CHANNELS]; // strip address 
  Short_t fDSADC[10][gMAX_CHANNELS]; // adc values of the strip
  Short_t fType[gMAX_CHANNELS]; // type of strip hit (0=signal or signal+bck, 1=bck)
  Float_t fCharge[gMAX_CHANNELS]; // total charge in strip
  Float_t fTime1[gMAX_CHANNELS]; // time of first sample since event started in target (TBC)

    public:
	ClassDef (TSolSimGEMDigitization, 1)
};

#endif

