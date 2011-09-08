// Class handling digitization of GEMs

#ifndef __TSolSimGEMDigitization__
#define __TSolSimGEMDigitization__

#include "Rtypes.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayS.h"

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

  Short_t fNOT;
  Short_t *fOverThr;

  TArrayS *fPStripADC;
  Short_t *fType;
  Short_t *fTotADC;

  Float_t *fCharge;
  Float_t *fTime;

  Short_t fNsample;
  Short_t fNstrip;

 public:
  // init and reset physics strips arrays
  TSolDigitizedPlane (Short_t nstrip,
		      Short_t nsample = 10);
  ~TSolDigitizedPlane();

  // cumulate hits (strips signals)
  void Cumulate (TSolGEMVStrip *vv, Int_t type) const;

  Float_t GetTime (Int_t n) const {return fTime[n];};
  Float_t GetCharge (Int_t n) const {return fCharge[n];};
  Short_t GetIdxOverThr (Int_t n) const {return fOverThr[n];};
  Short_t GetTypeOverThr (Int_t n) const {return fType[fOverThr[n]];};
  Short_t GetADCOverThr (Int_t ks, Int_t n) const {return fPStripADC->At(fOverThr[n]*fNsample+ks);};
  Float_t GetChargeOverThr (Int_t n) const {return fCharge[fOverThr[n]];};
  Float_t GetTimeOverThr (Int_t n) const {return fTime[fOverThr[n]];};
};



class TSolSimGEMDigitization
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

  void Digitize (const TSolGEMData& gdata,
		 const TSolSpec& spect); // digitize event 
  const TSolDigitizedPlane& GetDigitizedPlane (UInt_t ich, UInt_t ip) const {return *(fDP[ich][ip]);};

 private:

  Int_t ionModel (Int_t ic,
		  TVector3 *xi, 
		  TVector3 *xo, 
		  Double_t elost, 
		  TVector3 *xrout);   // used only to calculate distance drift-readout (to be removed in future version)
  
  TSolGEMVStrip ** avaModel (Int_t ic,
			     TVector3 *xi, 
			     TVector3 *xo,
			     Double_t time_off);

  

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
  Int_t    fEleThr;       // threshold for sparse readout
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
  const TSolSpec* fSpect; // the spectrometer
  TSolDigitizedPlane*** fDP; // 2D array of plane pointers indexed by chamber, plane #

  TRandom3 fTrnd;   // time randomizer
  UInt_t   fRNIon;  // number of ions
  Double_t *fRX;
  Double_t *fRY;
  Double_t *fRSNorm;
  Double_t *fRCharge;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;
};

#endif

