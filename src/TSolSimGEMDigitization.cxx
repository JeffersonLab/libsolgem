#include "TSolSimGEMDigitization.h"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "TCanvas.h"
#include "TF2.h"
#include "TH2F.h"
#include "TMath.h"

#include "TSolGEMData.h"
#include "TSolGEMVStrip.h"
#include "TSolSpec.h"
#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"
#include "TSolSimAux.h"

using namespace std;

// Auxiliary class

TSolDigitizedPlane::TSolDigitizedPlane (Short_t nstrip,
					Short_t nsample) 
{
  fNOT = 0;
  fNSample = nsample;
  fNStrip = nstrip;

  fType = new Short_t[nstrip];
  fCharge = new Float_t[nstrip];
  fTime = new Float_t[nstrip];
  fTotADC = new Short_t[nstrip];
  
  for (Int_t i=0;i<nstrip;i++) {
    fTotADC[i]=0.;
    fType[i]=0;
    fCharge[i] = 0.;
    fTime[i] = 9999.;
  }
  
  fPStripADC = new TArrayS(fNSample*nstrip);
  
  fPStripADC->Reset();
  
  if (fPStripADC==0) {
    cerr << __FUNCTION__ << " allocation failed" << endl;
  }
};

TSolDigitizedPlane::~TSolDigitizedPlane() 
{
  delete fPStripADC;
  delete[] fType;
  delete[] fCharge;
  delete[] fTime;
  delete[] fTotADC;
};

void 
TSolDigitizedPlane::Cumulate (TSolGEMVStrip *vv, Int_t type) const
{
  
  Int_t j,k;
  Int_t idx;

  Short_t ooo,nnn;

  if (vv) {
    for (j=0;j<vv->GetSize();j++) {
      idx = vv->GetIdx(j);
      fType[idx] |= (Short_t) type;
      fTime[idx] = (fTime[idx] < vv->GetTime()) ? fTime[idx] : vv->GetTime();
      fCharge[idx] += vv->GetCharge(j);
      for (k=0;k<fNSample;k++) {
	ooo=fPStripADC->At(idx*fNSample+k);
	nnn=vv->GetADC(j,k);
	fPStripADC->AddAt(ooo+nnn, idx*fNSample+k);
	fTotADC[idx] += nnn;
      }
    }
  }
};



TSolSimGEMDigitization::TSolSimGEMDigitization(const TSolSpec& spect)
{
  Initialize (spect);
}


TSolSimGEMDigitization::~TSolSimGEMDigitization()
{
  UInt_t nc = fSpect->GetNChambers();
  for (UInt_t ic = 0; ic < nc; ++ic)
    {
      UInt_t np = fSpect->GetChamber (ic).GetNPlanes();
      for (UInt_t ip = 0; ip < np; ++ip)
	delete fDP[ic][ip];
      delete[] fDP[ic];
    }
  delete[] fDP;
}

void 
TSolSimGEMDigitization::Initialize(const TSolSpec& spect)
{
  fSpect = &spect;

  UInt_t nc = fSpect->GetNChambers();
  fDP = new TSolDigitizedPlane**[nc];
  for (UInt_t ic = 0; ic < nc; ++ic)
    {
      UInt_t np = fSpect->GetChamber(ic).GetNPlanes();
      fDP[ic] = new TSolDigitizedPlane*[np];
      for (UInt_t ip = 0; ip < np; ++ip)
	fDP[ic][ip] = new TSolDigitizedPlane (fSpect->GetChamber(ic).GetPlane(ip).GetNStrips());
    }

  SetGasParams();
  SetEleParams();
  SetPulseShaping();
}

void 
TSolSimGEMDigitization::SetGasParams (Double_t wion, // eV 
				      Double_t diff,  // cm2/s
				      Double_t vDrift, // cm/s
				      Double_t avFidBand, // n sigma
				      Int_t    avChargeStats,  // 0 Furry, 1 Gaussian
				      Double_t gainMean,
				      Double_t gain0
				      )
{
  fGasWion = wion;
  fGasDiffusion = diff* 100;
  fGasDriftVelocity = vDrift*10.;
  fAvalancheFiducialBand = avFidBand; // number of sigma defining the band around the avalanche in readout plane
  fAvalancheChargeStatistics = avChargeStats;
  fGainMean = gainMean;
  fGain0 = gain0;
}

void 
TSolSimGEMDigitization::SetEleParams (Double_t off,  // trigger offset [ns] 
				      Double_t jit,   // trigger sigma jitter [ns]
				      Double_t sampleTime, // 25 ns
				      Int_t samplePoints,
				      Double_t pulseNoiseSigma,
				      Int_t thr,       // electronics sparse readout threshold (adc unit)
				      Double_t adc_off, 
				      Double_t adc_gain, 
				      Int_t adc_bits,
				      Double_t gateWidth
				      )
{
  fTriggerOffset = off;
  fTriggerJitter = jit;
  fEleSamplingPoints = samplePoints;
  fEleSamplingPeriod = sampleTime;
  fPulseNoiseSigma = pulseNoiseSigma;
  fADCoffset = adc_off;
  fADCgain = adc_gain;
  fADCbits = adc_bits;
  fGateWidth = gateWidth;
}

void 
TSolSimGEMDigitization::SetPulseShaping (Double_t tau0, // [ns] GEM model; = 50. in SiD model 
					 Double_t tau1 // [ns] GEM model only; if negative assume SiD model
					 ) 
{
  fPulseShapeTau0 = tau0;
  fPulseShapeTau1 = tau1;
}

Int_t 
TSolSimGEMDigitization::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  const DBRequest request[] = 
    {
      {"gasionwidth",               &fGasWion,                   kDouble, 0, 1, 0},
      {"gasdiffusion",              &fGasDiffusion,              kDouble, 0, 1, 0},
      {"gasdriftvelocity",          &fGasDriftVelocity,          kDouble, 0, 1, 0},
      {"avalanchefiducialband",     &fAvalancheFiducialBand,     kDouble, 0, 1, 0},
      {"avalanchechargestatistics", &fAvalancheChargeStatistics, kInt,    0, 1, 0},
      {"gainmean",                  &fGainMean,                  kDouble, 0, 1, 0},
      {"gain0",                     &fGain0,                     kDouble, 0, 1, 0},
      {"triggeroffset",             &fTriggerOffset,             kDouble, 0, 1, 0},
      {"triggerjitter",             &fTriggerJitter,             kDouble, 0, 1, 0},
      {"elesamplingpoints",         &fEleSamplingPoints,         kInt,    0, 1, 0},
      {"elesamplingperiod",         &fEleSamplingPeriod,         kDouble, 0, 1, 0},
      {"pulsenoisesigma",           &fPulseNoiseSigma,           kDouble, 0, 1, 0},
      {"adcoffset",                 &fADCoffset,                 kDouble, 0, 1, 0},
      {"adcgain",                   &fADCgain,                   kDouble, 0, 1, 0},
      {"adcbits",                   &fADCbits,                   kInt,    0, 1, 0},
      {"gatewidth",                 &fGateWidth,                 kDouble, 0, 1, 0},
      {"pulseshapetau0",            &fPulseShapeTau0,            kDouble, 0, 1, 0},
      {"pulseshapetau1",            &fPulseShapeTau1,            kDouble, 0, 1, 0},
      {0}
    };

  Int_t err = LoadDB (file, date, request, fPrefix);
  fclose(file);
  if (err)
    return kInitError;

  return kOK;
}

void 
TSolSimGEMDigitization::Digitize (const TSolGEMData& gdata) // digitize event 
{
  UInt_t nh = gdata.GetNHit();

  for (UInt_t ih = 0; ih < nh; ++ih)
    {
      Int_t igem = gdata.GetHitChamber (ih);
      if (igem >= (Int_t) fSpect->GetNChambers())
	continue;
      
      Short_t itype = (1 << gdata.GetParticleType(ih)); // signal = 1, bck = 2, 4, 8 ...
	
      TVector3 vv1 = gdata.GetHitEntrance (ih);
      TVector3 vv2 = gdata.GetHitExit (ih);
      TVector3 vv3 = gdata.GetHitReadout (ih);

      // These vectors are in the lab frame, we need them in the chamber frame
      // Also convert to mm

      TVector3 offset = fSpect->GetChamber(igem).GetOrigin() * 1000.0;
      Double_t angle = fSpect->GetChamber(igem).GetAngle();
      vv1 -= offset;
      vv2 -= offset;
      vv3 -= offset;
      vv1.RotateZ (angle);
      vv2.RotateZ (angle);
      vv3.RotateZ (angle);
	
      IonModel (vv1, vv2, gdata.GetHitEnergy(ih), vv3);
      if (fRNIon > 0) 
	{
	  Double_t time_zero = 
	    (itype == 1) ? 0.
	    : fTrnd.Uniform (fGateWidth + 75.) - fGateWidth; // randomization of the bck ( assume 3 useful samples at 25 ns)
	  TSolGEMVStrip **dh = AvaModel (igem, vv1, vv2, time_zero);
	  for (UInt_t j = 0; j < 2; j++) 
	    {
	      fDP[igem][j]->Cumulate (dh[j], itype);
	      delete dh[j];
	    }
	  delete[] dh;
	} 
    }
}
 

//-------------------------------------------------------
//.......................................................
// ionization Model
//

void
TSolSimGEMDigitization::IonModel(const TVector3& xi, 
				 const TVector3& xo, 
				 const Double_t elost, 
				 const TVector3& xrout)   // used only to calculate distance drift-readout (to be removed in future version)
{

  Double_t LL;
  Double_t deltaE=elost; // eV  MC

  TVector3 vseg = xo-xi;

  // DEBUG  TRandom3 rnd(0);
  TRandom3 rnd;

  // ---- extract primary ions from Poisson

  fRNIon = rnd.Poisson(deltaE/fGasWion);

  if (fRNIon <=0)
    return;

  if (fRNIon > 200) {
    cerr << __FUNCTION__ << ": WARNING: too many primary ions " << fRNIon << " limit to 200" << endl;
    fRNIon = 200;
  }

  Float_t lion;
  Float_t ttime;
  Float_t gnorm;

  fRSNorm  = new Double_t[fRNIon];
  fRCharge = new Double_t[fRNIon];

  fRX = new Double_t[fRNIon]; // position of the point on the projection
  fRY = new Double_t[fRNIon]; 

  fRSMax = 0.;
  fRTotalCharge = 0;

  fRTime0 = 999999.; // minimum time of drift

  for (UInt_t i=0;i<fRNIon;i++) { // first loop used to evaluate quantities

    lion  = rnd.Uniform(0.,1.); // position of the hit along the track segment (fraction)

    fRX[i]=vseg.X()*lion+xi.X(); 
    fRY[i]=vseg.Y()*lion+xi.Y();

    LL = TMath::Abs(xrout.Z() - (vseg.Z()*lion+xi.Z()));

    ttime = LL/fGasDriftVelocity; // traveling time from the drift gap to the readout

    fRTime0 = (ttime<fRTime0) ? ttime : fRTime0; // minimum traveling time 

    gnorm = fGainMean/sqrt(fGain0); // overall gain TBC

    fRSNorm[i] = sqrt(2.*fGasDiffusion*ttime); // spatial sigma on readout

    switch (fAvalancheChargeStatistics) {
    case 1:
      fRCharge[i]= rnd.Gaus(fGainMean, gnorm); // Gaussian distribution of the charge
      break;
    default: 
      fRCharge[i]= rnd.Exp(fGainMean); // Furry distribution
      break;
    }

    fRTotalCharge += fRCharge[i];

    fRSMax = (fRSNorm[i]>fRSMax) ? fRSNorm[i] : fRSMax;

  }

}

//.......................................................
// avalanche model
//

TSolGEMVStrip **
TSolSimGEMDigitization::AvaModel(const Int_t ic,
				 const TVector3& xi, 
				 const TVector3& xo,
				 const Double_t time_off)
{
  Double_t nsigma = fAvalancheFiducialBand; // coverage factor

  // DEBUG TRandom3 rnd(0);
  TRandom3 rnd;

  Double_t x0,y0,x1,y1; // active window lower and upper corners

  // check if track is in the active area of the sector

  if (xi.X()<xo.X()) {
    x0 = xi.X()-nsigma*fRSMax;
    x1 = xo.X()+nsigma*fRSMax;
  } else {
    x1 = xi.X()+nsigma*fRSMax;
    x0 = xo.X()-nsigma*fRSMax;
  }

  if (xi.Y()< xo.Y()) {
    y0 = xi.Y()-nsigma*fRSMax;
    y1 = xo.Y()+nsigma*fRSMax;
  } else {
    y1 = xi.Y()+nsigma*fRSMax;
    y0 = xo.Y()-nsigma*fRSMax;
  }

  // --- loop on sectors 

  Double_t glx = fSpect->GetChamber(ic).GetLowerEdgeX() * 1000.0;
  Double_t gly = fSpect->GetChamber(ic).GetLowerEdgeY() * 1000.0;
  Double_t gux = fSpect->GetChamber(ic).GetUpperEdgeX() * 1000.0;
  Double_t guy = fSpect->GetChamber(ic).GetUpperEdgeY() * 1000.0;

  if (x1<glx || x0>gux ||
      y1<gly || y0>guy) { // out of active area of the sector
    delete[] fRSNorm;
    delete[] fRCharge;
    delete[] fRX; 
    delete[] fRY; 
    return 0;
  }

  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes

  UInt_t np = fSpect->GetChamber(ic).GetNPlanes();
  TSolGEMVStrip **virs;
  virs = new TSolGEMVStrip *[np];
  for (UInt_t ipl = 0; ipl < np; ++ipl)
    {
      // Compute strips affected by the avalanche

      TSolGEMPlane* pl = &(fSpect->GetChamber(ic).GetPlane(ipl)); 
      Double_t z = (pl->GetOrigin())[2] * 1000.0;

      // Positions in strip frame
      TVector3 v0 = pl->LabToStrip (TVector3 (x0, y0, z) * 1e-3) * 1000.0;
      TVector3 v1 = pl->LabToStrip (TVector3 (x1, y1, z) * 1e-3) * 1000.0;
      Int_t iL = pl->GetStrip (v0[0] * 1e-3, v0[1] * 1e-3);
      Int_t iU = pl->GetStrip (v1[0] * 1e-3, v1[1] * 1e-3);
      if (iL > iU)
	{
	  Int_t t = iL;
	  iL = iU;
	  iU = t;
	}

      //
      // Bounds of rectangular avalanche region, in strip frame
      //

      // Limits in x are low edge of first strip to high edge of last
      Double_t xl = pl->GetStripLowerEdge (iL) * 1000.0;
      Double_t xr = pl->GetStripUpperEdge (iU) * 1000.0;

      // Limits in y are y limits of track plus some reasonable margin
      // We do this in units of strip pitch for convenience (even though
      // this is the direction orthogonal to the pitch direction)

      Double_t pitch = pl->GetSPitch() * 1000.0;
      Double_t yb = v0[1] - 10 * pitch;
      yb = pitch * TMath::Floor (yb / pitch);
      Double_t yt = v1[1] + 10 * pitch;
      yt = pitch * TMath::Floor (yt / pitch);
      if (yb > yt)
	{
	  Double_t t = yb;
	  yb = yt;
	  yt = t;
	}

      // # of coarse histogram bins: 1 per pitch each direction

      Int_t nx = iU - iL + 1;
      Int_t ny = (Int_t) ((yt - yb) / pitch + 0.5);

      // # of fine histogram bins
      Int_t nnx = (nx < 20) ? nx * 4 : nx; // TF2 bins cannot be <4 !!! and high precision with higher bins
      Int_t nny = (ny < 20) ? ny * 4 : ny;

      // define function, gaussian and sum of gaussian

      TF2 *sgaus = new TF2[fRNIon];
      TH2F *fsum = 0;
  
      for (UInt_t i = 0; i < fRNIon; i++) 
	{ 
	  // define gaussian functions (for spatial distribution)
	  sgaus[i] = TF2 (Form ("sgaus%d", i),
			  TSolSimAux::SimpleCircle,
			  xl, xr, yb, yt, 4);
	  sgaus[i].SetParNames ("Const", "X_{0}", "Y_{0}", "#sigma");
	  Double_t ggnorm = fRCharge[i] / 3.14 / 9. / fRSNorm[i] / fRSNorm[i]; // normalized to charge
	  TVector3 fr = pl->LabToStrip (TVector3 (fRX[i], fRY[i], z) * 1e-3) * 1000.0;
	  sgaus[i].SetParameters (ggnorm, fr.X(), fr.Y(), 3. * fRSNorm[i]);
	  sgaus[i].SetNpx (nnx); 
	  sgaus[i].SetNpy (nny);
	  if (i == 0)
	    {
	      fsum = (TH2F*) sgaus[i].CreateHistogram(); // on first loop, create histo
	    }
	  else
	    {
	      fsum->Add (&sgaus[i], 1.0); 
	    }
	}


      Double_t scale = ((Double_t) nnx)/((Double_t) nx) * ((Double_t) nny)/((Double_t) ny);

      fsum->Rebin2D (nnx/nx, nny/ny); // rebin (fix problem on TF2 that cannot have nbin<4) ; warning rebin does not rescale width

      Double_t *us = new Double_t[nx];
      for (Int_t j = 0; j < nx; j++) 
	us[j] = 0;

      for (Int_t i = iL; i <= iU; i++)
	{
	  Int_t kk = i - iL;
	  us[kk] += fsum->Integral (kk+1, kk+1, 1, ny, "width") / scale;
	}

      Float_t t0 = time_off + fRTime0
	+ rnd.Gaus(fTriggerOffset, fTriggerJitter); // ... NEED Time of flight from GEANT4

      virs[ipl] = new TSolGEMVStrip(nx,fEleSamplingPoints);

      virs[ipl]->SetTime(t0);
      virs[ipl]->SetHitCharge(fRTotalCharge);

      Double_t pulse=0.;
      Double_t noisy_pulse;

      Short_t *dadc = new Short_t[fEleSamplingPoints];
      Int_t ai=0;
      Int_t posflag=0;

      Double_t cccsssx = 0.0;

      for (Int_t i = iL; i <= iU; i++) 
	{
	  posflag = 0;
	  for (Int_t b = 0; b < fEleSamplingPoints; b++) 
	    { // sampling
	      pulse = TSolSimAux::PulseShape (t0 + fEleSamplingPeriod * b, 
					      us[i-iL], 
					      fPulseShapeTau0, 
					      fPulseShapeTau1);
	      noisy_pulse = rnd.Gaus (pulse, fPulseNoiseSigma); 
	      dadc[b] = TSolSimAux::ADCConvert (noisy_pulse, 
						fADCoffset, 
						fADCgain, 
						fADCbits);
	      posflag += (Int_t) dadc[b];
	    }
	  if (posflag > 0) 
	    { // store only strip with signal -- do not work yet
	      for (Int_t b = 0; b < fEleSamplingPoints; b++)
		virs[ipl]->AddSampleAt (dadc[b], b, ai);
	      virs[ipl]->AddStripAt (i, ai);
	      virs[ipl]->AddChargeAt (us[i-iL], ai);
	      cccsssx += us[i-iL];
	      ai++;
	    }
	}
      virs[ipl]->SetSize(ai);

      //

      delete fsum;
//       if (ipl == 0)
// 	{
// 	  fsum->SetName ("fsum");
// 	  TCanvas* c1;
// 	  fsum->DrawCopy();
// 	}
      
      delete[] dadc;
      delete[] us;
      
      delete[] sgaus;
    }

  delete[] fRSNorm;
  delete[] fRCharge;
  delete[] fRX; 
  delete[] fRY;
      
  return virs;

}

void 
TSolSimGEMDigitization::Print() const
{
  cout << "GEM digitization:" << endl;
  cout << "  Gas parameters:" << endl;
  cout << "    Gas ion width: " << fGasWion << endl;
  cout << "    Gas diffusion: " << fGasDiffusion << endl;
  cout << "    Gas drift velocity: " << fGasDriftVelocity << endl;
  cout << "    Avalanche fiducial band: " << fAvalancheFiducialBand << endl;
  cout << "    Avalanche charge statistics: " << fAvalancheChargeStatistics << endl;
  cout << "    Gain mean: " << fGainMean << endl;
  cout << "    Gain 0: " << fGain0 << endl;
  
  cout << "  Electronics parameters:" << endl;
  cout << "    Trigger offset: " << fTriggerOffset << endl;
  cout << "    Trigger jitter: " << fTriggerJitter << endl;
  cout << "    Sampling Period: " << fEleSamplingPeriod << endl;
  cout << "    Sampling Points: " << fEleSamplingPoints   << endl;
  cout << "    Pulse Noise width: " << fPulseNoiseSigma << endl;
  cout << "    ADC offset: " << fADCoffset << endl;
  cout << "    ADC gain: " << fADCgain << endl;
  cout << "    ADC bits: " << fADCbits << endl;
  cout << "    Gate width: " << fGateWidth << endl;
  
  cout << "  Pulse shaping parameters:" << endl;
  cout << "    Pulse shape tau0: " << fPulseShapeTau0 << endl;
  cout << "    Pulse shape tau1: " << fPulseShapeTau1 << endl;
}


void
TSolSimGEMDigitization::PrintResults() const
{
  cout << " Chb  Pln  Strip  Typ    ADC    Charge      Time\n";
  UInt_t nc = fSpect->GetNChambers();
  for (UInt_t ic = 0; ic < nc; ++ic)
    {
      UInt_t np = fSpect->GetChamber (ic).GetNPlanes();
      for (UInt_t ip = 0; ip < np; ++ip)
	for (UInt_t ist = 0; ist < (UInt_t) fSpect->GetChamber (ic).GetPlane (ip).GetNStrips(); ++ist)
	  {
	    if (fDP[ic][ip]->GetCharge (ist) > 0)
	      cout << setw(4) << ic
		   << " " << setw(4) << ip
		   << " " << setw(6) << ist
		   << " " << setw(4) << fDP[ic][ip]->GetType (ist)
		   << " " << setw(6) << fDP[ic][ip]->GetTotADC (ist)
		   << fixed << setprecision(1)
		   << " " << setw(9) << fDP[ic][ip]->GetCharge (ist)
		   << fixed << setprecision(3)
		   << " " << setw(9) << fDP[ic][ip]->GetTime (ist)
		   << endl;
	  }
    }
}
