#include "TSolSimGEMDigitization.h"

#include <cmath>
#include <iostream>

#include "TF2.h"
#include "TH2F.h"

#include "TSolGEMData.h"
#include "TSolGEMVStrip.h"
#include "TSolSpec.h"
#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"
#include "TSolSimAux.h"

using namespace std;

TSolSimGEMDigitization::TSolSimGEMDigitization()
{
  Initialize ();
}

void 
TSolSimGEMDigitization::Initialize()
{
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
  fEleThr = thr;
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

void 
TSolSimGEMDigitization::Digitize (const TSolGEMData& gdata,
				  const TSolSpec& spect)
{
  fSpect = &spect;

  UInt_t nh = gdata.GetNHit();
  TSolGEMVStrip ***dh = new TSolGEMVStrip **[nh];

  for (UInt_t ih = 0; ih < nh; ++ih)
    {
      Int_t igem = gdata.GetHitChamber (ih);
      if (igem >= (Int_t) fSpect->GetNChambers())
	continue;
      
      Short_t itype = (1 << gdata.GetParticleType(ih)); // signal = 1, bck = 2, 4, 8 ...
	
      TVector3 *vv1 = gdata.GetHitEntrance (ih);
      TVector3 *vv2 = gdata.GetHitExit (ih);
      TVector3 *vv3 = gdata.GetHitReadout (ih);

      Double_t angle = fSpect->GetChamber(igem).GetPlane(0).GetSAngleComp();
      vv1->RotateZ (angle);
      vv2->RotateZ (angle);
      vv3->RotateZ (angle);
	
      if (ionModel (igem, vv1, vv2, gdata.GetHitEnergy(ih), vv3) > 0) 
	{
	  Double_t time_zero = 
	    (itype == 0) ? 0.
	    : fTrnd.Uniform (fGateWidth + 75.) - fGateWidth; // randomization of the bck ( assume 3 useful samples at 25 ns)
	  dh[ih] = avaModel (igem, vv1, vv2, time_zero);
	}
      else
	dh[ih] = 0;
    }
      
  // -- cumulate hits on chamber planes

  // This needs rethinking! ===============================
      
//   for (ic = 0; ic < nchamb; ic++) 
//     { // chamber

//       for (j = 0; j < 2; j++) 
// 	{ // axis
// 	  gPlane[j] = 0;
// 	}

//       for (UInt_t ih = 0; ih < nh; ih++) 
// 	{ 
// 	  // hits
// 	  Int_t igem = gdata.GetHitChamber (ih);
// 	  if (igem >= (Int_t) fSpect->GetNChambers())
// 	    continue;
// 	  if (dh[ih]) 
// 	    {
// 	      if (igem == ic) 
// 		{
// 		  if (gPlane[0] == 0) { // first time of chamber ic
// 		    gPlane[0] = new TSolGEMPlane (fSpect->GetChamber(igem).GetPlane(0).GetNStrips(), fEleSamplingPoints);
// 		    gPlane[1] = new TSolGEMPlane (fSpect->GetChamber(igem).GetPlane(1).GetNStrips(), fEleSamplingPoints);
// 		  }
	      
// 		  for (j = 0; j < 2; j++) 
// 		    {
// 		      gPlane[j]->Cumulate (dh[ih][j], itype);
// 		    }

// 		}

// 	    }

// 	}
		
//     } // end loop on chamber

  for (UInt_t ih = 0; ih < nh; ih++) 
    {
      if (dh[ih]) 
	{
	  for (UInt_t j = 0; j < 2; j++) 
	    delete dh[ih][j];
	  delete dh[ih];
	  dh[ih] = 0;
	}
    }
      
  delete[] dh;
}


//-------------------------------------------------------
//.......................................................
// ionization Model
//

Int_t 
TSolSimGEMDigitization::ionModel(Int_t ic,
				 TVector3 *xi, 
				 TVector3 *xo, 
				 Double_t elost, 
				 TVector3 *xrout)   // used only to calculate distance drift-readout (to be removed in future version)
{

  if ((xi==0) || (xo==0) || (xrout==0)) {
    return 0;
  }

  Double_t LL;
  Double_t deltaE=elost; // eV  MC

  TVector3 vseg = *xo-*xi;

  // DEBUG  TRandom3 rnd(0);
  TRandom3 rnd;

  // ---- extract primary ions from Poisson

  fRNIon = rnd.Poisson(deltaE/fGasWion);

  if (fRNIon <=0) {
    return 0;
  }

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

    fRX[i]=vseg.X()*lion+xi->X(); 
    fRY[i]=vseg.Y()*lion+xi->Y();

    LL = TMath::Abs(xrout->Z() - (vseg.Z()*lion+xi->Z()));

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

  return fRNIon;

}

//.......................................................
// avalanche model
//

TSolGEMVStrip **
TSolSimGEMDigitization::avaModel(Int_t ic,
				 TVector3 *xi, 
				 TVector3 *xo,
				 Double_t time_off)
{
  Double_t nsigma = fAvalancheFiducialBand; // coverage factor

  // DEBUG TRandom3 rnd(0);
  TRandom3 rnd;

  Double_t x0,y0,x1,y1; // active window lower and upper corners

  // check if track is in the active area of the sector

  if (xi->X()<xo->X()) {
    x0 = xi->X()-nsigma*fRSMax;
    x1 = xo->X()+nsigma*fRSMax;
  } else {
    x1 = xi->X()+nsigma*fRSMax;
    x0 = xo->X()-nsigma*fRSMax;
  }

  if (xi->Y()< xo->Y()) {
    y0 = xi->Y()-nsigma*fRSMax;
    y1 = xo->Y()+nsigma*fRSMax;
  } else {
    y1 = xi->Y()+nsigma*fRSMax;
    y0 = xo->Y()-nsigma*fRSMax;
  }

  // --- loop on sectors 

  Double_t glx = fSpect->GetChamber(ic).GetLowerEdgeX();
  Double_t gly = fSpect->GetChamber(ic).GetLowerEdgeY();
  Double_t gux = fSpect->GetChamber(ic).GetUpperEdgeX();
  Double_t guy = fSpect->GetChamber(ic).GetUpperEdgeY();

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


  // Compute strips affected by the avalanche

  Int_t iL;
  Int_t iU;
  Int_t jB;
  Int_t jT;
  Double_t dxl;
  Double_t dxu;
  Double_t dyb;
  Double_t dyt;

  //

  Double_t ppx = fSpect->GetChamber(ic).GetPlane(0).GetPPitch();
  Double_t ppy = fSpect->GetChamber(ic).GetPlane(1).GetPPitch();
  UInt_t npx = fSpect->GetChamber(ic).GetPlane(0).GetNPixels();
  UInt_t npy = fSpect->GetChamber(ic).GetPlane(1).GetNPixels();
  dxl=TMath::Abs(x0-glx) / ppx;
  iL=floor(dxl);

  dxu=TMath::Abs(x1-glx) / ppx;
  iU = (ceil(dxu) >= npx) ? (npx-1) : ceil(dxu);

  dyb=TMath::Abs(y0-gly) / ppy;
  jB=floor(dyb);

  dyt=TMath::Abs(y1-gly) / ppy;
  jT = (ceil(dyt) >= npy) ? (npy-1) : ceil(dyt);

  //
  // new integration domain
  //

  Double_t xl;
  Double_t xr;
  Double_t yb;
  Double_t yt;

  xl=fSpect->GetChamber(ic).GetPlane(0).GetPixelLocation(iL); // Left-Lower corner of the pixel
  xr=fSpect->GetChamber(ic).GetPlane(0).GetPixelLocation(iU)+fSpect->GetChamber(ic).GetPlane(0).GetPPitch();
  yb=fSpect->GetChamber(ic).GetPlane(1).GetPixelLocation(jB);
  yt=fSpect->GetChamber(ic).GetPlane(1).GetPixelLocation(jT)+fSpect->GetChamber(ic).GetPlane(0).GetPPitch();

  Int_t nx=0;
  Int_t ny=0;

  Int_t nnx, nny; // bins for the function, which must be >=4!

  nx = iU-iL+1;
  ny = jT-jB+1;

  nnx = (nx<20) ? nx*4 : nx; // TF2 bins cannot be <4 !!! and high precision with higher bins
  nny = (ny<20) ? ny*4 : ny;

  // define function, gaussian and sum of gaussian

  TF2 *sgaus = new TF2[fRNIon];

  TH2F *fsum = 0;
  
  for (Int_t i=0;i<fRNIon;i++) { // define gaussian functions (for spatial distribution)

    sgaus[i] = TF2(Form("sgaus%d",i),TSolSimAux::SimpleCircle,xl,xr,yb,yt,4);

    sgaus[i].SetParNames("Const","X_{0}","Y_{0}","#sigma");
    Double_t ggnorm = fRCharge[i]/3.14/9./fRSNorm[i]/fRSNorm[i]; // normalized to charge

    sgaus[i].SetParameters(ggnorm,fRX[i],fRY[i],3.*fRSNorm[i]);

    sgaus[i].SetNpx(nnx); 
    sgaus[i].SetNpy(nny);

    if (i==0) {
      fsum=(TH2F*) sgaus[i].CreateHistogram(); // on first loop, create histo
    } else {
      fsum->Add(&sgaus[i],1.0); 
    }

  }

  Double_t scale = ((Double_t) nnx)/((Double_t) nx) * ((Double_t) nny)/((Double_t) ny);

  fsum->Rebin2D(nnx/nx,nny/ny); // rebin (fix problem on TF2 that cannot have nbin<4) ; warning rebin does not rescale width

  Double_t *us = new Double_t[nx];
  Double_t *vs = new Double_t[ny];
  
  for (Int_t j=0;j<nx;j++) {
    us[j]=0;
  }

  for (Int_t j=0;j<ny;j++) {
    vs[j]=0;
  }

  Int_t c1, c2;

  c1=0;
  c2=0;

  Int_t js=0;

  for (Int_t j=jB;j<=jT;j++) {
  
    Int_t mm=j-jB;
    
    if ((j%2)==1) {
      vs[js] = fsum->Integral(1,nx,mm+1,mm+1,"width")/scale; 
      js++;
      c1++;
    } else {
      for (Int_t i=iL;i<=iU;i++){
	Int_t kk=i-iL;
	us[kk] += fsum->Integral(kk+1,kk+1,mm+1,mm+1,"width")/scale;
      }
      c2++;
    }  

  }

  Float_t t0 = time_off + fRTime0+rnd.Gaus(fTriggerOffset, fTriggerJitter); // ... NEED Time of fly from GEANT4

  TSolGEMVStrip **virs;
  virs = new TSolGEMVStrip *[2];
  virs[0] = new TSolGEMVStrip(nx,fEleSamplingPoints);
  virs[1] = new TSolGEMVStrip(js,fEleSamplingPoints);

  virs[0]->setTime(t0);
  virs[0]->setHitCharge(fRTotalCharge);

  virs[1]->setTime(t0);
  virs[1]->setHitCharge(fRTotalCharge);

  Double_t pulse=0.;
  Double_t noisy_pulse;

  Short_t *dadc = new Short_t[fEleSamplingPoints];
  Int_t ai=0;
  Int_t posflag=0;

  Double_t cccsssx, cccsssy;
  cccsssx=0.;

  for (Int_t i=iL;i<=iU;i++) {
    posflag = 0;

    for (Int_t b=0;b<fEleSamplingPoints;b++) { // sampling

      pulse = TSolSimAux::PulseShape(t0+fEleSamplingPeriod*b, us[i-iL], fPulseShapeTau0, fPulseShapeTau1);
      noisy_pulse = rnd.Gaus(pulse, fPulseNoiseSigma); 
      dadc[b] = TSolSimAux::ADCConvert(noisy_pulse, fADCoffset, fADCgain, fADCbits);
      posflag += (Int_t) dadc[b];
    }

    if (posflag > 0) { // store only strip with signal -- do not work yet
      for (Int_t b=0;b<fEleSamplingPoints;b++) {
	virs[0]->AddSampleAt(dadc[b], b, ai);
      }

      virs[0]->AddStripAt(i,ai);

      virs[0]->AddChargeAt(us[i-iL], ai);

      cccsssx += us[i-iL];

      ai++;
    }

  }

  // y
  cccsssy=0.;

  virs[0]->setSize(ai);

  ai=0;
  for (Int_t j=0;j<js;j++) {
    posflag = 0;
    for (Int_t b=0;b<fEleSamplingPoints;b++) { // sampling

      pulse = TSolSimAux::PulseShape(t0+fEleSamplingPeriod*b, vs[j], fPulseShapeTau0, fPulseShapeTau1);
      noisy_pulse = rnd.Gaus(pulse, fPulseNoiseSigma); 
      dadc[b] = TSolSimAux::ADCConvert(noisy_pulse,fADCoffset, fADCgain, fADCbits);

      posflag += (Int_t) dadc[b];

    }

    if (posflag > 0) { // store only strip with signal -- do not work yet
      for (Int_t b=0;b<fEleSamplingPoints;b++) {
	virs[1]->AddSampleAt(dadc[b], b, ai);
      }

      virs[1]->AddStripAt(j+(Short_t) jB/2,ai);

      virs[1]->AddChargeAt(vs[j], ai);

      cccsssy += vs[j];

      ai++;

    }
  }

  virs[1]->setSize(ai);

  //

  delete fsum;

  delete[] dadc;
  delete[] us;
  delete[] vs;

  delete[] fRSNorm;
  delete[] fRCharge;
  delete[] fRX; 
  delete[] fRY;
 
  delete[] sgaus;

  return virs;

}

