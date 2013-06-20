#include "TSolSimGEMDigitization.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>

#include "TCanvas.h"
#include "TFile.h"
#include "TF2.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "TSolEVIOFile.h"  // needed for gendata class def
#include "TSolGEMData.h"
#include "TSolGEMVStrip.h"
#include "TSolSpec.h"
#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"
#include "TSolSimAux.h"
#include "TSolSimEvent.h"

using namespace std;

// Auxiliary class

TSolDigitizedPlane::TSolDigitizedPlane (Short_t nstrip,
					Short_t nsample) 
{
  fNOT = 0;
  fNSamples = nsample;

  fNStrips = nstrip;

  fType = new Short_t[nstrip];
  fCharge = new Float_t[nstrip];
  fTime = new Float_t[nstrip];
  fTotADC = new Short_t[nstrip];
  fOverThr = new UInt_t[nstrip];
  
  fPStripADC = new TArrayS(fNSamples*nstrip);
  
  if (fPStripADC==0) {
    cerr << __FUNCTION__ << " allocation failed" << endl;
    return;
  }

  fStripClusters.resize(fNStrips);

  Clear();
};

TSolDigitizedPlane::~TSolDigitizedPlane() 
{
  delete fPStripADC;
  delete[] fType;
  delete[] fCharge;
  delete[] fTime;
  delete[] fTotADC;
  delete[] fOverThr;
};

void
TSolDigitizedPlane::Clear()
{
  for (Int_t i = 0; i < fNStrips; i++) 
    {
      fTotADC[i]=0.;
      fType[i]=0;
      fCharge[i] = 0.;
      fTime[i] = 9999.;
    }
  fPStripADC->Reset();

  for( vector< std::vector<Short_t> >::iterator it = fStripClusters.begin();
       it != fStripClusters.end(); ++it ) {
    (*it).clear();
  }
}

void 
TSolDigitizedPlane::Cumulate (const TSolGEMVStrip *vv, Int_t type,
			      Short_t clusterID )
{
  Int_t j,k;
  Int_t idx;

  Short_t ooo,nnn;

  if (vv ) {
    for (j=0;j<vv->GetSize();j++) {
      idx = vv->GetIdx(j);
      if( idx < 0 ) continue; //  This is I believe should never happen?  SPR 5/30/2012
      fType[idx] |= (Short_t) type;
      fTime[idx] = (fTime[idx] < vv->GetTime()) ? fTime[idx] : vv->GetTime();
      fCharge[idx] += vv->GetCharge(j);
      for (k=0;k<fNSamples;k++) {
	ooo=fPStripADC->At(idx*fNSamples+k);
	nnn=vv->GetADC(j,k);
	fPStripADC->AddAt(ooo+nnn, idx*fNSamples+k);
	fTotADC[idx] += nnn;
      }
      fStripClusters[idx].push_back(clusterID);
    }
  }
};


UInt_t 
TSolDigitizedPlane::Threshold (Short_t thr) 
{

  fNOT = 0;

  for (UInt_t j = 0; j < (UInt_t) fNStrips; j++) 
    {
      if (fTotADC[j] > thr) 
	{
	  fOverThr[fNOT] = j;
	  fNOT++;	
	}
    }

  return fNOT;  
};


TSolSimGEMDigitization::TSolSimGEMDigitization( const TSolSpec& spect,
						const char* name )
  : THaAnalysisObject(name, "GEM simulation digitizer"),
    fDP(0), fNChambers(0), fNPlanes(0), fRX(0), fRY(0), fRSNorm(0), 
    fRCharge(0), fOFile(0), fOTree(0), fEvent(0)
{
  Init();
  Initialize (spect);

  fEvent = new TSolSimEvent(5);
}


TSolSimGEMDigitization::~TSolSimGEMDigitization()
{
  DeleteObjects();
}

void TSolSimGEMDigitization::DeleteObjects()
{
  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip)
	delete fDP[ic][ip];
      delete[] fDP[ic];
    }
  delete[] fDP;       fDP = 0;
  delete[] fNPlanes;  fNPlanes = 0;

  delete fOFile;      fOFile = 0;
  delete fOTree;      fOTree = 0;
  delete fEvent;      fEvent = 0;
}

void 
TSolSimGEMDigitization::Initialize(const TSolSpec& spect)
{
  // Initialize digitization structures based on parameters from given
  // spectrometer

  // Avoid memory leaks in case of reinitialization
  DeleteObjects();

  fNChambers = spect.GetNChambers();
  fDP = new TSolDigitizedPlane**[fNChambers];
  fNPlanes = new UInt_t[fNChambers];
  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      fNPlanes[ic] = spect.GetChamber(ic).GetNPlanes();
      fDP[ic] = new TSolDigitizedPlane*[fNPlanes[ic]];
      for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip) {
	fDP[ic][ip] =
	  new TSolDigitizedPlane( spect.GetChamber(ic).GetPlane(ip).GetNStrips());
      }
    }
  fEvCleared = false;
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
TSolSimGEMDigitization::Digitize (const TSolGEMData& gdata, const TSolSpec& spect) // digitize event 
{
  if (!fEvCleared)
    fEvent->Clear();
  UInt_t nh = gdata.GetNHit();

  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip)
	fDP[ic][ip]->Clear();
    }

  // ... NEED Time of flight from GEANT4
  Double_t trigger_time = fTrnd.Gaus(fTriggerOffset, fTriggerJitter);

  for (UInt_t ih = 0; ih < nh; ++ih)
    {
      UInt_t igem = gdata.GetHitChamber (ih);
      if (igem >= fNChambers)
	continue;
      
      Short_t itype = (1 << gdata.GetParticleType(ih)); // signal = 1, bck = 2, 4, 8 ...
	
      // Only use primary particle here
//      if( gdata.GetParticleID(ih) != 1 ) continue;

      TVector3 vv1 = gdata.GetHitEntrance (ih);
      TVector3 vv2 = gdata.GetHitExit (ih);
      //TVector3 vv3 = gdata.GetHitReadout (ih);
      //  FIXME:  Sometimes we don't always have the readoutplane data
      //          which is a problem when we try to digitize.  However,
      //          the distance is fixed from the initial hit plane so 
      //          I'm hardcoding it here for now.  It shoudl go in the 
      //          database
      TVector3 vv3 = gdata.GetHitEntrance(ih) + TVector3(0.0, 0.0, 4.5975+4.5925);

      // These vectors are in the lab frame, we need them in the chamber frame
      // Also convert to mm

      TVector3 offset = spect.GetChamber(igem).GetOrigin() * 1000.0;
      Double_t angle = spect.GetChamber(igem).GetAngle();
      vv1 -= offset;
      vv2 -= offset;
      vv3 -= offset;
      vv1.RotateZ (-angle);
      vv2.RotateZ (-angle);
      vv3.RotateZ (-angle);
	
      TSolGEMVStrip **dh = NULL;
      IonModel (vv1, vv2, gdata.GetHitEnergy(ih), vv3);
      if (fRNIon > 0) 
	{
	  // Time of this hit's avalance relative to the trigger
	  Double_t time_zero = fRTime0*1e9 - trigger_time;
	  itype = 1; // itype 1 is no time randomization, itype 2 has time randomization
	  time_zero += (itype == 1) ? 0.
	    // randomization of the bck ( assume 3 useful samples at 25 ns)
	    : fTrnd.Uniform (fGateWidth + 75.) - fGateWidth;

	  dh = AvaModel (igem, spect, vv1, vv2, time_zero);
	}
      Short_t id = SetTreeHit (ih, spect, dh, gdata);
      if (dh != NULL)
	{
	  for (UInt_t j = 0; j < 2; j++) 
	    {
	      fDP[igem][j]->Cumulate (dh[j], itype, id );
	    }
	  // TODO: make dh[2] a member variable & clear it here to avoid the constant
	  // construction and deletion
	  delete dh[0];
	  delete dh[1];
	  delete[] dh;
	} 
    }
  SetTreeStrips (gdata);
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
#define DBG_ION 0

  Double_t LL;
  Double_t deltaE=elost; // eV  MC

  TVector3 vseg = xo-xi; // mm

  // DEBUG  TRandom3 rnd(0);
  TRandom3& rnd = fTrnd;

  // ---- extract primary ions from Poisson

  fRNIon = rnd.Poisson(deltaE/fGasWion);

  if (fRNIon <=0)
    return;

#if DBG_ION > 0
  cout << "E lost = " << elost << ", " << fRNIon << " ions";
#endif
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

#if DBG_ION > 1
    printf("z coords %f %f %f %f lion %f LL %lf\n", 
	   xi.Z(), xo.Z(), xrout.Z(), vseg.Z(), lion, LL);
    printf("ttime = %e\n", ttime);
#endif

    fRTime0 = (ttime<fRTime0) ? ttime : fRTime0; // minimum traveling time [s]

    gnorm = fGainMean/sqrt(fGain0); // overall gain TBC

    fRSNorm[i] = sqrt(2.*fGasDiffusion*ttime); // spatial sigma on readout

    switch (fAvalancheChargeStatistics) {
    case 1:
      fRCharge[i]= rnd.Gaus(fGainMean, gnorm); // Gaussian distribution of the charge
#if DBG_ION > 0
      cout << "x, y = " << fRX[i] << ", " << fRY[i] << " snorm = "
	   << fRSNorm[i] << " charge " << fRCharge[i] << endl;
#endif
      break;
    default: 
      fRCharge[i]= rnd.Exp(fGainMean); // Furry distribution
#if DBG_ION > 0
      cout << "x, y = " << fRX[i] << ", " << fRY[i] << " snorm = " << fRSNorm[i] 
	   << " charge " << fRCharge[i] << endl;
#endif
      break;
    }
#if DBG_ION > 0
    cout << "fRTime0 = " << fRTime0 << endl;
#endif

    fRTotalCharge += fRCharge[i];

    fRSMax = (fRSNorm[i]>fRSMax) ? fRSNorm[i] : fRSMax;
  }

}

//.......................................................
// avalanche model
//

TSolGEMVStrip **
TSolSimGEMDigitization::AvaModel(const Int_t ic,
				 const TSolSpec& spect,
				 const TVector3& xi, 
				 const TVector3& xo,
				 const Double_t t0)
{
#define DBG_AVA 0

#if DBG_AVA > 0
  cout << "Chamber " << ic << "----------------------------------" << endl;
  cout << "In  " << xi.X() << " " << xi.Y() << " " << xi.Z() << endl;
  cout << "Out " << xo.X() << " " << xo.Y() << " " << xo.Z() << endl;
#endif

  // xi, xo are in chamber frame, in mm

  Double_t nsigma = fAvalancheFiducialBand; // coverage factor

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

  Double_t glx = spect.GetChamber(ic).GetLowerEdgeX() * 1000.0;
  Double_t gly = spect.GetChamber(ic).GetLowerEdgeY() * 1000.0;
  Double_t gux = spect.GetChamber(ic).GetUpperEdgeX() * 1000.0;
  Double_t guy = spect.GetChamber(ic).GetUpperEdgeY() * 1000.0;

   // fprintf(stderr, "%s %s:  out of sector, chamber %d\nFollowing relations should hold:\n(x1 %f>glx %f) (x0 %f<gux %f)\n(y1 %f>gly %f) (y0 %f<guy %f)\n\n", __FILE__,  __FUNCTION__, ic, x1, glx, x0, gux, y1, gly, y0, guy );

  if (x1<glx || x0>gux ||
      y1<gly || y0>guy) { // out of active area of the sector
    delete[] fRSNorm;
    delete[] fRCharge;
    delete[] fRX; 
    delete[] fRY; 
    cerr << "out of sector" << endl;
    fprintf(stderr, "%s %s:  out of sector, chamber %d\nFollowing relations should hold:\n(x1 %f>glx %f) (x0 %f<gux %f)\n(y1 %f>gly %f) (y0 %f<guy %f)\n", __FILE__,  __FUNCTION__, ic, x1, glx, x0, gux, y1, gly, y0, guy );
    fprintf(stderr, "\tplane %d r %f phi %f\n\n", (ic-1)/30+1, sqrt(x0*x0+y0*y0), atan(y0/x0)*180/3.14159  );
    return 0;
  }

  if(x0<glx) x0=glx;
  if(y0<gly) y0=gly;
  if(x1>gux) x1=gux;
  if(y1>guy) y1=guy;

  // Loop over chamber planes

  TSolGEMVStrip **virs;
  virs = new TSolGEMVStrip *[fNPlanes[ic]];
  for (UInt_t ipl = 0; ipl < fNPlanes[ic]; ++ipl)
    {
#if DBG_AVA > 0
      cout << "plane " << ipl << " =========================" << endl;
#endif	

      // Compute strips affected by the avalanche

      TSolGEMPlane* pl = &(spect.GetChamber(ic).GetPlane(ipl)); 

      // Positions in strip frame
      Double_t xs0 = x0 * 1e-3; Double_t ys0 = y0 * 1e-3;
      pl->PlaneToStrip (xs0, ys0);
      xs0 *= 1e3; ys0 *= 1e3;
      Double_t xs1 = x1 * 1e-3; Double_t ys1 = y1 * 1e-3;
      pl->PlaneToStrip (xs1, ys1);
      xs1 *= 1e3; ys1 *= 1e3;

#if DBG_AVA > 0
      cout << "xs0 ys0 xs1 ys1 " << xs0 << " " << ys0 << " " << xs1 << " " << ys1 << endl;
#endif

      Int_t iL = pl->GetStrip (xs0 * 1e-3, ys0 * 1e-3);
      Int_t iU = pl->GetStrip (xs1 * 1e-3, ys1 * 1e-3);
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
#if DBG_AVA > 0
      cout << "iL gsle " << iL << " " << pl->GetStripLowerEdge (iL) << endl;
      cout << "iU gsle " << iU << " " << pl->GetStripLowerEdge (iU) << endl;
#endif
      Double_t xl = pl->GetStripLowerEdge (iL) * 1000.0;
      Double_t xr = pl->GetStripUpperEdge (iU) * 1000.0;

      // Limits in y are y limits of track plus some reasonable margin
      // We do this in units of strip pitch for convenience (even though
      // this is the direction orthogonal to the pitch direction)

      Double_t pitch = pl->GetSPitch() * 1000.0;
      // Should the .1 be here?? ???????????????????????????????????????
      Double_t yq = pitch * .1;
      // Double_t yq = pitch;
      Double_t yb = yq * TMath::Floor (ys0 / yq);
      Double_t yt = yq * TMath::Floor (ys1 / yq);
      if (yb > yt)
	{
	  Double_t t = yb;
	  yb = yt;
	  yt = t;
	}

      // # of coarse histogram bins: 1 per pitch each direction

      Int_t nx = iU - iL + 1;
      Int_t ny = (Int_t) ((yt - yb + 1) / yq + 0.5);
#if DBG_AVA > 0
      cout << "xr xl yt yb nx ny " 
	   << xr << " " << xl << " " << yt << " " << yb 
	   << " " << nx << " " << ny << endl;
#endif
      // define function, gaussian and sum of gaussian

      Double_t* fsuma = new Double_t[nx*ny];
      Double_t xbw = (xr - xl) / nx;
      Double_t ybw = (yt - yb) / ny;
#if DBG_AVA > 0
      cout << "xbw ybw " << xbw << " " << ybw << endl;
#endif
      memset (fsuma, 0, nx * ny * sizeof (fsuma));
      for (UInt_t i = 0; i < fRNIon; i++) 
	{
	  Double_t ggnorm = fRCharge[i] / 3.14 / 9. / fRSNorm[i] / fRSNorm[i]; // normalized to charge
	  Double_t frxs = fRX[i] * 1e-3; Double_t frys = fRY[i] * 1e-3;
	  pl->PlaneToStrip (frxs, frys);
	  frxs *= 1e3; frys *= 1e3;
	  // bin containing center and # bins each side to process
	  Int_t ix = (frxs-xl) / xbw;
	  Int_t iy = (frys-yb) / ybw;
	  Int_t dx = 3 * fRSNorm[i] / xbw  + 1;
	  Int_t dy = 3 * fRSNorm[i] / ybw  + 1;
#if DBG_AVA > 1
	  cout << "ix dx iy dy " << ix << " " << dx << " " << iy << " " << dy << endl;
#endif
	  Double_t r2 = pow (3 * fRSNorm[i], 2);
	  // Loop over bins
	  // xc and yc are center of current bin
	  Double_t xc = xl + (ix - dx + 0.5) * xbw;
	  for (Int_t jx = ix-dx; jx <= ix+dx; ++jx)
	    {
	      Double_t xd2 = pow (frxs-xc, 2);
	      Double_t yc = yb + (iy - dy + 0.5) * ybw;
	      for (Int_t jy = iy-dy; jy <= iy+dy; ++jy)
		{
		  if (jx >= 0 && jx < nx && jy >= 0 && jy < ny && xd2 + pow(frys-yc, 2) <= r2)
		    fsuma[jx*ny+jy] += ggnorm;
		  yc += ybw;
	        }
	      xc += xbw;
	    }
	}

      TH2F *fsum = new TH2F ("fsum", "", nx, xl, xr, ny, yb, yt);
      for (Int_t jx = 1; jx <= nx; ++jx)
	for (Int_t jy = 1; jy <= ny; ++jy)
	  {
	    fsum->SetBinContent (jx, jy, fsuma[(jx-1)*ny+jy-1]);
	  }
      delete[] fsuma;

      Double_t *us = new Double_t[nx];
      for (Int_t j = 0; j < nx; j++) 
	us[j] = 0;

      for (Int_t i = iL; i <= iU; i++)
	{
	  Int_t kk = i - iL;
	  us[kk] += fsum->Integral (kk+1, kk+1, 1, ny, "width");
	}

#if DBG_AVA > 0      
      cout << "t0 = " << t0
	   << endl;
#endif
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
	      pulse = TSolSimAux::PulseShape (fEleSamplingPeriod * b - t0, 
					      us[i-iL], 
					      fPulseShapeTau0, 
					      fPulseShapeTau1);
	      noisy_pulse = fTrnd.Gaus (pulse, fPulseNoiseSigma); 
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
      
      delete[] dadc;
      delete[] us;
      
      // delete[] sgaus;
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
TSolSimGEMDigitization::PrintCharges() const
{
  cout << " Chb  Pln  Strip  Typ    ADC    Charge      Time\n";
  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip)
	for (UInt_t ist = 0; ist < (UInt_t) GetNStrips(ic, ip); ++ist)
	  {
	    if (fDP[ic][ip]->GetCharge (ist) > 0)
	      cout << setw(4) << ic
		   << " " << setw(4) << ip
		   << " " << setw(6) << ist
		   << " " << setw(4) << GetType (ic, ip, ist)
		   << " " << setw(6) << GetTotADC (ic, ip, ist)
		   << fixed << setprecision(1)
		   << " " << setw(9) << GetCharge (ic, ip, ist)
		   << fixed << setprecision(3)
		   << " " << setw(9) << GetTime (ic, ip, ist)
		   << endl;
	  }
    }
}


void
TSolSimGEMDigitization::PrintSamples() const
{
  cout << " Chb  Pln  Strip  Typ    ADC samples \n";
  for (UInt_t ic = 0; ic < fNChambers; ++ic)
    {
      for (UInt_t ip = 0; ip < fNPlanes[ic]; ++ip)
	for (UInt_t ist = 0; ist < (UInt_t) GetNStrips (ic, ip); ++ist)
	  if (GetCharge (ic, ip, ist) > 0)
	    {
	      cout << setw(4) << ic
		   << " " << setw(4) << ip
		   << " " << setw(6) << ist
		   << " " << setw(4) << GetType (ic, ip, ist);
	      for (UInt_t is = 0; is < (UInt_t) GetNSamples (ic, ip); ++is)
		cout << " " << setw(6) << GetADC (ic, ip, ist, is);
	      cout << endl;
	    }
    }
}

// Tree methods
void 
TSolSimGEMDigitization::InitTree (const TSolSpec& spect, const TString& ofile)
{
  fOFile = new TFile( ofile, "RECREATE");

  if (fOFile == 0 || fOFile->IsZombie() ) 
    {
      cerr << "Error: cannot open output file " << ofile << endl;
      delete fOFile; fOFile = 0;
      return;
    } 

  fOTree = new TTree( treeName, "Tree of digitized values");

  // create the tree variables

  fOTree->Branch( eventBranchName, "TSolSimEvent", &fEvent );
  
 }

inline
static void ChamberToSector( Short_t chamber, Short_t& sector, Short_t& plane )
{
  // Conversion from chamber index to sector/plane indices.
  // The meaning of the chamber index is not defined anywhere else this code.
  // The only requirement is that the database for TSolSpec matches whatever
  // is defined in the MC used to generate the input.
  // The database floating around at this time (February 2013) uses the definition
  // ich = is + nsectors*ipl (is = sector, ipl = plane).
  // The number of sectors is implied to be 30.

  static const int NSECTORS = 30;
  div_t d = div( chamber, NSECTORS );
  sector = d.rem;
  plane  = d.quot;
}

void
TSolSimGEMDigitization::SetTreeEvent (const TSolGEMData& tsgd,
				      const TSolEVIOFile& f, Int_t evnum )
{
  // Set overall event info.

  fEvent->Clear();
  fEvCleared = true;
  fEvent->fRunID = tsgd.GetRun();
  fEvent->fEvtID = (evnum < 0) ? tsgd.GetEvent() : evnum;
  for( UInt_t i=0; i<f.GetNGen(); ++i ) {
    const gendata* gd = f.GetGenData(i);
    //TODO: get GEANT id?
    fEvent->AddTrack( i+1, gd->GetPID(), gd->GetWeight(),
		      gd->GetV()*1e-3, // Vertex coordinates in [m]
		      gd->GetP()*1e-3  // Momentum in [GeV]
		      );
  }
  // FIXME: either only one GenData per event, or multiple weights per event
  if( f.GetNGen() > 0 ) 
    fEvent->fWeight = f.GetGenData(0)->GetWeight();
}

Short_t
TSolSimGEMDigitization::SetTreeHit (const UInt_t ih,
				    const TSolSpec& spect,
				    TSolGEMVStrip* const *dh,
				    const TSolGEMData& tsgd)
{
  // Sets the variables in fEvent->fGEMClust describing a hit
  // This is later used to fill the tree.

  TSolSimEvent::GEMCluster clust;

  cout << "here1" << endl;
  UInt_t igem = tsgd.GetHitChamber(ih);
  ChamberToSector( igem, clust.fSector, clust.fPlane );
  //  clust.fRefEntry = tsgd.GetEntryNumber(ih);  // Apparently never initialized
  clust.fType     = tsgd.GetParticleID(ih);   // GEANT particle counter
  clust.fPID      = tsgd.GetParticleType(ih); // PDG PID
  clust.fP        = tsgd.GetMomentum(ih)    * 1e-3; // [GeV]
  clust.fXEntry   = tsgd.GetHitEntrance(ih) * 1e-3; // [m]
  // The best estimate of the "true" hit position is the center of the
  // ionization region
  clust.fMCpos    = (tsgd.GetHitEntrance(ih)+tsgd.GetHitExit(ih)) * 5e-4; // [m]

  // Calculate hit position in the Tracker frame. This is fMCpos relative to
  // the origin of first plane of the sector, but rotated by the nominal
  // (non-offset) sector angle.
  // NB: assumes NSECTORS=30, even spacing, clockwise numbering
  Double_t sector_angle = 12.*clust.fSector*TMath::DegToRad();
  clust.fHitpos = clust.fMCpos - spect.GetChamber(clust.fSector).GetOrigin();
  clust.fHitpos.RotateZ(-sector_angle);

  const TSolGEMChamber& ch = spect.GetChamber(igem);
  for (UInt_t j = 0; j < 2; j++) {
    if (dh != NULL && dh[j] != NULL)
      {
	clust.fSize[j]  = dh[j]->GetSize();
	clust.fStart[j] = (clust.fSize[j] > 0) ? dh[j]->GetIdx(0) : -1;
      }
    else
      {
	clust.fSize[j] = 0;
	clust.fStart[j] = 0;
      }
    const TSolGEMPlane& pl = ch.GetPlane(j);
    Double_t proj_angle = pl.GetAngle() + pl.GetSAngle() - sector_angle;
    TVector3 hitpos(clust.fHitpos);
    hitpos.RotateZ(-proj_angle);
    clust.fXProj[j] = hitpos.X();
  }
  if (dh != NULL && dh[0] != NULL)
    {
      clust.fCharge = dh[0]->GetHitCharge();
      clust.fTime   = dh[0]->GetTime();
    }
  else
    {
      clust.fCharge = 0;
      clust. fTime = 0;
    }
  clust.fID     = fEvent->fGEMClust.size()+1;
  clust.fVertex = tsgd.GetVertex (ih);

  fEvent->fGEMClust.push_back( clust );

  if( clust.fPlane == 0 && clust.fType == 1 )
    fEvent->fNSignal++;
  cout << "here9" << endl;

  return clust.fID;
  // after the return??
  //  fEvCleared = false;
}

void
TSolSimGEMDigitization::SetTreeStrips (const TSolGEMData& tsgd)
{
  // Sets the variables in fEvent->fGEMStrips describing strip signals
  // This is later used to fill the tree.

  fEvent->fGEMStrips.clear();
  
  TSolSimEvent::DigiGEMStrip strip;
  for (UInt_t ich = 0; ich < GetNChambers(); ++ich)
    {
      ChamberToSector( ich, strip.fSector, strip.fPlane );

      // The "plane" here is actually the projection (= readout coordinate).
      // TSolGEMChamber::ReadDatabase associates the name suffix "x" with
      // the first "plane", and "y", with the second. However, strip angles can
      // be different for each chamber, so what's "x" in one chamber may very
      // well be something else, like "x'", in another. These angles don't even
      // have to match anything in the Monte Carlo.
      for (UInt_t ip = 0; ip < GetNPlanes (ich); ++ip)
	{
	  strip.fProj = (Short_t) ip; 
	  strip.fNsamp = TMath::Min((Short_t)MC_MAXSAMP,
				    (Short_t)GetNSamples(ich, ip));
	  UInt_t nover = Threshold (ich, ip, 0); // threshold is zero for now
	  for (UInt_t iover = 0; iover < nover; iover++) 
	    {
	      UInt_t idx = GetIdxOverThr(ich, ip, iover);
	      strip.fChan = (Short_t) idx;

	      for (Int_t ss = 0; ss < strip.fNsamp; ++ss)
		strip.fADC[ss] = (Short_t) GetADC(ich, ip, idx, ss);

	      strip.fSigType = (Short_t) GetType(ich, ip, idx);
	      strip.fCharge  = GetCharge(ich, ip, idx);
	      strip.fTime1   = GetTime(ich, ip, idx);

	      const vector<Short_t>& sc = GetStripClusters(ich, ip, idx);
	      strip.fClusters.Set( sc.size(), &sc[0] );

	      fEvent->fGEMStrips.push_back( strip );
	    }
	} 
    }
  fEvCleared = false;
}

void
TSolSimGEMDigitization::FillTree ()
{
  if (fOFile && fOTree 
      // added this line to not write events where there are no entries

      // Remove for background study
      //      && fEvent->fGEMStrips.size() > 0 && fEvent->fGEMClust.size() > 0

      )
    {
      fOFile->cd();
      fOTree->Fill();
    }
}

void 
TSolSimGEMDigitization::WriteTree () const
{
  if (fOFile && fOTree) {
    fOFile->cd();
    fOTree->Write();
  }
}

void 
TSolSimGEMDigitization::CloseTree () const
{
  if (fOFile) fOFile->Close();
}

