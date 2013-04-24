#include <iostream>

#include "TSolGEMData.h"

using namespace std;

TSolGEMData::TSolGEMData (UInt_t h)
{
    fGem  = 0;
    fEdep = 0;
    fTime = 0;
    fPID  = 0;
    fType = 0;
    fXi   = 0;
    fXo   = 0;
    fXr   = 0;
    fMom  = 0;
    fVert  = 0;
    fEntryNumber  = 0;

  if (h > 0) InitEvent (h);

  // Initialize variables
}

TSolGEMData::~TSolGEMData()
{
  ClearEvent();
}

void 
TSolGEMData::ClearEvent() 
{
  if (fGem)
  {
    delete[] fGem;
    delete[] fEdep;
    delete[] fTime;
    delete[] fPID;
    delete[] fType;
    delete[] fEntryNumber;
    
    for (UInt_t k = 0; k < fNHit; k++) 
    {
      delete fXi[k];
      delete fXo[k];
      delete fXr[k];
      delete fMom[k];
      delete fVert[k];
    }
    delete[] fXi;
    delete[] fXo;
    delete[] fXr;
    delete[] fMom;
    delete[] fVert;
    
    //    cerr << __FUNCTION__ << "fGEM deleted " << endl;
    
    fGem = 0;
    fNHit = 0;    
  }
};

void
TSolGEMData::InitEvent (UInt_t h)
{
  if (h <= 0)
    return;

  fNHit = h; // number of hits in event

  // Hits

  fGem = new UInt_t[h]; // chamber with hit
  fEdep = new Double_t[h]; // energy lost in drift
  fTime = new Double_t[h]; // Hit time
  fPID = new UInt_t[h]; // particle ID
  fType = new UInt_t[h]; // particle type
  fEntryNumber = new UInt_t[h]; // entry in file
  fXi = new TVector3*[h];
  fXo = new TVector3*[h];
  fXr = new TVector3*[h];
  fMom = new TVector3*[h];
  fVert = new TVector3*[h];
  for (UInt_t k = 0; k < fNHit; k++) 
  {
  	fXi[k] = new TVector3;
  	fXo[k] = new TVector3;
  	fXr[k] = new TVector3;
  	fMom[k] = new TVector3;
  	fVert[k] = new TVector3;
  }	  		
}

void 
TSolGEMData::AddGEMData(TSolGEMData *gd)
{
    // Create appropriately sized storage and copy

    UInt_t  Nnewhit = gd->GetNHit();

    UInt_t *newfGem    = new UInt_t[fNHit + Nnewhit];
    Double_t *newfEdep = new Double_t[fNHit + Nnewhit];
    Double_t *newfTime = new Double_t[fNHit + Nnewhit];
    UInt_t *newfPID    = new UInt_t[fNHit + Nnewhit];
    UInt_t *newfType   = new UInt_t[fNHit + Nnewhit];
    UInt_t *newfEntryNumber = new UInt_t[fNHit + Nnewhit];
    TVector3 **newfXi = new TVector3*[fNHit + Nnewhit];
    TVector3 **newfXo = new TVector3*[fNHit + Nnewhit];
    TVector3 **newfXr = new TVector3*[fNHit + Nnewhit];
    TVector3 **newfMom = new TVector3*[fNHit + Nnewhit];

    for (UInt_t k = 0; k < fNHit; k++){
	newfGem[k] = fGem[k];
	newfEdep[k] = fEdep[k];
	newfTime[k] = fTime[k];
	newfPID[k]  = fPID[k];
	newfType[k] = fType[k];
	newfEntryNumber[k] = fEntryNumber[k];

	newfXi[k]  = fXi[k];
	newfXo[k]  = fXo[k];
	newfXr[k]  = fXr[k];
	newfMom[k] = fMom[k];
    }	  		
    for (UInt_t k = 0; k < Nnewhit; k++){
	newfGem[k+fNHit] = gd->GetHitChamber(k);
	newfEdep[k+fNHit] = gd->GetHitEnergy(k);
	newfTime[k+fNHit] = gd->GetHitTime(k);
	newfPID[k+fNHit]  = gd->GetParticleID(k);
	newfType[k+fNHit] = gd->GetParticleType(k);
	newfEntryNumber[k+fNHit] = gd->GetEntryNumber(k);

	// Need to make new pointers so if se destroy the passed
	// GEMData we don't have dangling pointers
	newfXi[k+fNHit]  = new TVector3( gd->GetHitEntrance(k) );
	newfXo[k+fNHit]  = new TVector3( gd->GetHitExit(k) );
	newfXr[k+fNHit]  = new TVector3( gd->GetHitReadout(k) );
	newfMom[k+fNHit] = new TVector3( gd->GetMomentum(k) );
    }	  		

    delete[] fGem;
    delete[] fEdep;
    delete[] fTime;
    delete[] fPID;
    delete[] fType;
    delete[] fEntryNumber;
    delete[] fXi;
    delete[] fXo;
    delete[] fXr;
    delete[] fMom;
    
    fNHit += Nnewhit;
  fGem = newfGem; // chamber with hit
  fEdep = newfEdep; // energy lost in drift
  fTime = newfTime; // hit time
  fPID = newfPID; // particle ID
  fType = newfType; // particle type
  fEntryNumber = newfEntryNumber; // entry in file
  fXi = newfXi;
  fXo = newfXo;
  fXr = newfXr;
  fMom = newfMom;

  return;
}

void 
TSolGEMData::Print() const
{
  cout << "Run " << GetRun() << " Event " << GetEvent() << " " << GetNHit() << " hits" << endl;
}

void 
TSolGEMData::PrintHit (UInt_t k) const
{
  cout << "  Event " << GetEvent() << ":" << endl;
  cout << "    Momentum: " << GetMomentum(k).X()
       << " " << GetMomentum(k).Y() 
       << " " << GetMomentum(k).Z() 
       << " (mag. " << GetMomentum(k).Mag() << ")"
       << " MeV" << endl;
  cout << "    Hit time: " << GetHitTime(k) << " ns" << endl;
  cout << "    Hit entrance pos.: " << GetHitEntrance(k).X()
       << " " << GetHitEntrance(k).Y() 
       << " " << GetHitEntrance(k).Z() 
       << " mm" << endl;
  
  cout << "    Hit exit pos.: " << GetHitExit(k).X()
       << " " << GetHitExit(k).Y() 
       << " " << GetHitExit(k).Z() 
       << " mm" << endl;
  cout << "    Hit readout pos.: " << GetHitReadout(k).X()
       << " " << GetHitReadout(k).Y() 
       << " " << GetHitReadout(k).Z() 
       << " mm" << endl;

  cout << "    Hit vertex pos.: " << GetVertex(k).X()
       << " " << GetVertex(k).Y() 
       << " " << GetVertex(k).Z() 
       << " mm" << endl;
  cout << "    Hit energy: " << GetHitEnergy(k) << " eV" << endl;
  cout << "    Hit chamber: " << GetHitChamber(k) << endl;
  cout << "    Particle ID: " << GetParticleID(k) << endl;
  cout << "    Particle type: " << GetParticleType(k) << endl;
}
