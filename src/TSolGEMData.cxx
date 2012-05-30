#include <iostream>

#include "TSolGEMData.h"

using namespace std;

TSolGEMData::TSolGEMData (UInt_t h)
{
    fGem  = 0;
    fEdep = 0;
    fPID  = 0;
    fType = 0;
    fXi   = 0;
    fXo   = 0;
    fXr   = 0;
    fMom  = 0;
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
    delete[] fPID;
    delete[] fType;
    delete[] fEntryNumber;
    
    for (UInt_t k = 0; k < fNHit; k++) 
    {
      delete fXi[k];
      delete fXo[k];
      delete fXr[k];
      delete fMom[k];
    }
    delete[] fXi;
    delete[] fXo;
    delete[] fXr;
    delete[] fMom;
    
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
  fPID = new UInt_t[h]; // particle ID
  fType = new UInt_t[h]; // particle type
  fEntryNumber = new UInt_t[h]; // entry in file
  fXi = new TVector3*[h];
  fXo = new TVector3*[h];
  fXr = new TVector3*[h];
  fMom = new TVector3*[h];
  for (UInt_t k = 0; k < fNHit; k++) 
  {
  	fXi[k] = new TVector3;
  	fXo[k] = new TVector3;
  	fXr[k] = new TVector3;
  	fMom[k] = new TVector3;
  }	  		
}

void 
TSolGEMData::Print()
{
  cout << "Run " << GetRun() << " Event " << GetEvent() << " " << GetNHit() << " hits" << endl;
}

void 
TSolGEMData::PrintHit (UInt_t k)
{
  cout << "  Event " << GetEvent() << ":" << endl;
  cout << "    Momentum: " << GetMomentum(k).X()
       << " " << GetMomentum(k).Y() 
       << " " << GetMomentum(k).Z() 
       << " MeV" << endl;
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
  cout << "    Hit energy: " << GetHitEnergy(k) << " eV" << endl;
  cout << "    Hit chamber: " << GetHitChamber(k) << endl;
  cout << "    Particle ID: " << GetParticleID(k) << endl;
  cout << "    Particle type: " << GetParticleType(k) << endl;
}
