#include <iostream>

#include "TSolGEMData.h"

using namespace std;

TSolGEMData::TSolGEMData (UInt_t h)
{
  fHitData.reserve(h);
}

TSolGEMData::~TSolGEMData()
{
}

void
TSolGEMData::ClearEvent() 
{
  fHitData.clear();
};

void
TSolGEMData::InitEvent (UInt_t h)
{
  if (h <= 0)
    return;

  fHitData.resize(h);
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
