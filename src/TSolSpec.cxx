#include "TSolSpec.h"
#include "TSolGEMCluster.h"
#include "TSolGEMChamber.h"
#include <stdio.h>
#include <iostream>

using namespace std;

TSolSpec::TSolSpec(const char* name, const char* desc )
    :THaSpectrometer(name,desc) {

  // We don't need run db (not yet at least)
  fProperties &= ~kNeedsRunDB;
  return;
}

Int_t 
TSolSpec::AddGEM (TSolGEMChamber& pdet)
{
  // Add a detector to the internal lists of spectrometer detectors.
  // The detector object must be allocated and deleted by the caller.
  // Duplicate detector names are not allowed.

  fChambers.Add (&pdet);
  return 0;
}

Int_t TSolSpec::CoarseTrack(){
    // Needs work


//     int i,j;

//     i = j = 0;

//     // Assume decoding is done.  You can get the clustered hits with
//     // calles like

//     TSolGEMCluster *c;
//     Double_t z;

//     // loop over X planes
//     if( GetPlane(i)->GetDirection() == kGEMX ){
// 	c = (TSolGEMCluster *) GetPlane(i)->GetClusters()->AddrAt(j);
// 	z = GetPlane(i)->GetOrigin().Z();
	
// 	c->GetPos();
// 	c->GetE();
//     }

    return 0;
}

Int_t TSolSpec::CoarseReconstruct(){
    return 0;
}

Int_t TSolSpec::Track(){
    return 0;
}


Int_t TSolSpec::Reconstruct(){
    return 0;
}

Int_t TSolSpec::FindVertices(TClonesArray &a){
    return 0;
}

void
TSolSpec::Print()
{
  cout << "Hello, I'm a spectrometer named " << GetName() << endl;
	
  TIter next (&fChambers);
  while (TSolGEMChamber* theChamber =
	 static_cast<TSolGEMChamber*>(next()))
    theChamber->Print();
}
