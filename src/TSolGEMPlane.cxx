#include "TSolGEMPlane.h"

#include <iostream>

#include "TClonesArray.h"
#include "TSolGEMCluster.h"
#include "THaEvData.h"

using namespace std;

TSolGEMPlane::TSolGEMPlane()
  : THaSubDetector()
{
  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  return;
}

Int_t 
TSolGEMPlane::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  return kOK;
}

Int_t 
TSolGEMPlane::ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required)
{
  THaSubDetector::ReadGeometry (file, date, false);

  const DBRequest request[] = 
    {
      {"direction",   &fDir,         kInt,    0, 1},
      {"nstrips",     &fNStrips,     kInt,    0, 1},
      {"stripbegin",  &fSBeg,        kDouble, 0, 1},
      {"strippitch",  &fSPitch,      kDouble, 0, 1},
      {"pixelfactor", &fPixelFactor, kInt,    0, 1},
      {0}
    };
  Int_t err = LoadDB( file, date, request, fPrefix );

  if (err)
    return err;

  return kOK;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

    int i = 0;

    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}

Double_t 
TSolGEMPlane::GetSAngle()   const
{
  if (fDir == kGEMX)
    return 3.14159/2;
  else if (fDir == kGEMY)
    return 0.0;
  else
    {
      cerr << __FUNCTION__ << " Strip angle undefined" << endl;
      return 9999.0;
    }
}

void 
TSolGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2) << endl;

  const Float_t* s = GetSize();
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;
  cout << "  " << GetNStrips() << " strips, pitch " << GetSPitch() << endl;
  cout << "  " << GetNPixels() << " pixels, pitch " << GetPPitch() << endl;
  cout << "  Strips begin at " << GetSBeg() << ", angle is " << GetSAngle() << endl;
}
