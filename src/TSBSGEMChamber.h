#ifndef __TSBSGEMCHAMBER_H
#define __TSBSGEMCHAMBER_H

#include "THaDetector.h"
#include "TSBSBox.h"

class THaEvData;
class TSBSGEMPlane;

// In the SBS geometry, a chamber is a box (as defined it TSBSBox class). 
// Refer to the comment in the header of TSBSBox class for more info.

// It uses many of the methods from this class (namely the transformations methods),
// which it completes by offering the options to use these methods
// with a TVector3 object as an input instead of 3 doubles.

// TSBSGEMChamber also inherits form THaDetector, which grants it all the functions from its class
// (see http://hallaweb.jlab.org/podd/doc/html_v16/ClassIndex.html for more info).

class TSBSGEMChamber : public THaDetector {
 public:
  //Constructors and destructor
  TSBSGEMChamber(const char *name, const char *desc);//It is recommended to use this constructor
  TSBSGEMChamber() : fPlanes(0), fNPlanes(0), fBox(0) {} // for ROOT RTTI

  virtual ~TSBSGEMChamber();
  
  //Read the geometry for the TSBSBox in the data base
  Int_t ReadDatabase (const TDatime& date);
  Int_t ReadGeometry (FILE* file, const TDatime& date, Bool_t required);
  const char* GetDBFileName() const;

  Int_t Decode( const THaEvData & );

  // Return positions of chamber edges, in its own reference frame, in meters
  Double_t GetLowerEdgeX() const {return -(GetSize())[0]/2;}
  Double_t GetLowerEdgeY() const {return -(GetSize())[1]/2;}
  Double_t GetUpperEdgeX() const {return +(GetSize())[0]/2;}
  Double_t GetUpperEdgeY() const {return +(GetSize())[1]/2;}

  TSBSBox& GetBox() const {return *fBox;};
  
  // Frame conversions
  void LabToPlane (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->LabToBox (x, y, z);
  };  // input and output in mm
  void PlaneToLab (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->BoxToLab (x, y, z);
  };  // input and output in mm

  void LabToSpec (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->LabToSpec (x, y, z);
  };  // input and output in mm
  void SpecToPlane (Double_t& x, Double_t& y, Double_t& z) const {
    z = z-fBox->GetD0()*1.0e3;
    fBox->SpecToBox (x, y);
  };  // input and output in mm
  void PlaneToSpec (Double_t& x, Double_t& y, Double_t& z) const {
    z = z+fBox->GetD0()*1.0e3;
    fBox->BoxToSpec (x, y);
  };  // input and output in mm
  void SpecToLab (Double_t& x, Double_t& y, Double_t& z) const {
    fBox->SpecToLab (x, y, z);
  };  // input and output in mm

  //Frame conversions with TVector3 objects as inputs
  void LabToPlane (TVector3& X_) const;
  void PlaneToLab (TVector3& X_) const;
  
  void LabToSpec (TVector3& X_) const;
  void SpecToPlane (TVector3& X_) const;
  void PlaneToSpec (TVector3& X_) const;
  void SpecToLab (TVector3& X_) const;
  
  // Access to the info of TSBSGEMPlane which is regarded as a subdetector of TSBSGEMChamber.
  // (see comments in the code of class TSBSGEMPlane)
   UInt_t GetNPlanes() const {return fNPlanes;};
    
  TSBSGEMPlane& GetPlane (UInt_t i) const {return *(fPlanes[i]);};
  Int_t InitPlane (const UInt_t i, const char* name, const char* desc);
  void Print (const Bool_t printplanes = kTRUE);

 private:
  TSBSGEMPlane** fPlanes; // List of planes
  UInt_t fNPlanes;
  TSBSBox* fBox;  // Box geometry

  ClassDef(TSBSGEMChamber,0)

};

#endif//__TSBSGEMCHAMBER_H
