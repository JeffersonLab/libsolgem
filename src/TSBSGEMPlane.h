#ifndef __TSBSGEMPLANE_H
#define __TSBSGEMPLANE_H

#include <cmath>

#include "THaSubDetector.h"
#include "TSBSBox.h"

#include "types.h"

class TClonesArray;

// In the present implementation we assume strips in a plane are
// uniform pitch and parallel to one another.

// In the SBS geometry, a plane is a box (as defined it TSBSBox class). 
// Refer to the comment in the header of TSBSBox class for more info.

// In the code, TSBSGEMPlane contains the information about the strips and such 

// It uses many of the methods from this class (namely the transformations methods).
// In addition to this it bears transformations to go from lab, spec, or box to the "strip" frame,
// Which is rotated wrt the box frame by the strip angle.

// TSBSGEMPlane also inherits form THaSubDetector, which grants it all the functions from its class
// (see http://hallaweb.jlab.org/podd/doc/html_v16/ClassIndex.html for more info).

// The "strip frame" is the box frame additionally rotated by the
// strip angle, so the x strips are parallel to y and measure position in x.

// The "projection frame" is a 1D coordinate system parallel to the x axis
// of the strip frame where 0 coincides with the the strip frame x coordinate 
// that runs through the lab origin.  This frame is equivalent for all boxes
// of the same wire  orientation regardless of position, so it is the frame
// to do tracking in

// The origin is specified in the spectrometer frame. The size is in the
// box frame.

class TSBSGEMPlane : public THaSubDetector {
    public:
        TSBSGEMPlane ();
	TSBSGEMPlane(const char *name, const char *desc,
		     THaDetectorBase* parent);
	virtual ~TSBSGEMPlane();
	
	//Read the geometry for the TSBSBox AND the strips parameters in the data base
	Int_t ReadDatabase (const TDatime& date);
	Int_t ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required = kFALSE);
	TClonesArray *GetClusters() { return fClusters; }

	Int_t Decode( const THaEvData &);
	TSBSBox& GetBox() const {return *fBox;};
	
	Int_t    GetNStrips()  const { return fNStrips; }// number of strips
	Double_t GetSPitch()   const { return fSPitch; } // pitch: distance between the middle of two strips
	Double_t GetSAngle()   const; // Angle (rad) between horizontal axis
	                              // in wedge frame
                                      // and normal to strips in dir of
	                              // increasing strip position
	Double_t GetSAngleComp() const { return 2*atan(1) - GetSAngle(); }

	// Frame conversions
	void LabToStrip (Double_t& x, Double_t& y, Double_t& z) const;//done
	// input and output in meters
	void LabToPlane (Double_t& x, Double_t& y, Double_t& z) const {
	  fBox->LabToBox (x, y, z);
	};  // input and output in meters
	void LabToSpec (Double_t& x, Double_t& y, Double_t& z) const {
	  fBox->LabToSpec (x, y, z);
	};  // input and output in meters
	void SpecToPlane (Double_t& x, Double_t& y) const {
	  fBox->SpecToBox (x, y);
	};  // input and output in meters
	void PlaneToStrip (Double_t& x, Double_t& y) const; // input and output in meters
	void SpecToStrip (Double_t& x, Double_t& y) const;  // input and output in meters
	void StripToSpec (Double_t& x, Double_t& y) const;  // input and output in meters
	void StripToPlane (Double_t& x, Double_t& y) const; // input and output in meters
	void PlaneToSpec (Double_t& x, Double_t& y) const {
	  fBox->BoxToSpec (x, y);
	};  // input and output in meters
	void SpecToLab (Double_t& x, Double_t& y, Double_t& z) const {
	  fBox->SpecToLab (x, y, z);
	};  // input and output in meters
	void PlaneToLab (Double_t& x, Double_t& y, Double_t& z) const {
	  fBox->BoxToLab (x, y, z);
	};  // input and output in meters
	void StripToLab (Double_t& x, Double_t& y, Double_t& z) const;//done

	Double_t StripNumtoStrip( Int_t num );

	Double_t StriptoProj( Double_t s );
	Double_t StripNumtoProj( Int_t s );

	// Edges of strip, in strip frame, in meters
	Double_t GetStripLowerEdge (UInt_t is) const;
	Double_t GetStripUpperEdge (UInt_t is) const;

        // Strip number corresponding to x-coordinate
        Int_t GetStripUnchecked( Double_t x )  const;
	Int_t GetStripInRange( Double_t x )    const;

	// Strip number corresponding to coordinates x, y in 
	// strip frame, or -1 if outside (2-d) bounds
	Int_t GetStrip (Double_t x, Double_t y) const;
	
	void Print( Option_t* opt="P" ) const;
	void SetRotations();

    private:
	TClonesArray  *fClusters; // Clusters

	Double_t fSAngle;        // Strip angle (measurement direction)
	Int_t    fNStrips;  // Number of strips
	Double_t fSPitch;   // Strip pitch (m)
	Double_t fSBeg;     // X coordinate of lower edge of first strip
	TSBSBox* fBox;  // Box geometry

	/* // Trig functions for rotations */
	/* // matrices for rotations */
	/* TMatrixD* fRotMat_SL;  */
	/* TMatrixD* fRotMat_LS;  */
	
	Double_t fCBS; // ... box to strip
	Double_t fSBS;
	
    public:
	ClassDef(TSBSGEMPlane,0)

};

// NB: I ignore why this is here and the rest is in the cxx file. 
// This is "inherithed" from TSolGEMPlane.
inline void
TSBSGEMPlane::PlaneToStrip (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCBS * x - fSBS * y;
  y = fSBS * temp + fCBS * y;
  return;
}

#endif//__TSBSGEMPLANE_H
