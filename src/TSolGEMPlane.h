#ifndef __TSOLGEMPLANE_H
#define __TSOLGEMPLANE_H

#include <cmath>

#include "THaSubDetector.h"
#include "TSolWedge.h"

#include "types.h"

class TSolGEMCluster;
class TClonesArray;
class TSolWedge;

// In the present implementation we assume strips in a plane are
// uniform pitch and parallel to one another.

// A plane is a "wedge" (section of an annulus). It is characterized by
// the minimum and maximum radius (r0 and r1), the minimum angle (phi0),
// the angular width (dphi), the z coordinate of the front face (z0),
// the thickness in z (depth), strip angle (of normal to strip's long axis, 
// with respect to symmetry axis of wedge) and strip pitch (normal to strip's
// long axis).

// The inner and outer arcs are always centered on (x, y) = (0, 0) in the
// lab frame. 

// Typically and always in the present implementation, r0, r1, phi0, and
// dphi are the same as the parent chamber, if there is one.

// Derived from these quantities is a rectangular prism bounding box,
// one of whose sides is parallel to the symmetry axis of the wedge,
// described by a "size" which is half the transverse sizes and the
// full z size, and an "origin" which is the center of the front face of the 
// bounding box.

// The "wedge frame" is the frame whose x/y origin is the center of
// the bounding box and whose x axis lies along the symmetry axis of
// the wedge.

// The "strip frame" is the wedge frame additionally rotated by the
// strip angle, so the strips are parallel to y and measure position in x.

// The "projection frame" is a 1D coordinate system parallel to the x axis
// of the strip frame where 0 coincides with the the strip frame x coordinate 
// that runs through the lab origin.  This frame is equivalent for all wedges 
// of the same wire  orientation regardless of position, so it is the frame
// to do tracking in

// The origin and phi0 are specified in the lab frame. The size is in the
// wedge frame.

class TSolGEMPlane : public THaSubDetector {
    public:
        TSolGEMPlane ();
	TSolGEMPlane(const char *name, const char *desc,
		     THaDetectorBase* parent);
	virtual ~TSolGEMPlane();

	Int_t ReadDatabase (const TDatime& date);
	Int_t ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required = kFALSE);
	TClonesArray *GetClusters() { return fClusters; }

	Int_t Decode( const THaEvData &);
	TSolWedge& GetWedge() const {return *fWedge;};
	Double_t GetAngle() const {return fWedge->GetAngle();}; // rotation angle between lab and wedge frame

	Int_t    GetNStrips()  const { return fNStrips; }
	Double_t GetSPitch()   const { return fSPitch; } // in meters
	Double_t GetSAngle()   const; // Angle (rad) between horizontal axis
	                              // in wedge frame
                                      // and normal to strips in dir of
	                              // increasing strip position
	Double_t GetSAngleComp() const { return 2*atan(1) - GetSAngle(); }

	// Frame conversions
	void LabToPlane (Double_t& x, Double_t& y) const {fWedge->LabToWedge (x, y);};  // input and output in meters
	void PlaneToStrip (Double_t& x, Double_t& y) const; // input and output in meters
	void LabToStrip (Double_t& x, Double_t& y) const;  // input and output in meters
	void StripToLab (Double_t& x, Double_t& y) const;  // input and output in meters
	void StripToPlane (Double_t& x, Double_t& y) const;  // input and output in meters
	void PlaneToLab (Double_t& x, Double_t& y) const {fWedge->WedgeToLab (x, y);};  // input and output in meters

	Double_t StripNumtoStrip( Int_t num );

	Double_t StriptoProj( Double_t s );
	Double_t StripNumtoProj( Int_t s );

	// Edges of strip, in strip frame, in meters
	Double_t GetStripLowerEdge (UInt_t is) const;
	Double_t GetStripUpperEdge (UInt_t is) const;

	// Strip number corresponding to coordinates x, y in 
	// strip frame, or -1 if outside (2-d) bounds
	Int_t GetStrip (Double_t x, Double_t y) const;

	void Print() const;
	void SetRotations();

    private:
	TClonesArray  *fClusters; // Clusters

	Double_t fSAngle;        // Strip angle (measurement direction)
	Int_t    fNStrips;  // Number of strips
	Double_t fSPitch;   // Strip pitch (m)
	Double_t fSBeg;     // X coordinate of lower edge of first strip
	TSolWedge* fWedge;  // Wedge geometry

	// Trig functions for rotations
	Double_t fCLS; // cos lab to strip angle)
	Double_t fCWS; // ... wedge to strip
	Double_t fSLS; // sin...
	Double_t fSWS;
	
    public:
	ClassDef(TSolGEMPlane,1)

};

#endif//__TSOLGEMPLANE_H
