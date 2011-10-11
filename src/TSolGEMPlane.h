#ifndef __TSOLGEMPLANE_H
#define __TSOLGEMPLANE_H

#include "THaSubDetector.h"

#include "types.h"

class TSolGEMCluster;
class TClonesArray;

// In the present implementation we assume strips in a plane are
// uniform pitch and parallel to one another.

// The geometry of a plane is characterized by: Origin (in lab frame;
// typically and always in the present implementation x and y are same
// as the parent chamber, if there is one), frame rotation angle (with 
// respect to lab x axis, same as parent chamber if there is one), size 
// (in rotated frame, same as parent chamber if there is one), 
// strip angle (of normal to strip's long axis, with 
// respect to x axis of rotated frame) and strip pitch (normal to strip's
// long axis).

// In the present implementation we assume strip angle is pi/2 (measuring
// along rotated x direction) or 0 (measuring along rotated y direction).

// The "chamber frame" (whether the plane is in a TSolGEMChamber or not)
// is the frame where the center of the plane is on the z axis and the
// x axis is rotated by the frame rotation angle with respect to the
// lab frame, so the edges of the frame are parallel to x and y.

// The "strip frame" is the chamber frame additionally rotated by the
// strip angle, so the strips are parallel to y and measure position in x.

// For planes that are in a 2-plane chamber, the paired plane is the 
// other plane in that chamber.

class TSolGEMPlane : public THaSubDetector {
    public:
        TSolGEMPlane ();
	TSolGEMPlane(const char *name, const char *desc,
		     THaDetectorBase* parent);
	virtual ~TSolGEMPlane() {;}

	Int_t ReadDatabase (const TDatime& date);
	Int_t ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required = kFALSE);
	TClonesArray *GetClusters() { return fClusters; }

	Int_t Decode( const THaEvData &);
	Double_t GetAngle() const {return fAngle;};
	GEMDir_t GetDirection() const { return fDir; }
	TSolGEMPlane *GetPairedPlane() { return fPairPlane; }

	Int_t    GetNStrips()  const { return fNStrips; }
	Double_t GetSPitch()   const { return fSPitch; } // in meters
	Double_t GetSAngle()   const; // Angle (rad) between horizontal axis
                                      // and normal to strips in dir of
	                              // increasing strip position
	Double_t GetSAngleComp() const { return 3.14159/2 - GetSAngle(); }

	// Frame conversions
	TVector3 LabToChamber (TVector3 v) const;  // input and output in meters
	TVector3 ChamberToStrip (TVector3 v) const; // input and output in meters
	TVector3 LabToStrip (TVector3 v) const {return ChamberToStrip (LabToChamber (v));}  // input and output in meters
	TVector3 StripToLab (TVector3 v) const {return ChamberToLab (StripToChamber (v));}  // input and output in meters
	TVector3 StripToChamber (TVector3 v) const;  // input and output in meters
	TVector3 ChamberToLab (TVector3 v) const;  // input and output in meters

	// Return positions of plane edges, in chamber frame, in meters
	Double_t GetLowerEdgeX() const {return (GetOrigin())[0] - (GetSize())[0];}
	Double_t GetLowerEdgeY() const {return (GetOrigin())[1] - (GetSize())[1];}
	Double_t GetUpperEdgeX() const {return (GetOrigin())[0] + (GetSize())[0];}
	Double_t GetUpperEdgeY() const {return (GetOrigin())[1] + (GetSize())[1];}

	// Edges of strip, in strip frame, in meters
	Double_t GetStripLowerEdge (UInt_t is) const {return (-(GetSize())[fDir] + is * GetSPitch());}
	Double_t GetStripUpperEdge (UInt_t is) const {return GetStripLowerEdge (is) + GetSPitch();}
	// Ends of strip, in strip frame, in meters
	Double_t GetStripLeftEdge (UInt_t is) const {return -(GetSize())[1-fDir];}
	Double_t GetStripRightEdge (UInt_t is) const {return (GetSize())[1-fDir];}

	// Strip number corresponding to coordinates x, y in 
	// strip frame, or -1 if outside (2-d) bounds
	Int_t GetStrip (Double_t x, Double_t y) const;

	void Print() const;

    private:
	TClonesArray  *fClusters; // Clusters
	Double_t fAngle;         // Angle of orientation (of frame)
	GEMDir_t fDir;		 // Strip orientation, x or y
	TSolGEMPlane *fPairPlane; // Paired plane

	Int_t    fNStrips;  // Number of strips
	Double_t fSPitch;   // Strip pitch (m)
	Double_t fSBeg;     // X coordinate of lower edge of first strip

    public:
	ClassDef(TSolGEMPlane,1)

};

#endif//__TSOLGEMPLANE_H
