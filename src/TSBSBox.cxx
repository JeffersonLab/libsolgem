#include "TSBSBox.h"

#include <iostream>
using namespace std;

#include <TDatime.h>
#include <TMath.h>

TSBSBox::TSBSBox (Double_t dmag, Double_t d0, Double_t xoffset, 
		  Double_t dx, Double_t dy, 
		  Double_t thetaH, Double_t thetaV)
  : fDMag(1),
    fD0(1),
    fXOffset(1),
    fDX(1),
    fDY(1),
    fThetaH(0),
    fThetaV(0)
{
  SetGeometry (dmag, d0, xoffset, dx, dy, thetaH, thetaV);
}

Bool_t
TSBSBox::Contains (Double_t x, Double_t y) const
{
  if(fabs(x)>fDX/2 || fabs(y)>fDY/2){
    return false;
  }else{
    return true;
  }
  
}

Bool_t
TSBSBox::Contains (Double_t x, Double_t y, Double_t z) const
{
  //tarnsform the point from the Lab to the Box coordinates
  LabToBox(x, y, z);
  
  return Contains(x, y);
}


void
TSBSBox::SetGeometry (const Double_t dmag,
		      const Double_t d0,
		      const Double_t xoffset,
		      const Double_t dx,
		      const Double_t dy,
		      const Double_t thetaH,
		      const Double_t thetaV)
{
  fDMag = dmag;
  fD0 = d0;
  fXOffset = xoffset;
  fDX = dx;
  fDY = dy;
  fThetaH = thetaH;
  fThetaV = thetaV;

  SetRotations();
   
  Double_t x0 = xoffset; 
  Double_t y0 = 0.0; 
  Double_t z0 = fD0; 
  
  SpecToLab(x0, y0, z0);
  
  // Evaluate the central point (in the lab) and the size of the Box
  fOrigin = TVector3(x0, y0, z0);
  //fSize = TVector3(fDX, fDY, 0.015955);
}

void
TSBSBox::HallCenterToBox (Double_t& x, Double_t& y, Double_t& z) const
{
  HallCenterToLab(x, y, z);
  LabToBox(x, y, z);
  return;
}

void
TSBSBox::HallCenterToSpec (Double_t& x, Double_t& y, Double_t& z) const
{
  HallCenterToLab(x, y, z);
  LabToSpec(x, y, z);
  return;
}

void
TSBSBox::LabToBox (Double_t& x, Double_t& y, Double_t& z) const
{
  LabToSpec(x, y, z);
  z = z-fD0*1.0e3;
  SpecToBox(x, y);
  
  return;
}

void
TSBSBox::HallCenterToLab (Double_t& x, Double_t& y, Double_t& z) const
{
  x = x + fDMag*sin(fThetaH)*1.0e3;
  // neutral in y
  z = z - fDMag*cos(fThetaH)*1.0e3;
  
  return;
}

void
TSBSBox::LabToSpec (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x, y, z};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_LB), m_temp);
  
  x = m_res(0, 0);
  y = m_res(1, 0);
  z = m_res(2, 0);
  
  return;
}

void
TSBSBox::SpecToBox (Double_t& x, Double_t& y) const
{
  x = x-fXOffset*1.0e3;
  //neutral for y
  
  return;
}

void
TSBSBox::BoxToSpec (Double_t& x, Double_t& y) const
{
  x = x+fXOffset*1.0e3;
  //neutral for y
  return;
}

void
TSBSBox::SpecToLab (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x, y, z};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_BL), m_temp);
  
  x = m_res(0, 0);
  y = m_res(1, 0);
  z = m_res(2, 0);
  
  return;
}

void
TSBSBox::LabToHallCenter (Double_t& x, Double_t& y, Double_t& z) const
{
  x = x - fDMag*sin(fThetaH)*1.0e3;
  // neutral in y
  z = z + fDMag*cos(fThetaH)*1.0e3;
  
  return;
}

void
TSBSBox::BoxToLab (Double_t& x, Double_t& y, Double_t& z) const
{
  BoxToSpec(x, y);
  z = z + fD0*1.0e3;
  SpecToLab(x, y, z);
  return;
}

void
TSBSBox::SpecToHallCenter (Double_t& x, Double_t& y, Double_t& z) const
{
  SpecToLab(x, y, z);
  LabToHallCenter(x, y, z);
  return;
}

void
TSBSBox::BoxToHallCenter (Double_t& x, Double_t& y, Double_t& z) const
{
  BoxToLab(x, y, z);
  LabToHallCenter(x, y, z);
  return;
}

void
TSBSBox::SetRotations()
{
  // arrays of variables for the matrices
  Double_t arr_roty0[9] = {cos(fThetaH),  0, sin(fThetaH),
			  0,             1,            0,
			  -sin(fThetaH), 0, cos(fThetaH)};
  Double_t arr_rotx1[9] = {1,            0,            0,
			   0, cos(fThetaV), -sin(fThetaV),
			   0, sin(fThetaV),  cos(fThetaV)};
  Double_t arr_rotz2[9] = {0, -1,  0,
  			   1,  0,  0,
  			   0,  0,  1};
  
  // the three following rotations are described in the lonc comment section 
  // in the class header file. 
  // Note that to obtain the rotation matrix from the box to the lab, 
  // the following matrices are multiplied in the reverse order they are declared.
  TMatrixD Roty0(3,3,arr_roty0);// rotation along hall pivot (y): spectrometer theta
  TMatrixD Rotx1(3,3,arr_rotx1);// rotation along x': spectrometer bending
  TMatrixD Rotz2(3,3,arr_rotz2);// rotation along z": box rotation
  TMatrixD Rotzx(3,3,arr_rotz2);
  fRotMat_LB = new TMatrixD(3,3, arr_roty0);
  
  Rotzx.Mult(Rotz2, Rotx1);
  fRotMat_LB->Mult(Rotzx, Roty0);// Box to Lab transformation
  
  // cout << " Lab -> Box " << endl;
  // fRotMat_LB->Print();
  
  fRotMat_BL = fRotMat_LB;
  fRotMat_BL->Invert();// Lab to Box transformation
  
  // cout << " Box -> Lab " << endl;
  // fRotMat_BL->Print();
}
