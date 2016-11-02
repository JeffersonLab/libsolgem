#include "TSBSBox.h"

#include <iostream>
using namespace std;

#include <TDatime.h>
#include <TMath.h>

TSBSBox::TSBSBox (Double_t d0, Double_t dx, Double_t dy, Double_t thetaH, Double_t thetaV)
  : fD0(1),
    fDX(1),
    fDY(1),
    fThetaH(0),
    fThetaV(0)
{
  SetGeometry (d0, dx, dy, thetaH, thetaV);
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
  // Is (x, y) within the box?
  
  LabToBox(x, y, z);
  
  if(fabs(x)>fDX/2 || fabs(y)>fDY/2){
    return false;
  }else{
    return true;
  }
  
}

void
TSBSBox::SetGeometry (const Double_t d0,
		      const Double_t dx,
		      const Double_t dy,
		      const Double_t thetaH,
		      const Double_t thetaV)
{
  fD0 = d0;
  fDX = dx;
  fDY = dy;
  fThetaH = thetaH;
  fThetaV = thetaV;

  SetRotations();
   
  fOrigin = TVector3(0.0, 0.0, 0.0);
  
  Double_t x0 = 0.0; 
  Double_t y0 = 0.0; 
  Double_t z0 = fD0; 
  
  LabToBox(x0, y0, z0);
  
  fOrigin = TVector3(x0, y0, z0);
  fSize = TVector3(fDX, fDY, 0.016);
}

void
TSBSBox::LabToBox (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x, y, z};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_LB), m_temp);
  
  x = m_res(0, 0) - fOrigin.X();
  y = m_res(1, 0) - fOrigin.Y();
  z = m_res(2, 0) - fOrigin.Z();
  
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

// void
// TSBSBox::SpecToBox (Double_t& x, Double_t& y) const
// {
//   Double_t r_temp[3] = {x, y, z};
//   TMatrixD m_temp(3, 1, r_temp);
  
//   TMatrixD m_res(3, 1, r_temp);
//   m_res.Mult((*fRotMat_LB), m_temp);
  
//   x = m_res(0, 0);
//   y = m_res(1, 0);
//   z = m_res(2, 0);
  
//   return;
// }

// void
// TSBSBox::BoxToSpec (Double_t& x, Double_t& y) const
// {
//   Double_t r_temp[3] = {x, y, z};
//   TMatrixD m_temp(3, 1, r_temp);
  
//   TMatrixD m_res(3, 1, r_temp);
//   m_res.Mult((*fRotMat_LB), m_temp);
  
//   x = m_res(0, 0);
//   y = m_res(1, 0);
//   z = m_res(2, 0);
  
//   return;
// }

void
TSBSBox::SpecToLab (Double_t& x, Double_t& y, Double_t& z) const
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
TSBSBox::BoxToLab (Double_t& x, Double_t& y, Double_t& z) const
{
  Double_t r_temp[3] = {x+fOrigin.X(), y+fOrigin.Y(), z+fOrigin.Z()};
  TMatrixD m_temp(3, 1, r_temp);
  
  TMatrixD m_res(3, 1, r_temp);
  m_res.Mult((*fRotMat_LB), m_temp);
  
  x = m_res(0, 0);
  y = m_res(1, 0);
  z = m_res(2, 0);
  
  return;
}

void
TSBSBox::SetRotations()
{
  Double_t arr_roty0[9] = {cos(fThetaH),  0, sin(fThetaH),
			  0,             1,            0,
			  -sin(fThetaH), 0, cos(fThetaH)};
  Double_t arr_rotx1[9] = {1,            0,            0,
			   0, cos(fThetaV), -sin(fThetaV),
			   0, sin(fThetaV),  cos(fThetaV)};
  Double_t arr_rotz2[9] = {0, -1,  0,
			   1,  0,  0,
			   0,  0,  1};
  
  TMatrixD Roty0(3,3,arr_roty0);// rotation along hall pivot (spectrometer theta)
  TMatrixD Rotx1(3,3,arr_rotx1);// ratotion along x': spectrometer bending
  TMatrixD Rotz2(3,3,arr_rotz2);// ratotion along x': spectrometer bending
  TMatrixD Rotyx(3,3,arr_roty0);
  fRotMat_BL = new TMatrixD(3,3, arr_roty0);

  // Set rotation angle trig functions
  Rotyx.Mult(Roty0, Rotx1);
  fRotMat_BL->Mult(Rotyx, Rotz2);
  
  fRotMat_LB = fRotMat_BL;
  fRotMat_LB->Invert();
}