#ifndef __TYPES_
#define __TYPES_

enum GEMDir_t { kGEMX, kGEMY, kGEMR, kGEMPhi };

enum StripSignalType{
  kPrimaryStrip = 0,
  kSecondaryStrip,
  kInducedStrip
};

struct GeoInfo{
    double r0;
    double r1;
    double phi0;
    double dphi;
    double z;
    double depth;
    double stripangle_u;
    double stripangle_v;
    double pitch_u;
    double pitch_v;

};
#endif//__TYPES_
