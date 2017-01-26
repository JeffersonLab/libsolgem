#ifndef __TYPES_
#define __TYPES_

enum GEMDir_t { kGEMX, kGEMY, kGEMR, kGEMPhi };

enum StripSignalType{
  kPrimaryStrip = 0,
  kSecondaryStrip,
  kInducedStrip
};

struct GeoInfo{
    double d0;
    double xoffset;
    double dx;
    double dy;
    double thetaH;
    double thetaV;
    double z;
    double depth;
    double stripangle_u;
    double stripangle_v;
    double pitch_u;
    double pitch_v;
};

#endif//__TYPES_
