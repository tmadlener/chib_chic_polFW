#ifndef CHIBCHICPOLFW_GENERAL_CALCANGLES_DEFS_H__
#define CHIBCHICPOLFW_GENERAL_CALCANGLES_DEFS_H__

#include "TVector3.h"

/** Reference frame enum. */
enum class RefFrame {
  HX, /**< helicity frame. */
  PX, /**< perpendicular helicity frame. */
  CS /**< collins soper frame. */
};

/** Small helper struct containing the angles of the positive muon in a reference frame */
struct Angles {
  double costh;
  double phi;
  double cosalpha;
};

/** Small helper struct for defining the axis of a reference frame. */
struct ReferenceAxis {
  TVector3 x;
  TVector3 y;
  TVector3 z;
};

#endif
