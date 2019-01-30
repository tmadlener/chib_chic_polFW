#ifndef CHIBCHICPOLFW_GENERAL_CALCANGLES_H__
#define CHIBCHICPOLFW_GENERAL_CALCANGLES_H__

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRotation.h"

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

/** calculate the angles of the lepInDilepFrame where the frame is defined by the passed axes. */
Angles calcAngles(const ReferenceAxis& refAxis,
                  const TLorentzVector& lepInDilepFrame)
{
  constexpr double overpi = 1 / M_PI;

  TRotation rotation;
  rotation.RotateAxes(refAxis.x, refAxis.y, refAxis.z);
  rotation.Invert(); // now transformin from "xyz" frame to the "new" frame

  auto lepInDilepRotated = lepInDilepFrame.Vect();
  lepInDilepRotated.Transform(rotation);

  const double cosTh = lepInDilepRotated.CosTheta();
  const double phi = lepInDilepRotated.Phi() * 180 * overpi;
  const double cosAlpha = std::sqrt(1 - cosTh*cosTh) * std::sin(lepInDilepRotated.Phi());

  return Angles{cosTh, phi, cosAlpha};
}

/**
 * Get the axis of the reference frame in the restframe of the onium.
 * b1 and b2 are the two beam directions in the onium rest frame.
 * oniaLab is the onium direction in the lab frame.
 */
ReferenceAxis determineRefFrameAxis(TVector3 b1, TVector3 b2, TVector3 oniaLab,
                                    const double oniaRap, const RefFrame refFrame)
{
  // y axis is the same in all frames
  const TVector3 yAxis = (b1.Cross(b2)).Unit() * (oniaRap < 0 ? -1 : 1); // rap depndent change of sign
  TVector3 zAxis;
  TVector3 xAxis;
  switch (refFrame) {
  case RefFrame::CS:
    zAxis = (b1 - b2).Unit(); // bisector of beams is z-axis
    break;
  case RefFrame::HX:
    zAxis = oniaLab;
    break;
  case RefFrame::PX:
    zAxis = (b1 - b2).Cross(yAxis).Unit(); // perpendicular to CS frame
    break;
  }

  xAxis = yAxis.Cross(zAxis);
  return ReferenceAxis{xAxis, yAxis, zAxis};
}

/**
 * calculate costh and phi in the given reference frame.
 * NOTE: taking the 4 momenta by value is done deliberately, since cloning would be necessary
 * otherwise.
 */
Angles calcAnglesInFrame(TLorentzVector muMinus, TLorentzVector muPlus,
                         const RefFrame refFrame)
{
  constexpr double pbeam = 4000; // GeV
  constexpr double Mprot = 0.9382720; // GeV
  constexpr double Ebeam = std::sqrt(pbeam*pbeam + Mprot*Mprot);
  TLorentzVector beam1(0, 0, pbeam, Ebeam);
  TLorentzVector beam2(0, 0, -pbeam, Ebeam);

  auto onia = muMinus + muPlus;
  const double oniaRap = onia.Rapidity();

  auto boostVecLabToOnia = -onia.BoostVector();
  beam1.Boost(boostVecLabToOnia);
  beam2.Boost(boostVecLabToOnia);
  muPlus.Boost(boostVecLabToOnia);

  auto beam1Dir = beam1.Vect().Unit();
  auto beam2Dir = beam2.Vect().Unit();
  auto oniaDir = onia.Vect().Unit();

  auto refFrameAxis = determineRefFrameAxis(beam1Dir, beam2Dir, oniaDir, oniaRap, refFrame);
  return calcAngles(refFrameAxis, muPlus);
}


/**
 * fold the angles as following:
 *  if phi >= -90 and phi < 0:
 *       phi = -1 * phi
 *  elif phi >= 90 and phi < 180::
 *      phi = 180 - phi
 *      costh = -1 * costh
 *  elif phi >= -180 and phi < -90:
 *      phi = 180 + phi
 *      costh = -1 * costh
 *
 */
Angles calcFoldAngles(const Angles& unfoldAngles)
{
  double foldPhi = unfoldAngles.phi;
  double foldCosth = unfoldAngles.costh;
  if (unfoldAngles.phi >= -90 && unfoldAngles.phi < 0) {
    foldPhi *= -1;
  } else if (unfoldAngles.phi >= 90 && unfoldAngles.phi < 180) {
    foldPhi = 180 - unfoldAngles.phi;
    foldCosth *= -1;
  } else if (unfoldAngles.phi >= -180 && unfoldAngles.phi < -90) {
    foldPhi = 180 + unfoldAngles.phi;
    foldCosth *= -1;
  }

  return Angles{foldCosth, foldPhi, unfoldAngles.cosalpha};
}

/**
 * Overload that takes only costh and phi as double values
 */
Angles calcFoldAngles(const double costh, const double phi) {
  return calcFoldAngles(Angles{costh, phi, 0});
}


#endif
