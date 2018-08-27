#include "Riostream.h"

#include "RooDoubleCB.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <cmath>
#include "TMath.h"

ClassImp(RooDoubleCB)

RooDoubleCB::RooDoubleCB(const char *name,  const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _mu,
                        RooAbsReal& _sig,
                        RooAbsReal& _a1,
                        RooAbsReal& _n1,
                        RooAbsReal& _a2,
                        RooAbsReal& _n2) :
   RooAbsPdf(name, title),
   x("x",  "x",  this,  _x),
   mu("mu", "mu", this, _mu),
   sig("sig", "sig", this, _sig),
   a1("a1", "a1", this, _a1),
   n1("n1", "n1", this, _n1),
   a2("a2", "a2", this, _a2),
   n2("n2", "n2", this, _n2)
 { }


 RooDoubleCB::RooDoubleCB(const RooDoubleCB& other,  const char* name) :
   RooAbsPdf(other, name),
   x("x", this, other.x),
   mu("mu", this, other.mu),
   sig("sig", this, other.sig),
   a1("a1", this, other.a1),
   n1("n1", this, other.n1),
   a2("a2", this, other.a2),
   n2("n2", this, other.n2)
 { }



Double_t RooDoubleCB::evaluate() const
{
  const double u   = (x - mu) / sig;

  if (u > -a1 && u < a2) {
    return TMath::Exp(-0.5 * u*u);
  }

  if (u < -a1) {
    const double absAlpha1 = std::abs(a1);
    const double n1OverA1 = n1 / absAlpha1;
    const double A1 = TMath::Power(n1OverA1, n1) * TMath::Exp(-0.5 * a1*a1);
    const double B1 = n1OverA1 - absAlpha1;

    return A1 * TMath::Power(B1 - u, -n1);
  }

  if (u > a2) {
    const double absAlpha2 = std::abs(a2);
    const double n2OverA2 = n2 / absAlpha2;
    const double A2 = TMath::Power(n2OverA2, n2) * TMath::Exp(-0.5 * a2*a2);
    const double B2 = n2OverA2 - absAlpha2;

    return A2 * TMath::Power(B2 + u, -n2);
  }

  return 0;
}
