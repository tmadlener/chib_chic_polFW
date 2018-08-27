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
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(n1/TMath::Abs(a1), n1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(n2/TMath::Abs(a2), n2)*TMath::Exp(-a2*a2/2);
  double B1  = n1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = n2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u, -n1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u, -n2);
  return result;
}
