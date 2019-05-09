#ifndef CHIBCHICPOLFW_GENERAL_ROOPOWERLAWEXPONENTIAL_H__
#define CHIBCHICPOLFW_GENERAL_ROOPOWERLAWEXPONENTIAL_H__

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include <cmath>
 
class RooPowerlawExponential : public RooAbsPdf {
public:
  RooPowerlawExponential() {} ;
  RooPowerlawExponential(const char *name, const char *title,
                    RooAbsReal& _x,
                    RooAbsReal& _beta,
                    RooAbsReal& _alpha,
                    RooAbsReal& _q0) :
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    beta("beta", "beta", this, _beta),
    alpha("alpha", "alpha", this, _alpha),
    q0("q0", "q0", this, _q0) {}

  RooPowerlawExponential(const RooPowerlawExponential &other, const char *name=0) :
    RooAbsPdf(other, name),
    x("x", this, other.x),
    beta("beta", this, other.beta),
    alpha("alpha", this, other.alpha),
    q0("q0", this, other.q0) {}

  virtual TObject* clone(const char *newname) const { return new RooPowerlawExponential(*this, newname); }

  inline virtual ~RooPowerlawExponential() { }

protected:

  RooRealProxy x;
  RooRealProxy beta;
  RooRealProxy alpha;
  RooRealProxy q0;

  Double_t evaluate() const {
    const double delta = x - q0;
    // branchless determination of the sign of delta
    const int delta_sign = (0 < delta) - (delta < 0); // to force the pdf to negative values below cut off
    return delta_sign * std::pow(std::abs(delta), alpha) * exp(beta * delta);

  }

private:

  ClassDef(RooPowerlawExponential, 1)

};


#endif
