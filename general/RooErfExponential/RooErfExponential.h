#ifndef CHIBCHICPOLFW_GENERAL_ROOERFEXPONENTIAL_H__
#define CHIBCHICPOLFW_GENERAL_ROOERFEXPONENTIAL_H__

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include <cmath>
 
class RooErfExponential : public RooAbsPdf {
public:
  RooErfExponential() {} ;
  RooErfExponential(const char *name, const char *title,
                    RooAbsReal& _x,
                    RooAbsReal& _lambda,
                    RooAbsReal& _mu,
                    RooAbsReal& _sigma) :
    RooAbsPdf(name, title),
    x("x", "x", this, _x),
    lambda("lambda", "lambda", this, _lambda),
    mu("mu", "mu", this, _mu),
    sigma("sigma", "sigma", this, _sigma) {}

  RooErfExponential(const RooErfExponential &other, const char *name=0) :
    RooAbsPdf(other, name),
    x("x", this, other.x),
    lambda("lambda", this, other.lambda),
    mu("mu", this, other.mu),
    sigma("sigma", this, other.sigma) {}

  virtual TObject* clone(const char *newname) const { return new RooErfExponential(*this, newname); }

  inline virtual ~RooErfExponential() { }

protected:

  RooRealProxy x;
  RooRealProxy lambda;
  RooRealProxy mu;
  RooRealProxy sigma;

  Double_t evaluate() const {
    return exp(x * lambda) * (1 + TMath::Erf( (x - mu) / sigma) );
  }

private:

  ClassDef(RooErfExponential, 1)

};


#endif
