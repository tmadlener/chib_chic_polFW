#ifndef CHIBCHICPOLFW_GENERAL_ROODOUBLECB_H__
#define CHIBCHICPOLFW_GENERAL_ROODOUBLECB_H__

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooDoubleCB : public RooAbsPdf {
public:
  RooDoubleCB() {};
  RooDoubleCB(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _mu,
              RooAbsReal& _sig,
              RooAbsReal& _a1,
              RooAbsReal& _n1,
              RooAbsReal& _a2,
              RooAbsReal& _n2);
  RooDoubleCB(const RooDoubleCB& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooDoubleCB(*this, newname); }
  inline virtual ~RooDoubleCB() { }

protected:

  RooRealProxy x;
  RooRealProxy mu;
  RooRealProxy sig;
  RooRealProxy a1;
  RooRealProxy n1;
  RooRealProxy a2;
  RooRealProxy n2;

  Double_t evaluate() const;

private:

  ClassDef(RooDoubleCB,1) // Your description goes here...
};

#endif
