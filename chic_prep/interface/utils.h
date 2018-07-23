#ifndef CHIBCHICPOLFW_CHICPREP_UTILS_H__
#define CHIBCHICPOLFW_CHICPREP_UTILS_H__

#include "general/interface/calcAngles.h"

#include "TTree.h"
#include "TLorentzVector.h"

#include <string>
#include <functional>
#include <utility>
#include <array>
#include <iostream>

// template<typename T>
// bool alwaysTrue(T) { return true; }


template<typename T, size_t N, typename U>
class Wrapper {
public:
  Wrapper(const std::string& name, const std::string& prefix,
          std::array<std::pair<std::string, std::function<U(const T)>>, N> funcs) :
    _wrapped(T{}), _name(name), _outPrefix(prefix), _conversionFunctions(std::move(funcs))
  {}

  virtual void Init(TTree *t) {
    t->SetBranchAddress(_name.c_str(), &_wrapped);
  }

  virtual void Create(TTree *t) {
    for (size_t i = 0; i < N; ++i) {
      const auto& nameFunc = _conversionFunctions[i];
      const auto storeName = _outPrefix + nameFunc.first;
      t->Branch(storeName.c_str(), &_outValues[i]);
    }
  }

  virtual void Convert() {
    for (size_t i = 0; i < N; ++i) {
      _outValues[i] = (_conversionFunctions[i].second)(_wrapped);
    }
  }

  // virtual bool Filter(std::function<bool(T)> filterFunc=alwaysTrue) const { return filterFunc(_wrapped); }

  typedef T WrappedType;

  const T& wrapped() const { return _wrapped; }

private:
  T _wrapped;
  std::string _name;
  std::string _outPrefix;
  std::array<std::pair<std::string, std::function<U(const T)>>, N> _conversionFunctions;
  std::array<U, N> _outValues{};
};


template<typename T>
class IdentityWrapper : public Wrapper<T, 1, T> {
public:
  IdentityWrapper() = delete;

  IdentityWrapper(const std::string& name, const std::string outname) :
    Wrapper<T, 1, T>(name, outname, { std::make_pair("", [] (const T v) { return v; }) }) {}
};


class TLorentzVectorWrapper : public Wrapper<TLorentzVector*, 9, double> {
public:
  TLorentzVectorWrapper(const std::string name, const std::string prefix) :
    Wrapper<TLorentzVector*, 9, double>(name, prefix, {
        // all information that is necessary to fully reconstruct TLorentzVector plus some more
        std::make_pair("Pt", [] (const TLorentzVector* t) { return t->Pt(); }),
          std::make_pair("Rap", [] (const TLorentzVector* t) { return t->Rapidity(); }),
          std::make_pair("Mass", [] (const TLorentzVector* t) { return t->M(); }),
          std::make_pair("Eta", [] (const TLorentzVector* t) { return t->Eta(); }),
          std::make_pair("Phi", [] (const TLorentzVector* t) { return t->Phi(); }),
          std::make_pair("P", [] (const TLorentzVector* t) { return t->P(); }),
          std::make_pair("Px", [] (const TLorentzVector* t) { return t->Px(); }),
          std::make_pair("Py", [] (const TLorentzVector* t) { return t->Py(); }),
          std::make_pair("Pz", [] (const TLorentzVector* t) { return t->Pz(); })
          }) {}
};


template<typename W1, typename W2, typename U>
class FunctionWrapper {
public:
  FunctionWrapper(const std::string& outname,
                  std::function<U(typename W1::WrappedType, typename W2::WrappedType)> func,
                  const W1& w1, const W2& w2) :
    _outName(outname), _func(std::move(func)), _outValue(U{}), _w1(w1), _w2(w2) {}

  virtual void Create(TTree *t) {
    t->Branch(_outName.c_str(), &_outValue);
  }

  virtual void Convert() {
    _outValue = _func(_w1.wrapped(), _w2.wrapped());
  }

  // virtual bool Filter(std::function<bool(typename W1::WrappedType, typename W2::WrappedType)> filterFunc=alwaysTrue) const { return filterFunc(_w1.wrapped(), _w2.wrapped()); }

  // virtual void Print() {
  //   std::cout << _outName << ": " << _w1.wrapped() << " " << _w2.wrapped() << std::endl;
  // }

protected:
  std::string _outName;
  std::function<U(typename W1::WrappedType, typename W2::WrappedType)> _func;
  U _outValue;
  const W1& _w1;
  const W2& _w2;
};


template<RefFrame Frame>
class CosthPhiWrapper : public FunctionWrapper<TLorentzVectorWrapper, TLorentzVectorWrapper, Angles> {
public:
  CosthPhiWrapper(const std::string& outname, const TLorentzVectorWrapper& muMinus, const TLorentzVectorWrapper& muPlus, const std::string& frame) :
    FunctionWrapper<TLorentzVectorWrapper, TLorentzVectorWrapper, Angles>(outname, [] (const TLorentzVector *muN, const TLorentzVector *muP)
                                                                          { return calcAnglesInFrame(*muN, *muP, Frame); },
                                                                          muMinus, muPlus),
    _frame(frame) {}

  virtual void Create(TTree *t) override
  {
    t->Branch((_outName + "costh" + _frame).c_str(), &_costh);
    t->Branch((_outName + "phi" + _frame).c_str(), &_phi);
    t->Branch((_outName + "cosalpha" + _frame).c_str(), &_cosalpha);
  }

  virtual void Convert() override
  {
    auto angles = _func(_w1.wrapped(), _w2.wrapped());
    _costh = angles.costh;
    _phi = angles.phi;
    _cosalpha = angles.cosalpha;
  }

private:
  double _costh;
  double _phi;
  double _cosalpha;
  std::string _frame;
};

#endif
