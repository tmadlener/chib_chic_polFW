#ifndef CHIBCHICPOLFW_CHICPREP_CHICTUPLEEVENT_H__new // TODO: remove new
#define CHIBCHICPOLFW_CHICPREP_CHICTUPLEEVENT_H__new // TODO: remove new

#include "utils.h"

#include "TTree.h"

#include <vector>
#include <string>
#include <utility>
#include <tuple>

using StringPair = std::pair<std::string, std::string>;
using StringTriple = std::tuple<std::string, std::string, std::string>;
using VStringPair = std::vector<StringPair>;


// to contain everything that is necessary for the calculation of the angles
// i.e. single muons plus associated functionality
struct AngleContainer {
  AngleContainer(const std::string& muPlusName, const std::string& muMinusName, const std::string& outPrefix="") :
    muPlus(TLorentzVectorWrapper(muPlusName, outPrefix + "muP")), muMinus(TLorentzVectorWrapper(muMinusName, outPrefix + "muN")),
    _anglesHX(outPrefix, muMinus, muPlus, "_HX"), _anglesCS(outPrefix, muMinus, muPlus, "_CS"), _anglesPX(outPrefix, muMinus, muPlus, "_PX") {}

  TLorentzVectorWrapper muPlus;
  TLorentzVectorWrapper muMinus;

  CosthPhiWrapper<RefFrame::HX> _anglesHX;
  CosthPhiWrapper<RefFrame::CS> _anglesCS;
  CosthPhiWrapper<RefFrame::PX> _anglesPX;
};

class ChicTupleEvent {
public:
  ChicTupleEvent();

  // ChicTupleEvent(const VStringPair& intVars, const VStringPair& doubleVars, const VStringPair& momentumVars,
  //                const std::vector<StringTriple>& singleMuons);

  void Init(TTree *t);

  void Create(TTree *t);

  void Convert();

  bool Accept();

  // void Print();

  void Add(const VStringPair& intVars, const VStringPair& doubleVars, const VStringPair& momentumVars,
           const std::vector<StringTriple>& singleMuons);

private:
  std::vector<IdentityWrapper<int>> _intVars;
  std::vector<IdentityWrapper<double>> _doubleVars;
  std::vector<TLorentzVectorWrapper> _fourMomentaVars;

  std::vector<AngleContainer> _angles;
};

ChicTupleEvent::ChicTupleEvent()
{
  // reserve enough space
  _intVars.reserve(20);
  _doubleVars.reserve(20);
  _fourMomentaVars.reserve(10);
  _angles.reserve(4);
}

// ChicTupleEvent::ChicTupleEvent(const VStringPair& intVars, const VStringPair& doubleVars, const VStringPair& momentumVars,
//                                const std::vector<StringTriple>& singleMuons)
// {
//   Add(intVars, doubleVars, momentumVars, singleMuons);
// }


void ChicTupleEvent::Init(TTree *t)
{
  for (auto& var : _intVars) {
    var.Init(t);
  }

  for (auto& var : _doubleVars) {
    var.Init(t);
  }

  for (auto& var : _fourMomentaVars) {
    var.Init(t);
  }

  for (auto& muons : _angles) {
    muons.muPlus.Init(t);
    muons.muMinus.Init(t);
  }
}




void ChicTupleEvent::Add(const VStringPair &intVars, const VStringPair &doubleVars, const VStringPair &momentumVars,
                         const std::vector<StringTriple> &singleMuons)
{
  for (const auto& intVar : intVars) {
    _intVars.emplace_back(intVar.first, intVar.second);
  }

  for (const auto& doubleVar : doubleVars) {
    _doubleVars.emplace_back(doubleVar.first, doubleVar.second);
  }

  for (const auto& fourMomentaVar : momentumVars) {
    _fourMomentaVars.emplace_back(fourMomentaVar.first, fourMomentaVar.second);
  }

  for (const auto& singleMuon : singleMuons) {
    _angles.emplace_back(std::get<0>(singleMuon), std::get<1>(singleMuon), std::get<2>(singleMuon));
  }
}

// void ChicTupleEvent::Print()
// {
//   std::cout << "Print() " << "\n";
//   for (auto& singleMuon : _angles) {
//     std::cout << singleMuon.muPlus.wrapped() << " " << singleMuon.muMinus.wrapped() << "\n";
//     singleMuon._anglesHX.Print();
//     singleMuon._anglesPX.Print();
//     singleMuon._anglesCS.Print();
//   }
//   std::cout << "end Print()" << std::endl;
// }

void ChicTupleEvent::Convert()
{
  for (auto& var : _intVars) {
    var.Convert();
  }

  for (auto& var : _doubleVars) {
    var.Convert();
  }

  for (auto& var : _fourMomentaVars) {
    var.Convert();
  }

  for (auto& singleMuon : _angles) {
    singleMuon.muPlus.Convert();
    singleMuon.muMinus.Convert();
    singleMuon._anglesHX.Convert();
    singleMuon._anglesPX.Convert();
    singleMuon._anglesCS.Convert();
  }
}

void ChicTupleEvent::Create(TTree *t)
{
  for (auto& var : _intVars) {
    var.Create(t);
  }

  for (auto& var : _doubleVars) {
    var.Create(t);
  }

  for (auto& var : _fourMomentaVars) {
    var.Create(t);
  }

  for (auto& angle : _angles) {
    angle.muPlus.Create(t);
    angle.muMinus.Create(t);
    angle._anglesHX.Create(t);
    angle._anglesPX.Create(t);
    angle._anglesCS.Create(t);
  }

}


bool ChicTupleEvent::Accept()
{
  // for (const auto& var : _intVars) {
  //   if (!var.Filter()) return false;
  // }

  // for (const auto& var : _doubleVars) {
  //   if (!var.Filter()) return false;
  // }

  // for (const auto& var : _fourMomentaVars) {
  //   if (!var.Filter()) return false;
  // }

  return true;
}


#endif
