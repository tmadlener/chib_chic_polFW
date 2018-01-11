#ifndef CHIBCHICPOLFW_UNBINNEDPOLFW_BKGWEIGHTCALCCHIC_H__
#define CHIBCHICPOLFW_UNBINNEDPOLFW_BKGWEIGHTCALCCHIC_H__

#include "Region.h"

#include "general/root_utils.h"

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"

#include "TF1.h"
#include "TH1D.h"

#include <utility>
#include <string>
#include <cmath>

///////////////////////////////////////////////////////////////////////////
// SOME CONST VARIABLES THAT SHOULD BE MOVED INTO A CONFIG FILE SOMETIME //
///////////////////////////////////////////////////////////////////////////
constexpr double mChic1Lo = 0.05;
constexpr double mChic1Hi = 0.99;
constexpr double mChic2Lo = 0.05;
constexpr double mChic2Hi = 0.99;
constexpr double mChicBkgLo = 0.05; // in terms of cdf of chic1
constexpr double mChicBkgHi = 0.999; // in terms of cdf of chic2
constexpr double mLSBLo = 3.325; // GeV
constexpr double mRSBHi = 3.725; // GeV
constexpr double nSigmaPR = 3.0; // width of prompt region (sigmas of resolution)
constexpr double ctMax = 1.00625; // maximum value of ctau
const std::string mBkgName = "M_background";
const std::string ltBkgName = "L_background";
const std::string mVarName = "chicMass";
const std::string ltVarName = "Jpsict";
const std::string mlBkgName = "ML_background";
const std::string mlCombBkgName = "ML_comb_background";
///////////////////////////////////////////////////////////////////////////

struct MassRegions {
  Region<Boundary::TwoSided> LSB;
  Region<Boundary::TwoSided> RSB;
  Region<Boundary::TwoSided> SR1;
  Region<Boundary::TwoSided> SR2;
};

MassRegions calcMassRegions(RooWorkspace *ws, const std::string& snapName = "",
                            const std::string& namePrefix = "")
{
  // using namespace RooFit;
  auto *mass = getVar(ws, "chicMass");
  if (!snapName.empty()) {
    ws->loadSnapshot(snapName.c_str());
  }

  // get the cdf from the pdf (stored in the workspace) as TF1 to be able to easily evaluate it
  // at arbitray x-values (haven't found a way to do this in RooFit directly)
  const TF1 *cdf1 = ws->pdf("M_chic1")->createCdf(RooArgSet(*mass))->asTF(RooArgSet(*mass));
  const TF1 *cdf2 = ws->pdf("M_chic2")->createCdf(RooArgSet(*mass))->asTF(RooArgSet(*mass));

  const double lowChic1 = cdf1->GetX(mChic1Lo);
  const double highChic1 = cdf1->GetX(mChic1Hi);
  double lowChic2 = cdf2->GetX(mChic2Lo);
  const double highChic2 = cdf2->GetX(mChic2Hi);

  if (lowChic2 < highChic1) lowChic2 = highChic1;

  const double lsbHigh = cdf1->GetX(mChicBkgLo);
  const double rsbLow  = cdf2->GetX(mChicBkgHi);

  return MassRegions{Region<Boundary::TwoSided>(mLSBLo, lsbHigh, namePrefix + "LSB"),
      Region<Boundary::TwoSided>(rsbLow, mRSBHi, namePrefix + "RSB"),
      Region<Boundary::TwoSided>(lowChic1, highChic1, namePrefix + "SR1"),
      Region<Boundary::TwoSided>(lowChic2, highChic2, namePrefix + "SR2")};
}

double getIntegralInRegion(RooWorkspace *ws, const Region<Boundary::TwoSided>& region,
                           const std::string& pdfname, const std::string& varname,
                           const std::string& snapName = "")
{
  using namespace RooFit;
  if (!snapName.empty()) {
    ws->loadSnapshot(snapName.c_str());
  }

  auto *var = getVar(ws, varname);
  var->setRange(region.name().c_str(), region.min(), region.max());

  auto *pdf = ws->pdf(pdfname.c_str());
  double rangeInt = pdf->createIntegral(RooArgSet(*var), NormSet(RooArgSet(*var)),
                                        Range(region.name().c_str()))->getVal();

  return rangeInt;
}

std::pair<double, double> calcMassWeights(RooWorkspace* ws, const MassRegions& mr, const std::string& pdfname,
                                          const std::string& varname, const std::string& snapname)
{
  const auto intLSB = getIntegralInRegion(ws, mr.LSB, pdfname, varname, snapname);
  const auto intRSB = getIntegralInRegion(ws, mr.RSB, pdfname, varname, snapname);
  const auto intSR1 = getIntegralInRegion(ws, mr.SR1, pdfname, varname, snapname);
  const auto intSR2 = getIntegralInRegion(ws, mr.SR2, pdfname, varname, snapname);

  const double wSR1 = - intSR1 / (intLSB + intRSB);
  const double wSR2 = - intSR2 / (intLSB + intRSB);

  return {wSR1, wSR2};
}

double getIntegralInRegion(RooWorkspace *ws, const std::vector<Region<Boundary::TwoSided>> &regions,
                           const std::string &pdfname, const std::string &varname)
{
  // COULDDO: make this the main method, since RooFit in principle allows to have more
  // than one regions in the integral statement
  double regionSum = 0;
  for (const auto &r : regions) {
    regionSum += getIntegralInRegion(ws, r, pdfname, varname);
  }
  return regionSum;
}


struct LifeTimeRegions {
  Region<Boundary::TwoSided> PR;
  Region<Boundary::TwoSided> NP;
};

LifeTimeRegions calcLifetimeRegions(RooWorkspace *ws, const Region<Boundary::TwoSided> &massSR,
                                    const std::string &dataname)
{
  std::cout << "Calculating lifetime regions\n";
  using namespace RooFit;

  auto *ct = getVar(ws, "Jpsict");
  auto *ctErr = getVar(ws, "JpsictErr");

  auto* promptPdf = static_cast<RooAbsPdf*>(ws->pdf("L_TotalPromptLifetime"));

  // NOTE: this is currently hardcoded for running the first tests
  // TODO: either detect this automatically via some iterator on th workspace
  //       or pass rap/pt bin into this function
  auto *data = static_cast<RooDataSet*>(ws->data(dataname.c_str()));

  auto *dataSR = data->reduce(massSR.getCutStr("chicMass").c_str());
  auto *dataJpsictErr = static_cast<RooDataSet*>(dataSR->reduce(SelectVars(RooArgSet(*ctErr)), Name("dataJpsictErr")));
  RooDataSet *PromptPseudoData = promptPdf->generate(*ct, ProtoData(*dataJpsictErr));

  constexpr int nbinsSigmaDef = 200;
  TH1D* hist1D = static_cast<TH1D*>(PromptPseudoData->createHistogram("hist1D", *ct,
                                                                      Binning(nbinsSigmaDef)));

  const double ctres = hist1D->GetRMS();

  return LifeTimeRegions{Region<Boundary::TwoSided>(-ctres * nSigmaPR, ctres * nSigmaPR, "PR"),
      Region<Boundary::TwoSided>(ctres * nSigmaPR, ctMax, "NP")};
}

LifeTimeRegions calcLifeTimeRegionsJpsi(RooWorkspace* ws, const std::string& binname)
{
  using namespace RooFit;

  auto *ct = getVar(ws, "Jpsict");
  auto *ctErr = getVar(ws, "JpsictErr");

  auto *promptPdf = static_cast<RooAbsPdf*>(ws->pdf("jpsi_TotalPromptLifetime"));

  // auto *data = static_cast<RooDataSet*>(ws->data(binname.c_str()));
  auto *data = ws->data(binname.c_str());
  auto *dataCtErr = static_cast<RooDataSet*>(data->reduce(SelectVars(RooArgSet(*ctErr)),
                                                          Name("dataJpsictErr")));

  RooDataSet *promptPseudoData = promptPdf->generate(*ct, ProtoData(*dataCtErr));

  constexpr int nBinsSigDef = 200; // bins to use for the histogram to define the regions
  TH1D *hist = static_cast<TH1D*>(promptPseudoData->createHistogram("hist", *ct, Binning(nBinsSigDef)));

  const double ctres = hist->GetRMS();
  return LifeTimeRegions{Region<Boundary::TwoSided>(-ctres * nSigmaPR, ctres * nSigmaPR, "PR"),
      Region<Boundary::TwoSided>(ctres * nSigmaPR, ctMax, "NP")};
}

struct BkgWeights {
  double wMBkg1;
  double wMBkg2;
  double wNP;
};

BkgWeights calculateRegionWeights(RooWorkspace *ws, const MassRegions &mr,
                                  const LifeTimeRegions &ltr)
{
  const double intMBkg = getIntegralInRegion(ws, {mr.LSB, mr.RSB}, mBkgName, mVarName);
  const double intMSR1 = getIntegralInRegion(ws, mr.SR1, mBkgName, mVarName);
  const double intMSR2 = getIntegralInRegion(ws, mr.SR2, mBkgName, mVarName);

  const double intNP = getIntegralInRegion(ws, ltr.NP, ltBkgName, ltVarName);
  const double intPR = getIntegralInRegion(ws, ltr.PR, ltBkgName, ltVarName);

  const double fMBkg1 = (intMBkg + intMSR1) / intMBkg;
  const double fMBkg2 = (intMBkg + intMSR2) / intMBkg;
  const double fNP = (intNP + intPR) / intNP;

  std::cout << intNP << " " << intPR << " -> " << fNP << " " << 1 - fNP << "\n";

  return BkgWeights{std::abs(1 - fMBkg1), std::abs(1 - fMBkg2), std::abs(1 - fNP)};
}

struct BkgWeights2D {
  double wChic1;
  double wChic2;
};

/**
 * Integrate pdf in 2D using the passed mass and lifetime region as boundaries.
 */
double get2DIntegralInRegions(RooWorkspace *ws, const Region<Boundary::TwoSided> &mR,
                              const Region<Boundary::TwoSided> &lR, const std::string &pdfname,
                              const std::string &mName, const std::string &lName,
                              const std::string &snapName = "")
{
  using namespace RooFit;
  if (!snapName.empty()) {
    ws->loadSnapshot(snapName.c_str());
  }

  const std::string rangeName = mR.name() + "_" + lR.name();

  auto *mVar = getVar(ws, mName);
  mVar->setRange(rangeName.c_str(), mR.min(), mR.max());
  auto *lVar = getVar(ws, lName);
  lVar->setRange(rangeName.c_str(), lR.min(), lR.max());
  RooArgSet mlSet(*mVar, *lVar);

  auto *pdf = ws->pdf(pdfname.c_str());

  return pdf->createIntegral(mlSet, NormSet(mlSet), Range(rangeName.c_str()))->getVal();
}


/**
 * calculate bkg weight 1 - N_all / (N_bkg) by integrating the 2D mass-lifetime
 * distribution in appropriate ranges.
 */
BkgWeights2D calculate2DWeights(RooWorkspace *ws, const MassRegions &mr,
                                const LifeTimeRegions &ltr, const std::string &snapname = "")
{
  // calculate the total background by integrating the LSB, RSB and signal regions
  // over the full lifetime region.
  // To get to the weight w = 1 - N_bkg,all / N_bkg,SBs, N_bkg,SBs is calculated as
  // N_bkg,all - N_bkg,SR
  // as background use weighted sum of ML_background and ML_comb_background
  const Region<Boundary::TwoSided> fullLT(ltr.PR.min(), ltr.NP.max(), "fullL");

  const double combWeight = getVarVal(ws, "jpsi_fBkg") / getVarVal(ws, "fracBackground");

  const double fullLSB = get2DIntegralInRegions(ws, mr.LSB, fullLT, mlCombBkgName,
                                                mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.LSB, fullLT, mlBkgName, mVarName, ltVarName, snapname);

  const double fullRSB = get2DIntegralInRegions(ws, mr.RSB, fullLT, mlCombBkgName,
                                                mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.RSB, fullLT, mlBkgName, mVarName, ltVarName, snapname);

  const double fullSR1 = get2DIntegralInRegions(ws, mr.SR1, fullLT, mlCombBkgName,
                                                mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.SR1, fullLT, mlBkgName, mVarName, ltVarName, snapname);

  const double fullSR2 = get2DIntegralInRegions(ws, mr.SR2, fullLT, mlCombBkgName,
                                                mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.SR2, fullLT, mlBkgName, mVarName, ltVarName, snapname);

  // get the background in the signal region(s)
  const double intPRSR1 = get2DIntegralInRegions(ws, mr.SR1, ltr.PR, mlCombBkgName,
                                                 mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.SR1, ltr.PR, mlBkgName, mVarName, ltVarName, snapname);

  const double intPRSR2 = get2DIntegralInRegions(ws, mr.SR2, ltr.PR, mlCombBkgName,
                                                 mVarName, ltVarName, snapname) * combWeight
    + get2DIntegralInRegions(ws, mr.SR2, ltr.PR, mlBkgName, mVarName, ltVarName, snapname);

  // N_bkg = N_all - N_SR
  // 1 - N_all / (N_bkg) = - N_SR / N_bkg
  const double wChic1 = - intPRSR1 / (fullLSB + fullRSB + fullSR1 - intPRSR1);
  const double wChic2 = - intPRSR2 / (fullLSB + fullRSB + fullSR2 - intPRSR2);

  return BkgWeights2D{wChic1, wChic2};
}

double calcLTWeight(RooWorkspace* ws, const LifeTimeRegions& ltr, const std::string& snapname = "")
{
  const double combWeight = getVarVal(ws, "jpsi_fBkg") / getVarVal(ws, "fracBackground");

  const double NP = getIntegralInRegion(ws, ltr.NP, mlBkgName, ltVarName, snapname)
    + combWeight * getIntegralInRegion(ws, ltr.NP, mlCombBkgName, ltVarName, snapname);

  const double PR = getIntegralInRegion(ws, ltr.PR, mlBkgName, ltVarName, snapname)
    + combWeight * getIntegralInRegion(ws, ltr.PR, mlCombBkgName, ltVarName, snapname);

  return - PR / NP;
}

#endif
