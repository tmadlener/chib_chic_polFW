// combined fit of two ratios imposing the same lambda_theta1 and lambda_theta2
// onto all passed pT bins
// compile:
// g++ -std=c++14 -Wall -Wextra $(root-config --cflags --libs) -I../../ -o two_dim_chi2_shape_fit two_dim_chi2_shape_fit.C

#include "general/interface/ArgParser.h"
#include "general/interface/misc_utils.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "HFitInterface.h"
#include "TGraphAsymmErrors.h"
#include "Math/WrappedMultiTF1.h"
#include "TFile.h"
#include "TF1.h"
#include "Math/Minimizer.h"
#include "TMath.h"
#include "TH2D.h"

#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <array>
#include <numeric>

constexpr auto GRAPHNAME = "r_chic2_chic1_v_costh_HX_fold_bin_0";
constexpr auto POSSCHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";


std::string rand_str(int len=12) {
  std::mt19937 generator{std::random_device{}()};
  std::uniform_int_distribution<int> distribution(0, 51); // Assuming that there are 52 elements in POSSCHARS

  std::string str(len, '\0');
  for (auto& c: str) {
    c = POSSCHARS[distribution(generator)];
  }

  return str;
}

TF1 Rcosth(const std::string& name) {
  auto func = TF1(name.c_str(),
                  // "[0] * (1 + ([1] + [2]) * x[0]*x[0]) / (1 + [1] * x[0]*x[0])",
                  "[0] * (1 + [2] * x[0]*x[0]) / (1 + [1] * x[0]*x[0])",
                  0, 1);
  func.SetParameters(0.5, 0, 0);
  return func;
}

class CombinedChi2 {
public:
  CombinedChi2() = delete;
  CombinedChi2(const std::vector<TGraphAsymmErrors*>& ratios, std::vector<TF1>& functions) {

    // Set the data options (these are the same for all bins)
    ROOT::Fit::DataOptions options;
    options.fIntegral = false; // evaluate the function at the coordinate
    options.fUseEmpty = false; // Ignore empty bins
    options.fCoordErrors = false; // Do not use the errors in the x-directions
    options.fAsymErrors = true; // Use asymmetric errors if present

    m_wrappedFuncs.reserve(ratios.size());
    m_binData.reserve(ratios.size());
    m_functions.reserve(ratios.size());

    for (size_t i = 0; i < ratios.size(); ++i) {
      m_wrappedFuncs.push_back(ROOT::Math::WrappedMultiTF1(functions[i], functions[i].GetNdim()));
      auto& wrappedF1 = m_wrappedFuncs.back();

      m_binData.emplace_back(options);
      auto& binData = m_binData.back();
      ROOT::Fit::FillData(binData, ratios[i]);

      m_functions.push_back(new ROOT::Fit::Chi2Function(binData, wrappedF1));
    }
  }

  double operator() (const double* par) const {
    double sum_chi2 = 0;
    for (size_t i = 0; i < m_functions.size(); ++i) {
      double params[3] = {par[i+2], par[0], par[1]};
      // std::cout << "p = {" << params[0] << "," << params[1] << "," << params[2]
      //           << "} -> chi2 = " << (*m_functions[i])(params) << "\n";
      sum_chi2 += (*m_functions[i])(params);
    }

    // std::cout << sum_chi2 << "\n";
    return sum_chi2;
  }

  int Size() const {
    int size = 0;
    for (const auto & data : m_binData) size += data.Size();
    return size;
  }

  int nPars() const {
    // n normalizations plus two shape parameters
    return m_binData.size() + 2;
  }

private:
  std::vector<ROOT::Fit::BinData> m_binData;
  std::vector<ROOT::Math::WrappedMultiTF1> m_wrappedFuncs;
  std::vector<ROOT::Math::IMultiGenFunction*> m_functions;
};



struct CombinedFitter {
  CombinedFitter() {
    ROOT::Math::MinimizerOptions minOpt;
    minOpt.SetMinimizerType("Minuit2");
    // minOpt.SetPrintLevel(1);

    fitter.Config().SetMinimizerOptions(minOpt);
  }

  bool Fit(CombinedChi2& chi2) {
    const int nPar = chi2.nPars();
    // initialize everything to the expected normalization in the beginning
    std::vector<double> startPar(nPar, 0.45);

    const auto ret = fitter.FitFCN(nPar, chi2, startPar.data(), chi2.Size(), true);

    fitter.GetMinimizer()->PrintResults();
    return ret;
  }

  TGraph get2DContour(double confLevel=0.683, unsigned nPoints=50) {
    auto *minimizer = fitter.GetMinimizer();
    if (!minimizer) {
      std::cout << "Could not get minimizer from fitter\n";
    }

    std::cout << "Recalculating the errors for the desired confLevel = " << confLevel << "\n";
    const double oldErrDef = minimizer->ErrorDef();
    minimizer->SetErrorDef(oldErrDef * TMath::ChisquareQuantile(confLevel, 2));
    const bool hesse = fitter.CalculateHessErrors();
    minimizer->PrintResults();
    const bool minos = fitter.CalculateMinosErrors();
    minimizer->PrintResults();

    if (!hesse || !minos) {
      std::cout << "HESSE (" << hesse << ") or MINOS (" << minos << ") errors could not be determined\n";
      return TGraph();
    }

    std::vector<double> x_c(50, 0.0);
    std::vector<double> y_c(50, 0.0);

    if (!minimizer->Contour(0, 1, nPoints, x_c.data(), y_c.data())) {
      std::cout << "Could not get contour\n";
      return TGraph();
    }

    minimizer->SetErrorDef(oldErrDef); // reset

    return TGraph(nPoints, x_c.data(), y_c.data());
  }

  TGraph getCentralValues() const {
    const auto vals = fitter.Result().Parameters();
    double x[1] = {vals[0]};
    double y[1] = {vals[1]};

    return TGraph(1, x, y);
  }


  TH2D getChi2Scan(const std::vector<double>& xBinning, const std::vector<double>& yBinning,
                   CombinedChi2& chi2) {
    TH2D scanhist("h_chi2scan", "",
                  xBinning.size() - 1, xBinning.data(),
                  yBinning.size() - 1, yBinning.data());

    if (!fitter.Result().IsValid()) {
      std::cout << "Fit result is not valid\n";
      return scanhist;
    }

    std::vector<double> params = fitter.Result().Parameters();
    const double minChi2 = fitter.Result().Chi2();

    for (size_t i = 0; i < xBinning.size() - 1; ++i) {
      params[0] = 0.5 * (xBinning[i] + xBinning[i + 1]);
      for (size_t j = 0; j < yBinning.size() - 1; ++j) {
        params[1] = 0.5 * (yBinning[j] + yBinning[j + 1]);

        // when filling account for under and overflow bins
        scanhist.SetBinContent(i + 1, j + 1, chi2(params.data()) - minChi2);
      }
    }

    return scanhist;
  }

  ROOT::Fit::Fitter fitter{};
};


#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  const auto parser = ArgParser(argc, argv);
  const auto infiles = parser.getOptionVal<std::vector<std::string>>("--infiles");
  const auto outfile = parser.getOptionVal<std::string>("--outfile", "fit_results.root");


  // open files, load graphs, setup functions
  std::vector<TFile*> files;
  std::vector<TGraphAsymmErrors*> graphs;
  std::vector<TF1> funcs;
  for (const auto& name: infiles) {
    files.push_back(TFile::Open(name.c_str()));
    auto *file = files.back();
    graphs.push_back(static_cast<TGraphAsymmErrors*>(file->Get(GRAPHNAME)));
    funcs.push_back(Rcosth(rand_str()));
  }

  CombinedChi2 globalChi2(graphs, funcs);

  std::cout << globalChi2.nPars() << "\n";

  CombinedFitter fitter;
  fitter.Fit(globalChi2);

  auto ofile = new TFile(outfile.c_str(), "recreate");
  ofile->cd();

  auto central = fitter.getCentralValues();
  central.SetName("min_lambda1_lambda2");
  central.Write();

  auto oneSigmaContour = fitter.get2DContour(0.683, 35);
  oneSigmaContour.SetName("cont_2d_lambda1_lambda2_0p683");
  if (oneSigmaContour.GetN() > 0) {
    oneSigmaContour.Write();
  }

  // auto twoSigmaContour = fitter.get2DContour(0.954, 35);
  // twoSigmaContour.SetName("cont_2d_lambda1_lambda2_0p95");
  // if (twoSigmaContour.GetN() > 0) {
  //   twoSigmaContour.Write();
  // }

  // auto threeSigmaContour = fitter.get2DContour(0.997, 35);
  // threeSigmaContour.SetName("cont_2d_lambda1_lambda2_0p99");
  // if (threeSigmaContour.GetN() > 0) {
  //   threeSigmaContour.Write();
  // }

  // auto scanHist = fitter.getChi2Scan(linspace(-5.0, 16.0, 500), linspace(-4.0, 20.0, 500), globalChi2);
  // scanHist.Write();

  ofile->Write("", TObject::kWriteDelete);
  ofile->Close();

  return 0;
}
#endif
