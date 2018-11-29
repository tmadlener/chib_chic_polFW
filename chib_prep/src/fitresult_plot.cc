#include "FitAnalyser.h"
#include "ArgParser.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"

#include <string>
#include <iostream>

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/LumiPublicResults#2016_proton_proton_13_TeV_collis
std::string lumi_2016 = "37.76 (13 TeV) fb^{-1}"; 
std::string lumi_2017 = "44.98 (13 TeV) fb^{-1}";

int main(int argc, char **argv) {

  ArgParser parser(argc, argv);

  auto filename = parser.getOptionVal < std::string >("--filename", "");
  auto modelname = parser.getOptionVal < std::string >("--modelname","model");
  auto fitvar = parser.getOptionVal < std::string >("--fitvar","chi_mass_rf1S");
  auto wsname = parser.getOptionVal < std::string >("--wsname", "ws-dimuon_mass_8p7_11p1-dimuon_pt_10_50_chib1P1S");
  auto yearstring = parser.getOptionVal < std::string >("--year","");
  auto dsname = parser.getOptionVal < std::string >("--dsname", "");
  auto snapshot = parser.getOptionVal < std::string >("--snapshot", "");
  auto frame = parser.getOptionVal < std::string >("--frame", "HX");
  auto fitresultname = parser.getOptionVal < std::string >("--fitresultname","");
  auto saveas = parser.getOptionVal < std::string >("--saveas", "");
  auto dimuonfit = parser.getOptionVal < bool >("--dimuonfit", false);
  auto databins = parser.getOptionVal < int >("--databins", 100);

  std::cout << "Plotting settings:\n";
  std::cout << "yearstring: " << yearstring << std::endl;
  int year = std::stoi(yearstring);
  std::cout << "filename: " << filename << std::endl;
  std::cout << "year: " << year << std::endl;
  std::cout << "wsname: " << wsname << std::endl;
  std::cout << "modelname: " << modelname << std::endl;
  std::cout << "dsname: " << dsname << std::endl;
  std::cout << "snapshot: " << snapshot << std::endl;
  std::cout << "databins: " << databins << std::endl;
  std::cout << "frame: " << frame << std::endl;

  FitAnalyser ana(filename, modelname, fitvar, wsname, dsname, snapshot, fitresultname);

  ana.CustomPlot1(saveas, year, dimuonfit, databins, frame);



}