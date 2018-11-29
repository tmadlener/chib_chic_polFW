#include "FitAnalyser.h"

#include "TFile.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TFormula.h"
#include "RooStats/SPlot.h"

#include <sstream>
#include <memory>
#include <iomanip>
#include <regex>
#include <cmath>


FitAnalyser::FitAnalyser(const Fitter &f) :
  m_filename(f.m_ofname)
{
  SetWorkspaceName(f.m_wsname);
  SetModelName(f.m_modelname);
  SetFitVariableName(f.m_fitvarname);
}

FitAnalyser::FitAnalyser(csr filename, csr model_name, csr fit_variable_name, csr workspace_name, csr dataset_name, csr snapshot_name, csr fitresult_name) :
  m_filename(filename)
{
  // Check if file is ok
  if (!get_file()) {
    std::cout << "EXITING PROGRAM: [FitAnalyser] could not open file. \n";
    exit(EXIT_FAILURE);
  }

  SetModelName(model_name);
  SetFitVariableName(fit_variable_name);
  SetWorkspaceName(workspace_name);
  if (!dataset_name.empty()) SetDatasetName(dataset_name);
  if (!snapshot_name.empty()) SetSnapshotName(snapshot_name);
  if (!fitresult_name.empty()) SetFitResultName(fitresult_name);
}

void FitAnalyser::SetWorkspaceName(csr workspace_name)
{
  if (workspace_name.empty()) {
    // If there is only one workspace in the file take that one
    auto f = get_file();
    auto keys = f->GetListOfKeys();
    std::vector<RooWorkspace*> wspaces;
    for (auto k : *keys) {
      auto ws = dynamic_cast<RooWorkspace*>(f->Get(k->GetName()));
      if (ws) wspaces.emplace_back(ws);
    }
    if (wspaces.size() == 1) m_wsname = wspaces[0]->GetName();
    else if (wspaces.empty()) {
      std::cout << "No workspace found in file '" << m_filename << "'." << std::endl;
      return;
    }
    else if (wspaces.size() > 1) {
      std::cout << "Too many workspaces found in file '" << m_filename << "':\n";
      for (auto ws : wspaces) std::cout << '\t' << ws->GetName() << "; " << ws->GetTitle() << '\n';
      std::cout << "You have to provide a workspace name." << std::endl;
    }
  }
  else m_wsname = workspace_name;
}

RooWorkspace * FitAnalyser::GetWorkspace(bool lsnap)
{
  if (!m_workspace) {
    m_workspace = dynamic_cast<RooWorkspace*>(get_file()->Get(m_wsname.c_str()));
    if (!m_workspace) std::cout << "Workspace '" << m_wsname << "' not found in file '" << m_filename << "'." << std::endl;

  }

  if (m_workspace && lsnap && !m_snapshotname.empty()) m_workspace->loadSnapshot(m_snapshotname.c_str());

  return m_workspace;
}

RooDataSet * FitAnalyser::GetDataset(csr dataset_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  auto dsname = dataset_name.empty() ? m_datasetname : dataset_name;
  auto ds = dynamic_cast<RooDataSet*>(ws->data(dsname.c_str()));
  if (!ds) std::cout << "No RooDataSet named '" << dsname << "' in RooWorkspace '" << m_wsname << "'." << std::endl;

  return ds;
}

RooRealVar * FitAnalyser::GetVariable(csr variable_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  auto varname = variable_name.empty() ? m_fitvarname : variable_name;
  auto var = ws->var(varname.c_str());
  if (!var) std::cout << "No RooRealVar named '" << varname << "' in RooWorkspace '" << m_wsname << "'." << std::endl;

  return var;
}

RooFitResult * FitAnalyser::GetFitResult(csr fitresult_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  if (!fitresult_name.empty()) {
    auto r = dynamic_cast<RooFitResult*>(ws->genobj(fitresult_name.c_str()));
    if (!r) std::cout << "Could not find RooFitResult '" << fitresult_name << "' in RooWorkspace '" << m_wsname << "'." << std::endl;
    return r;
  }
  if (!m_fitresultname.empty()) {
    auto r = dynamic_cast<RooFitResult*>(ws->genobj(m_fitresultname.c_str()));
    if (!r) std::cout << "Could not find RooFitResult '" << m_fitresultname << "' in RooWorkspace '" << m_wsname << "'." << std::endl;
    return r;
  }
  // else return fitresult if only one exists in workspace
  auto objs = ws->allGenericObjects();
  std::vector<RooFitResult*> allfitresults;
  for (auto o : objs) {
    auto r = dynamic_cast<RooFitResult*>(o);
    if (r) allfitresults.emplace_back(r);
  }
  if (allfitresults.size() == 1) return allfitresults[0];
  if (allfitresults.empty()) std::cout << "No RooFitResult found in RooWorkspace '" << m_wsname << "'." << std::endl;
  if (allfitresults.size() > 1) {
    std::cout << "Too many RooFitResults found in RooWorkspace '" << m_wsname << "':\n";
    for (auto res : allfitresults) std::cout << '\t' << res->GetName() << "; " << res->GetTitle() << '\n';
    std::cout << "You have to provide a RooFitResult name." << std::endl;
  }

  return nullptr;
}

RooAbsPdf * FitAnalyser::GetPdf(csr pdf_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  return ws->pdf(pdf_name.empty() ? m_modelname.c_str() : pdf_name.c_str());
}

double FitAnalyser::GetVariableValue(csr variable_name, bool &ok)
{
  ok = false;

  auto var = GetVariable(variable_name);
  if (!var) return 0;

  return ok = true, var->getValV();
}

void FitAnalyser::PlotFitResult(csr output_file, const std::vector<double> & lines, bool addResultBox, bool noTitle)
{
  if (!output_file.empty()) gROOT->SetBatch(kTRUE);
  auto canv = new TCanvas("fit_results", "Fit Results", 1200, 1000);

  auto var = GetVariable();
  auto ds = GetDataset();
  auto model = dynamic_cast<RooAddPdf*>(GetPdf());
  auto ws = GetWorkspace(true);
  if (!ds || !var || !model || !ws) return;

  auto plot = var->frame();
  auto pull = var->frame();
  std::string titel = noTitle ? "" : "Fit results from " + m_wsname;
  plot->SetTitle((titel + ";" + var->GetTitle() + ";").c_str());
  pull->SetTitle(";;Pull");
  plot->SetMinimum(0);

  // Legend - TODO: find a solution
  //https://root-forum.cern.ch/t/pdf-color-in-tlegend/26621

  // Print variables on the plot
  auto ibox = new TPaveText(0.7, 0.5, 0.99, 0.95, "BRNDC");
  ibox->SetTextSize(ibox->GetTextSize()*0.6);
  {
    auto chi2_ndof = get_chi2_ndof();
    std::stringstream ss;
    ss << std::fixed << "#chi^{2}/N_{dof}: "
      << std::setprecision(2) << chi2_ndof.first << "/"
      << std::setprecision(0) << chi2_ndof.second;
    ibox->AddText(ss.str().c_str());
  }
  ibox->AddText("");

  {
    auto vars = ws->allVars();
    auto it = vars.createIterator();
    TObject *var;
    while ((var = it->Next())) {
      auto v = dynamic_cast<RooRealVar*>(var);
      if (!v) continue;

      std::stringstream ss;
      ss << std::fixed << std::setprecision(4);
      auto val = v->getVal();
      ss << v->GetName() << ": ";
      if (!v->isConstant()) {
        auto min = v->getMin(), max = v->getMax();
        bool val_is_min = (min == val), val_is_max = (max == val);

        if (val_is_min || val_is_max) ss << "#color[2]{ " << val << "}";
        else ss << val;
        ss << ", [ ";
        if (val_is_min) ss << "#color[2]{" << min << "}";
        else ss << min;
        ss << ", ";
        if (val_is_max) ss << "#color[2]{" << max << " }";
        else ss << max;
        ss << "]";
      }
      else {
        ss << val;
      }

      ibox->AddText(ss.str().c_str());
    }
  }
  ibox->Print();



  // Pads for the two plots
  double pull_size = 0.15;
  double big_marg = 0.15;
  double small_marg = 0.05;
  TPad *plot_pad = new TPad("plot", "", 0, pull_size, 1, 1);
  TPad *pull_pad = new TPad("pull", "", 0, 0, 1, pull_size);
  plot_pad->SetMargin(big_marg, small_marg, 0.1, noTitle ? small_marg : 0.1);
  pull_pad->SetMargin(big_marg, small_marg, small_marg, small_marg);
  plot_pad->SetGrid();
  pull_pad->SetGrid();
  pull_pad->SetTicks(1, 1);


  //auto fr = GetFitResult();
  //if (fr) model->plotOn(plot, RooFit::LineColor(kRed), RooFit::VisualizeError(*fr));

  // Plot model and dataset and get pull
  ds->plotOn(plot);
  model->plotOn(plot, RooFit::LineColor(kRed));
  auto pull_histo = plot->pullHist();

  auto y = pull->GetYaxis();
  y->SetLabelSize(big_marg*0.8);
  y->SetTitleOffset(0.3);
  y->SetTitleSize(0.2);
  y->SetNdivisions(6, 2, 0, true);
  y->SetLabelSize(0.2);
  y->SetTickSize(0.01);
  auto x = pull->GetXaxis();
  x->SetLabelSize(0);
  x->SetTickSize(0.15);
  pull->addPlotable(pull_histo, "P");
  pull->addObject(new TLine(var->getMin(), 0, var->getMax(), 0));


  // Plot single pdfs
  auto pdfs = model->pdfList();
  for (auto i = 0, s = pdfs.getSize(); i < s; ++i) {
    auto p = model->plotOn(plot, RooFit::Components(*pdfs.at(i)), RooFit::LineColor(kBlue + 2 - 2 * i), RooFit::LineWidth(2), RooFit::LineStyle(kDashed));
  }

  //plot dataset again, to be on top
  ds->plotOn(plot);

  for (double x : lines) plot->addObject(new TLine(x, plot->GetMinimum(), x, plot->GetMaximum()));

  if (addResultBox) plot->addObject(ibox);

  // Draw to canvas
  canv->cd();
  pull_pad->Draw();
  plot_pad->Draw();
  plot_pad->cd();
  plot->Draw();
  pull_pad->cd();
  pull->Draw();

  if (!output_file.empty()) canv->SaveAs(output_file.c_str());
}

void FitAnalyser::CustomPlot1(csr output_file, int year, bool is_dimuonfit, int data_bins, csr frame)
{
  // 2016:HLT_Dimuon8_Upsilon_Barrel, 2017:HLT_Dimuon10_Upsilon_Barrel_Seagulls
  double lumi = (year == 2016 ? 34.88 : 37.15);
  TString lumiText = (year == 2016 ? "34.88 fb^{-1} (13 TeV)" : "37.15 fb^{-1} (13 TeV)");

  double width = 1200;
  double height = 1000;

  std::vector < std::pair< std::string, std::string > > pdfs{ {"background", "Background" } };

  if (is_dimuonfit) {
    pdfs.push_back({ "ups1s","#Upsilon(1S)" });
    pdfs.push_back({ "ups2s","#Upsilon(2S)" });
    pdfs.push_back({ "ups3s","#Upsilon(3S)" });
  }
  else {
    pdfs.push_back({ "chib1","#chi_{b1}(1P)" });
    pdfs.push_back({ "chib2","#chi_{b2}(1P)" });
  }

  if (!output_file.empty()) gROOT->SetBatch(kTRUE);
  auto canv = new TCanvas("fit_results", "Fit Results", width, height);

  auto var = GetVariable();
  auto ds = GetDataset();
  auto model = dynamic_cast<RooAddPdf*>(GetPdf());
  auto ws = GetWorkspace(true);
  if (!ds || !var || !model || !ws) return;


  // Settings
  // (in NDC, i.e. in percent of Canvas)

  // Pdf colors
  auto model_color = kBlue;
  std::vector<int> colors = {
    kBlack, // background
    kGreen, // peak1
    kRed,   // peak2
    kViolet // peak3
  };
  std::vector<int> lines = {
    kDashed,     // background
    kDashDotted, // peak1
    kDashDotted, // peak2
    kDashDotted  // peak3
  };

  // Label
  std::string xaxis_label = is_dimuonfit ? "m^{#mu#mu}" : "m^{#chi}";
  xaxis_label += " (GeV)";
  double label_size = 0.03;
  double axis_title_size = 0.03;
  double tick_size = 0.01;

  // Margins 
  double mt = 0.1, mb = 0.03, ml = 0.13, mr = 0.05;
  double mbb = 0.08, mtt = 0.02; // between pull and plot (x2)

  //Pull
  double pull_range = 3.5;
  double pull_pad_height = 0.15;

  double pull_scale_y = 1. / pull_pad_height;
  double plot_scale_y = 1. / (1. - pull_pad_height);
  double scale = height / width;

  std::cout << "pull_scale_y:" << pull_scale_y << "\nplot_scale_y:" << plot_scale_y
    << "\nx-scale:" << scale << std::endl;

  //
  // End settings

  // Plots and Pads
  TPad *plot_pad = new TPad("plot", "", 0, pull_pad_height, 1, 1);
  plot_pad->SetMargin(ml* scale, mr* scale, mbb*plot_scale_y, mt*plot_scale_y);
  TPad *pull_pad = new TPad("pull", "", 0, 0, 1, pull_pad_height);
  pull_pad->SetMargin(ml* scale, mr* scale, mb*pull_scale_y, mtt*pull_scale_y);
  canv->cd();

  auto plot = var->frame();
  auto pull = var->frame();

  plot->SetTitle((";" + xaxis_label + ";Events").c_str());
  plot->SetTitleSize(axis_title_size*plot_scale_y);
  plot->SetMinimum(0);

  pull->SetTitle(";;Pull");
  pull_pad->SetTicks(1, 1);
  {
    // Plot
    auto y = plot->GetYaxis();
    y->SetLabelSize(label_size*plot_scale_y);
    y->SetTickSize(tick_size*scale);
    auto x = plot->GetXaxis();
    x->SetTickSize(tick_size*plot_scale_y);
    x->SetTitleSize(axis_title_size*plot_scale_y);
  }
  {
    // Pull
    auto y = pull->GetYaxis();
    y->SetLabelSize(label_size*pull_scale_y);
    y->SetNdivisions(6, 2, 0, true);
    y->SetTitleSize(axis_title_size*pull_scale_y);
    y->SetTickSize(tick_size*scale);
    y->SetLimits(-pull_range, pull_range);
    y->SetTitleOffset(pull_pad_height);
    auto x = pull->GetXaxis();
    x->SetLabelSize(0);
    x->SetTickSize(tick_size*pull_scale_y);
  }


  TLegend *leg = new TLegend(1. - 7.*mr*scale, 1. - 3.5*mt*plot_scale_y, 1. - 1.2*mr*scale, 1. - 1.2*mt*plot_scale_y);
  leg->SetFillColor(0);
  leg->SetLineColor(0);


  // Plot costh range
  if (!is_dimuonfit) {
    // costh var only in costh binned workspaces
    std::string costhvarname = "costh_" + frame;
    auto costh = ws->var(costhvarname.c_str());
    double abs_costh_min = 0, abs_costh_max = 1;
    if (costh) {
      auto low_part = ds->reduce((costhvarname + " <0").c_str());
      auto high_part = ds->reduce((costhvarname + " >0").c_str());
      double min_low, max_low, min_high, max_high, tmp_minlow;
      low_part->getRange(*costh, min_low, max_low);
      high_part->getRange(*costh, min_high, max_high);
      std::cout << "low: " << min_low << "," << max_low << "\n"
        << "high: " << min_high << "," << max_high << "\n";
      tmp_minlow = min_low;
      min_low = abs(max_low);
      max_low = abs(tmp_minlow);
      abs_costh_min = min_low;
      if (min_low > min_high) abs_costh_min = min_high;
      abs_costh_max = max_low;
      if (max_low < max_high) abs_costh_max = max_high;
      abs_costh_max = round(abs_costh_max * 1000) / 1000.;
      abs_costh_min = round(abs_costh_min * 1000) / 1000.;
      if (year == 2016 && std::string(ds->GetName()) == std::string("data_costh_bin_5")) abs_costh_max = 0.86;
      if (year == 2017 && std::string(ds->GetName()) == std::string("data_costh_bin_4")) abs_costh_max = 0.82;

    }
    auto ibox = new TPaveText(0.1 + plot_pad->GetLeftMargin(), 0.7 - plot_pad->GetTopMargin(),
      0.25 + plot_pad->GetLeftMargin(), 0.9 - plot_pad->GetTopMargin(), "NBNDC");
    ibox->SetFillColor(0);
    ibox->SetTextSize(0.042);
    {
      std::stringstream ss;
      ss <</* std::setprecision(2) <<*/ abs_costh_min << " < |cos#vartheta^{" << frame << "}| < " << abs_costh_max;
      ibox->AddText(ss.str().c_str());
      std::cout << ss.str() << std::endl;
    }
    //ibox->AddText("");
    plot->addObject(ibox);
  }

  // Plot model and dataset and get pull
  RooCmdArg data_binning(RooFit::Binning(data_bins, var->getMin(), var->getMax()));
  ds->plotOn(plot, RooFit::Name("dataplot"), data_binning);
  model->plotOn(plot, RooFit::Name("fullmodel"), RooFit::LineColor(model_color));
  auto pull_histo = plot->pullHist();
  leg->AddEntry(plot->findObject("fullmodel"), "Full model", "L");
  leg->AddEntry(plot->findObject("dataplot"), "Data", "PE");

  pull->addPlotable(pull_histo, "P");
  pull->addObject(new TLine(var->getMin(), 0, var->getMax(), 0));


  // Plot single pdfs

  for (int i = 0, s = pdfs.size(); i < s; ++i) {
    auto pdfname = pdfs.at(i).first.c_str();
    auto pdftitle = pdfs.at(i).second.c_str();
    auto pdf = ws->pdf(pdfname);
    if (!pdf) continue;
    auto name = ("pdf" + std::to_string(i)).c_str();
    model->plotOn(plot, RooFit::Name(name), RooFit::Components(*pdf), RooFit::LineColor(colors.at(i)), RooFit::LineStyle(lines.at(i)), RooFit::LineWidth(2));
    leg->AddEntry(plot->findObject(name), pdftitle, "L");
  }

  //plot dataset again, to be on top
  ds->plotOn(plot, data_binning);

  //for (double x : lines) plot->addObject(new TLine(x, plot->GetMinimum(), x, plot->GetMaximum()));

  plot->addObject(leg);

  if (is_dimuonfit) {
    auto f = get_file();
    auto min = (RooRealVar*)f->Get("dimuon_cut_min");
    auto max = (RooRealVar*)f->Get("dimuon_cut_max");
    if (min && max) for (double x : {min->getValV(), max->getValV()}) {
      auto line = new TLine(x, plot->GetMinimum(), x, plot->GetMaximum());
      line->SetLineColor(kRed + 1);
      line->SetLineWidth(2);
      line->Print();
      plot->addObject(line);
      double xndc = (x - 8.7) / 2.4 - 0.1;
      auto box = new TPaveText(xndc + plot_pad->GetLeftMargin(), 0.90 - plot_pad->GetTopMargin(),
        xndc + plot_pad->GetLeftMargin() + 0.1, 0.95 - plot_pad->GetTopMargin(), "NDC");
      std::stringstream ss;
      ss << std::setprecision(3) << x;
      box->AddText(ss.str().c_str());
      box->SetFillColor(0);
      box->SetTextColor(kRed + 1);
      box->Print();
      plot->addObject(box);
    }
  }
  // Draw to canvas
  canv->cd();
  pull_pad->Draw();
  plot_pad->Draw();
  plot_pad->cd();
  plot->Draw();

  pull_pad->cd();
  pull->Draw();

  plot_pad->cd();
  // Add CMS and Lumi
  add_cms_lumi(plot_pad, lumi);

  canv->Update();


  std::cout << "PULL TITLE SIZE:" << pull->GetYaxis()->GetTitleSize()
    << "\nPLOT pad width:" << plot_pad->GetWw() << ", " << plot_pad->GetWNDC()
    << "\nPULL pad width:" << pull_pad->GetWw() << ", " << pull_pad->GetWNDC() << std::endl;


  if (!output_file.empty()) canv->SaveAs(output_file.c_str());
}

void FitAnalyser::CustomPlot(csr output_file)
{
  std::vector<double> lines;
  bool noTitle = true;
  bool addResultBox = true;

  if (!output_file.empty()) gROOT->SetBatch(kTRUE);
  auto canv = new TCanvas("fit_results", "Fit Results", 1200, 1000);

  auto var = GetVariable();
  auto ds = GetDataset();
  auto model = dynamic_cast<RooAddPdf*>(GetPdf());
  auto ws = GetWorkspace(true);
  if (!ds || !var || !model || !ws) return;

  std::string xachsentitel = "m_{#mu^{+}#mu^{-}";
  if (m_fitvarname == "dimuon_mass") xachsentitel += "}[GeV]";
  else xachsentitel += "#gamma}[GeV]";

  var->SetTitle(xachsentitel.c_str());

  auto plot = var->frame();
  auto pull = var->frame();
  std::string titel = noTitle ? "" : "Fit results from " + m_wsname;
  plot->SetTitle((titel + ";" + var->GetTitle() + ";").c_str());
  pull->SetTitle(";;Pull");
  plot->SetMinimum(0);

  // Pads for the two plots
  double pull_size = 0.15;
  double big_marg = 0.15;
  double small_marg = 0.05;
  TPad *plot_pad = new TPad("plot", "", 0, pull_size, 1, 1);
  TPad *pull_pad = new TPad("pull", "", 0, 0, 1, pull_size);
  // Werte beziehen sich immer auf aktuelles Pad, also 0.1 bedeuted etwas anderes im pull_pad als im plot_pad
  plot_pad->SetMargin(big_marg, small_marg, 0.1, noTitle ? small_marg : 0.1);
  pull_pad->SetMargin(big_marg, small_marg, 0.2, 0.1);
  plot_pad->SetGrid();
  pull_pad->SetGrid();
  pull_pad->SetTicks(1, 1);


  // Legend - 
  // This works: https://root-forum.cern.ch/t/roofit-tlegend/6356/5

  // Print variables on the plot
  double iboxy1 = 0.7;
  if (m_fitvarname == "dimuon_mass") iboxy1 = 0.8;

  auto ibox = new TPaveText(0.7, iboxy1, 0.98 - plot_pad->GetTopMargin(), 0.98 - plot_pad->GetRightMargin(), "NBNDC");
  ibox->SetFillColor(0);
  ibox->SetTextSize(ibox->GetTextSize()*0.9);
  ibox->AddText("8 < p_{T}^{#mu#mu} < 50 GeV");
  ibox->AddText("");

  {
    auto chi2_ndof = get_chi2_ndof();
    std::stringstream ss;
    ss << std::fixed << "#chi^{2}/N_{dof}: "
      << std::setprecision(2) << chi2_ndof.first << "/"
      << std::setprecision(0) << chi2_ndof.second;
    ibox->AddText(ss.str().c_str());
  }
  ibox->AddText("");

  {
    auto vars = { ws->var("N1"), ws->var("N2"), ws->var("Nbg") };
    for (TObject* var : vars) {
      auto v = dynamic_cast<RooRealVar*>(var);
      if (!v) continue;

      std::stringstream ss;
      ss << std::fixed << std::setprecision(0);
      ss << v->GetName() << ": " << v->getVal();

      ibox->AddText(ss.str().c_str());
    }
  }
  TPaveText *mintext = nullptr, *maxtext = nullptr;
  if (m_fitvarname == "dimuon_mass") {

    bool ok = false;
    double rangemin = EvaluateFormula("mu1s-3*TMath::Sqrt(tail_left*sigma1s1*sigma1s1+(1-tail_left)*sigma1s2*sigma1s2)", ok);
    double rangemax = EvaluateFormula("mu1s+3*TMath::Sqrt(tail_left*sigma1s1*sigma1s1+(1-tail_left)*sigma1s2*sigma1s2)", ok);

    //std::stringstream ss;
    //ss << std::fixed << std::setprecision(2);
    //ss << "#Upsilon(1S) 3#sigma region: [" << rangemin << ", " << rangemax << "]";
    //ibox->AddText(ss.str().c_str());

    double rangemaxndc = (rangemax - 8.7) / 2.4;
    double rangeminndc = (rangemin - 8.7) / 2.4;
    auto ibox = new TPaveText(rangemaxndc + 0.07, 0.8, rangemaxndc + 0.18, 0.85, "NBNDC");
    ibox->SetFillColor(0);
    ibox->SetTextColor(kRed + 2);
    ibox->AddText("#mu(#Upsilon(1S)) + 3#sigma");
    {
      std::stringstream ss;
      ss << std::setprecision(3) << rangemax;
      ibox->AddText(ss.str().c_str());
    }
    mintext = ibox;

    ibox = new TPaveText(rangeminndc, 0.8, rangeminndc + 0.11, 0.85, "NBNDC");
    ibox->SetFillColor(0);
    ibox->SetTextColor(kRed + 2);
    ibox->AddText("#mu(#Upsilon(1S)) - 3#sigma");
    {
      std::stringstream ss;
      ss << std::setprecision(3) << rangemin;
      ibox->AddText(ss.str().c_str());
    }
    maxtext = ibox;

    lines.push_back(rangemin);
    lines.push_back(rangemax);
  }

  ibox->Print();


  // Plot model and dataset and get pull
  ds->plotOn(plot);
  model->plotOn(plot, RooFit::LineColor(kRed));
  auto pull_histo = plot->pullHist();

  auto y = pull->GetYaxis();
  y->SetLabelSize(big_marg*0.8);
  y->SetTitleOffset(0.3);
  y->SetTitleSize(0.2);
  y->SetLabelSize(0.2);
  y->SetTickSize(0.01);
  y->SetNdivisions(4, 0, 0, false);
  auto x = pull->GetXaxis();
  x->SetLabelSize(0);
  x->SetTickSize(0.15);

  double pullmin = pull_histo->getYAxisMin();
  double pullmax = pull_histo->getYAxisMax();
  double pullrange = ceil(abs(pullmax));
  if (abs(pullmin) > abs(pullmax)) pullrange = ceil(abs(pullmin)) / 1.1;
  pull_histo->setYAxisLimits(-pullrange, pullrange);
  std::cout << "PULLRANGE: " << pullrange << ", " << pullmin << ", " << pullmax << ", " << ceil(abs(pullmin)) << ", " << ceil(abs(pullmax)) << std::endl;

  pull->addPlotable(pull_histo, "P");
  pull->addObject(new TLine(var->getMin(), 0, var->getMax(), 0));

  // Plot single pdfs
  auto pdfs = model->pdfList();
  for (auto i = 0, s = pdfs.getSize(); i < s; ++i) {
    auto p = model->plotOn(plot, RooFit::Components(*pdfs.at(i)), RooFit::LineColor(kBlue + 2 - 2 * i), RooFit::LineWidth(2), RooFit::LineStyle(kDashed));
  }

  //plot dataset again, to be on top
  ds->plotOn(plot);

  for (double x : lines) {
    auto line = new TLine(x, plot->GetMinimum(), x, plot->GetMaximum());
    line->SetLineColor(kRed + 2);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    plot->addObject(line);
    //plot->addObject(new TLine(x, plot->GetMinimum(), x, plot->GetMaximum()));
  }

  if (addResultBox) plot->addObject(ibox);
  if (m_fitvarname == "dimuon_mass") {
    plot->addObject(mintext);
    plot->addObject(maxtext);
  }

  // Draw to canvas
  canv->cd();
  pull_pad->Draw();
  plot_pad->Draw();
  plot_pad->cd();
  plot->Draw();
  pull_pad->cd();
  pull->Draw();

  if (!output_file.empty()) canv->SaveAs(output_file.c_str());
}

double FitAnalyser::EvaluateFormula(csr formulastring, bool &ok)
{
  ok = false;
  auto ws = GetWorkspace(true);
  if (!ws) return 0;
  std::string formulaexpr(formulastring);
  auto params = ws->allVars();
  auto it = params.createIterator();

  while (auto param = dynamic_cast<RooRealVar*>(it->Next())) {
    std::string param_name(param->GetName());
    if (param_name != m_fitvarname) {
      std::regex reg("\\b" + param_name + "\\b");
      std::string param_val(std::to_string(param->getVal()));

      formulaexpr = std::regex_replace(formulaexpr, reg, param_val);
    }
  }

  // Create a TFormula and check if it's valid
  TFormula formula("temp_formula", formulaexpr.c_str());
  if (!formula.IsValid()) {
    std::cout << "Could not evaluate " << formulaexpr << " from workspace " << ws->GetName() << std::endl;
    return 0;
  }

  return ok = true, formula.Eval(0.0);
}

double FitAnalyser::GetXForQuantile(csr pdf_name, double q, bool &ok)
{
  ok = false;
  if (q > 1 || q < 0) {
    std::cout << "GetXForQuantile: The quantile has to be between 0 and 1, but is " << q << ".";
    return 0;
  }
  //TODO: see https://github.com/tmadlener/phys_utils/blob/unbinned-ml-fit/UnbinnedPolFW/interface/BkgWeightCalcChic.h

}

RooDataSet * FitAnalyser::AddSWeights(const std::vector<std::string>& yield_names)
{
  if (yield_names.empty()) return GetDataset();

  auto ws = GetWorkspace(true);
  auto ds = GetDataset();
  auto pdf = GetPdf();

  if (!ws || !ds || !pdf) return GetDataset();

  RooArgList yield_list;
  for (const auto &yld : yield_names) yield_list.add(*GetVariable(yld));
  yield_list.Print();

  //{
  //  // Set all variables, but yields and fitvar constant
  //  auto it = vars->createIterator();
  //  TObject *var;
  //  while ((var = it->Next())) {
  //    auto v = dynamic_cast<RooRealVar*>(var);
  //    if (!v) continue;
  //    std::string v_name(v->GetName());

  //    if (v_name == m_modelname || v_name == m_fitvarname) continue;
  //    bool yieldcheck_passed = true;
  //    for (const auto &y_name : yield_names) if (v_name == y_name) {
  //      yieldcheck_passed = false;
  //      break;
  //    }
  //    if (!yieldcheck_passed) continue;
  //    const_vars.push_back({ v,v->isConstant() });
  //    std::cout << "Setting variable " << v_name << " constant." << std::endl;
  //    v->setConstant();
  //  }
  //}

  RooStats::SPlot sData("sData", "sData", *ds, pdf, yield_list);

  get_file()->cd();
  ws->Write(0, TObject::kWriteDelete);



  return GetDataset();
  //
  // Next steps:
  // TreeMerger: RooDataset as argument and do the following:
  // TODO in calling program:
  // Set defaultstoragetype to ttree
  // Copy this dataset into a new one
  // get tree from dataset and merge it with the input tree
}

TFile * FitAnalyser::get_file()
{
  if (!m_file)  m_file.reset(new TFile(m_filename.c_str(), "update"));
  if (!m_file) std::cout << "FitAnalyser: Could not open file  '" << m_filename << "'." << std::endl;
  if (m_file->IsZombie()) {
    std::cout << "FitAnalyser: File  '" << m_filename << "' is a zombie." << std::endl;
    m_file.reset(nullptr);
  }

  return m_file.get();
}


std::pair<double, double> FitAnalyser::get_chi2_ndof()
{
  auto var = GetVariable();
  auto fr = GetFitResult();
  auto ds = GetDataset();
  auto model = GetPdf();
  if (!fr || !var || !ds || !model) return{ -1,-1 };

  auto frame = var->frame();
  ds->plotOn(frame);
  model->plotOn(frame);

  int fit_parameter_count = (fr->floatParsFinal()).getSize();
  int bin_count = frame->GetNbinsX();
  double reduced_chi2 = frame->chiSquare(fit_parameter_count);  //reduced chi-squared = chi2/ndof

  auto ndof = bin_count - fit_parameter_count;  //number of degrees of freedom
  auto chi2 = reduced_chi2*ndof;

  return{ chi2, ndof };


}

void FitAnalyser::add_cms_lumi(TPad* c, double lumi, std::string cmsExtraText, double TeV) {
  auto t = c->GetTopMargin();
  auto r = c->GetRightMargin();
  auto l = c->GetLeftMargin();

  // Add CMS text
  std::string cmsText = "CMS";
  int cmsTextFont = 61;
  double cmsTextSize = 0.35;
  double cmsTextOffset = 0.1;
  int extraTextFont = 52; // for cmsExtraText
  double extraTextSize = cmsTextSize*0.8;

  std::stringstream ss;
  ss << lumi << " fb^{-1} (" << TeV << " TeV)";
  std::string lumitext = ss.str();
  std::cout << lumitext << std::endl;
  int  lumiTextFont = 42;
  double lumiTextSize = 0.3;
  double  lumiTextOffset = 0.1;

  auto latex = TLatex();
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);

  latex.SetTextFont(lumiTextFont);
  latex.SetTextAlign(31);
  latex.SetTextSize(lumiTextSize * t);
  latex.DrawLatex(1 - r, 1 - t + lumiTextOffset * t, lumitext.c_str());

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(11);
  latex.SetTextSize(cmsTextSize * t);
  auto cms_text_width = latex.DrawLatex(l, 1 - t + cmsTextOffset*t, cmsText.c_str())->GetXsize();

  latex.SetTextFont(extraTextFont);
  latex.SetTextSize(extraTextSize*t);
  latex.DrawLatex(l + cms_text_width, 1 - t + cmsTextOffset*t, (" " + cmsExtraText).c_str());
}