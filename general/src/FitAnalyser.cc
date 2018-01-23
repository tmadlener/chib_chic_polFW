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

#include <sstream>
#include <memory>
#include <iomanip>


FitAnalyser::FitAnalyser(const Fitter &f)
{
  SetWorkspaceName(f.m_wsname);
  SetFileName(f.m_ofname);
  SetFitVariableName(f.m_fitvarname);
  SetModelName(f.m_modelname);
}

RooWorkspace * FitAnalyser::GetWorkspace(bool lsnap)
{
  auto f = open_file();
  auto ws = dynamic_cast<RooWorkspace*>(f->Get(m_wsname.c_str()));
  if (ws && lsnap) ws->loadSnapshot(m_snapshotname.c_str());
  return ws;
}

RooAbsData * FitAnalyser::GetDataset(csr dataset_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  return ws->data(dataset_name.empty() ? m_datasetname.c_str() : dataset_name.c_str());
}

RooRealVar * FitAnalyser::GetVariable(csr variable_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  return ws->var(variable_name.empty() ? m_fitvarname.c_str() : variable_name.c_str());
}

RooFitResult * FitAnalyser::GetFitResult(csr fitresult_name)
{
  auto ws = GetWorkspace(true);
  if (!ws) return nullptr;

  if (!fitresult_name.empty()) return dynamic_cast<RooFitResult*>(ws->genobj(fitresult_name.c_str()));

  auto objs = ws->allGenericObjects();
  for (auto o : objs) {
    auto r = dynamic_cast<RooFitResult*>(o);
    if (r) return r;
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

void FitAnalyser::PlotFitResult(csr output_file)
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
  plot->SetTitle(("Fit results from " + m_wsname + ";" + var->GetTitle() + ";").c_str());
  pull->SetTitle(";;Pull");

  // Legend - TODO: find a solution
  //https://root-forum.cern.ch/t/pdf-color-in-tlegend/26621

  // Print variables on the plot
  auto ibox = new TPaveText(0.7, 0.5, 0.99, 0.95, "BRNDC");
  ibox->SetTextSize(ibox->GetTextSize()*0.6);
  {
    auto chi2_ndof = get_chi2_ndof();
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << "#chi^{2}/N_{dof}: " << chi2_ndof.first << "/" << chi2_ndof.second;
    ibox->AddText(ss.str().c_str());
  }
  ibox->AddText(("Number of Events: " + std::to_string(ds->numEntries())).c_str());
  ibox->AddText("");

  {
    auto vars = ws->allVars();
    auto it = vars.createIterator();
    TObject *var;
    while ((var = it->Next())) {
      auto v = dynamic_cast<RooRealVar*>(var);
      if (!v) continue;

      std::stringstream ss;
      ss << std::fixed << std::setprecision(2);
      ss << v->GetName() << ": " << v->getVal() << ", [ " << v->getMin() << ", " << v->getMax() << "]";

      ibox->AddText(ss.str().c_str());
    }
  }
  ibox->Print();

  plot->addObject(ibox);


  // Pads for the two plots
  double pull_size = 0.15;
  double big_marg = 0.15;
  double small_marg = 0.05;
  TPad *plot_pad = new TPad("plot", "", 0, pull_size, 1, 1);
  TPad *pull_pad = new TPad("pull", "", 0, 0, 1, pull_size);
  plot_pad->SetMargin(big_marg, small_marg, 0.1, 0.1);
  pull_pad->SetMargin(big_marg, small_marg, small_marg, 0);
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

TFile * FitAnalyser::open_file()
{
  if (!m_file || (m_file && std::string(m_file->GetName()) != m_filename)) {
    m_file.reset(TFile::Open(m_filename.c_str()));
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