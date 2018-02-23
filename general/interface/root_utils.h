#ifndef CHIBCHICPOLFW_GENERAL_ROOTUTILS_H__
#define CHIBCHICPOLFW_GENERAL_ROOTUTILS_H__

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

#include <string>
#include <vector>
#include <iostream>

/** get from Root file. Using dynamic cast to get throught the inheritance stuff correctly. */
template<typename T> inline
T* getFromFile(TFile* f, const std::string& name)
{
  return dynamic_cast<T*>(f->Get(name.c_str()));
}

/** Try to open TFile with passed filename. */
TFile* checkOpenFile(const std::string& filename)
{
  TFile* f = TFile::Open(filename.c_str());
  if (f) return f;

  std::cerr << "Could not open file: \'" << filename << "\'" << std::endl;
  return nullptr;
}

/** Redirect TObject->Print() from std::cout to arbitrary std::ostream (e.g. a file) j.n. 01-2018 */
void print2stream(TObject* object_to_Print, std::ostream & stream)
{
  // From https://stackoverflow.com/questions/10150468/how-to-redirect-cin-and-cout-to-files
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(stream.rdbuf()); //redirect std::cout
  object_to_Print->Print();
  std::cout.rdbuf(coutbuf); //reset to standard output again
}

/** Try to get object with name from TFile f. */
template<typename T>
T* checkGetFromFile(TFile* f, const std::string& name)
{
  T* t = static_cast<T*>(f->Get(name.c_str()));
  if (t) return t;

  std::cerr << "Could not get \'" << name << "\' from TFile \'" << f->GetName() << "\'" << std::endl;
  return nullptr;
}

inline bool checkGetEntry(TTree* t, const int event)
{
  if (t->GetEntry(event) < 0) {
    std::cerr << "I/O error while reading event " << event << " in TTree \'" << t->GetName() << "\'" << std::endl;
    return false;
  }
  return true;
}


RooRealVar* getVar(RooWorkspace* ws, const std::string& name)
{
  auto* var = static_cast<RooRealVar*>(ws->var(name.c_str()));
  if (var) return var;
  var = static_cast<RooRealVar*>(ws->function(name.c_str()));
  if (var) return var;

  std::cerr << "Could not get " << name << " from workspace\n";
  return nullptr;
}

/** get value of variable with name from workspace. */
double getVarVal(RooWorkspace* ws, const std::string& name)
{
  if (auto* var = getVar(ws, name)) {
    return var->getVal();
  }
  return 0; // returning zero to make this bugs a bit more subtle and harder to detect ;)
}

/** get value error of variable with name from workspace. */
double getVarError(RooWorkspace* ws, const std::string& name)
{
  if (auto* var = getVar(ws, name)) {
    return var->getError();
  }
  return 0; // returning zero to make this bugs a bit more subtle and harder to detect ;)
}

/** set workspace variable constant to the passed value. */
inline void setVarConstant(RooWorkspace* ws, const std::string& name, const double val)
{
  getVar(ws, name)->setVal(val);
  getVar(ws, name)->setConstant(true);
}

/** Create TChain from passed file names and name of TTree. */
TChain* createTChain(const std::vector<std::string>& fileNames, const std::string& treename)
{
  TChain* inChain = new TChain(treename.c_str());
  for (const auto& name : fileNames) {
    inChain->Add(name.c_str());
  }
  return inChain;
}

/** Create TChain from passed file name and name of TTree. */
TChain* createTChain(const std::string& filename, const std::string& treename)
{
  TChain* inChain = new TChain(treename.c_str());
  inChain->Add(filename.c_str());
  return inChain;
}




#endif
