#include "ParallelTreeLooper.h"
#include "utils.h"

#include <iomanip>
#include <chrono>
#include <thread>
#include <cstdio>


#ifndef __CLING__
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TROOT.h"
#endif

void ParallelTreeLooper::loop(int max_events, int nThreads) {

  if (already_looped) {
    std::cout << "Sorry, but at the moment loop() can be called only once for the same instance of a ParallelTreeLooper." << std::endl;
    return;
  }
  already_looped = true;

  if (hasError()) {
    std::cout << "ERROR TreeLooper::loop() : not starting because there are problems with the trees." << std::endl;
    return;
  }

  if (!(nThreads > 0 && nThreads < 9)) nThreads = 4;

  const long long nInputEvents = (m_in_tree->GetEntries());
  const size_t nEvents = (max_events < 0 || max_events > nInputEvents) ? nInputEvents : max_events;

  std::cout << "Looping over " << nEvents << " events..." << std::endl;

  ///////////////////
  // START LOOPING //
  ///////////////////
  long long count = 0;
  progress = 0;
  stopwatch w("ParallelTreeLooper");

  if (nThreads > 1) {
    count = loop_multithreaded(nEvents, nThreads);
  }
  else {
    count = loop_singlethreaded(nEvents);
  }

  ////////////////

  std::cout << "number of reconstructed events: " << count << " of a total of " << nEvents << " events\n";

  auto of = m_out_tree->GetCurrentFile();
  if (of) of->cd();
  else std::cout << "ERROR TreeLooper: Could not write output file." << std::endl;
  m_out_tree->Write("", TObject::kWriteDelete);
  if (of) of->Write("", TObject::kWriteDelete);

  auto outfilename(of ? of->GetName() : "NONE");
  auto infilename(m_in_tree->GetCurrentFile() ? m_in_tree->GetCurrentFile()->GetName() : "NONE"); //TODO: chain has more than one file!

  std::cout << "Created file " << outfilename << /*", input from " << infilename << "." <<*/ std::endl;

}

void ParallelTreeLooper::printVariables(cs file, cs tree) {
  //open file
  TFile f(file.c_str());
  if (f.IsZombie()) {
    std::cout << "printVariables: Could not open file " << file << std::endl;
    return;
  }

  //open tree
  auto t = (TTree *)f.Get(tree.c_str());
  if (!t) {
    std::cout << "printVariables: Could not open tree " << tree << std::endl;
    return;
  }

  std::stringstream var_ss;
  std::stringstream in_branch_ss;
  std::stringstream out_branch_ss;

  auto branch_list = t->GetListOfBranches();
  for (auto obj : *branch_list) {
    auto b = dynamic_cast<TBranch*>(obj);
    if (!b) continue;

    bool isObject = true;
    std::string type = b->GetClassName();
    if (type.empty()) {
      isObject = false;
      auto l = b->GetLeaf(b->GetName());
      if (!l) type = "UNKNOWN TYPE";
      else type = l->GetTypeName();
    }
    if (isObject) type.append("*");
    auto name = b->GetName();
    var_ss << type << " " << name << (isObject ? " = nullptr;\n" : ";\n");
    in_branch_ss << "in->SetBranchAddress(\"" << name << "\", &" << name << ");\n";
    out_branch_ss << "out->Branch(\"" << name << "\", &" << name << ");\n";
  }
  std::cout << var_ss.rdbuf() << "\n" << in_branch_ss.rdbuf() << '\n' << out_branch_ss.rdbuf() << std::endl;
}


TTree* ParallelTreeLooper::openOutTree(cs filename, cs treename) {
  TFile *f = nullptr;
  TTree *t = nullptr;

  //open file
  f = new TFile(filename.c_str(), "RECREATE");
  if (f->IsZombie()) {
    std::cout << "Could not open file " << filename << std::endl;
    SafeDelete(f);
    return nullptr;
  }

  //open tree
  t = new TTree(treename.c_str(), treename.c_str());
  if (!t) {
    std::cout << "Could not open tree " << treename << std::endl;
    SafeDelete(f);
    return nullptr;
  }
  t->SetDirectory(f);

  delete m_out_file;
  m_out_file = f;
  return t;
}

TTree* ParallelTreeLooper::openChain(const std::vector<std::string> & file_names, cs tree_name) {
  //delete m_chain;
  m_chain = new TChain(tree_name.c_str(), tree_name.c_str());
  for (auto f : file_names) m_chain->Add(f.c_str());
  if (m_chain->GetEntries() < 1) {
    SafeDelete(m_chain);
  }
  return m_chain;
}

long long ParallelTreeLooper::loop_multithreaded(long long nEvents, long long nThreads)
{

  std::mutex cout_lock;
  std::mutex clone_lock;

  ROOT::EnableThreadSafety();

  std::cout << "Processing with " << nThreads << " threads." << std::endl;

  std::vector<std::thread> threads;
  std::vector<std::string> worker_filenames;
  std::vector<long long> counters;
  worker_filenames.reserve(nThreads);
  counters.reserve(nThreads);
  long long events_per_worker = nEvents / nThreads;
  long long update_every = events_per_worker / (100. / nThreads);

  for (long long start = 0, end = events_per_worker, i = 0; i < nThreads; ++i, start = end, end += events_per_worker) {

    if (i == nThreads - 1) end = nEvents;
    counters.emplace_back(0);
    worker_filenames.emplace_back("");

    threads.emplace_back(&ParallelTreeLooper::worker, this, start, end, (int)i, std::ref(counters.back()), std::ref(worker_filenames.back()), std::ref(cout_lock), std::ref(clone_lock), update_every);
  }

  //join all threads, then sum up worker_count variables into count
  for (auto &t : threads) t.join();

  unsigned long long count = 0;
  for (const auto &c : counters) count += c;

  //merge trees from worker loopers and then delete the temporary worker file:
  //https://root-forum.cern.ch/t/merging-ttrees-on-the-fly/1833
  {
    std::cout << "Merging worker trees..." << std::endl;
    TList *trees = new TList;
    std::vector<TFile*> files;
    const std::string out_treename = m_out_tree->GetName();
    for (int i = 0, s = worker_filenames.size(); i < s; ++i) {
      auto filename = worker_filenames.at(i);
      files.emplace_back(TFile::Open(filename.c_str()));
      auto f = files.back();
      if (!f || (f && f->IsZombie())) std::cout << "Problems opening file for merging: " << filename << std::endl;
      auto tree = (TTree*)f->Get(out_treename.c_str());
      if (!tree) std::cout << "Problems reading tree " << out_treename << " in workers file " << filename << std::endl;
      trees->Add(tree);
    }
    m_out_tree->GetCurrentFile()->cd();
    init(); //call init AFTER cloning the Looper for the threads
    m_out_tree->Merge(trees, "fast");
    //m_out_tree->Print();
    m_out_tree->Write("", TObject::kWriteDelete);
    m_out_tree->GetCurrentFile()->Write("", TObject::kWriteDelete);
    for (auto f : files) delete f;
    for (auto &s : worker_filenames) if (std::remove(s.c_str())) std::cout << "Could not delete worker file " << s << " from disk." << std::endl;
    std::cout << "Merging worker trees completed." << std::endl;
  }
  return count;
}

long long ParallelTreeLooper::loop_singlethreaded(long long nEvents)
{
  const long long onepercent = nEvents / 100.;
  unsigned long long count = 0;
  init();
  for (long long i = 0; i < nEvents; ++i) {
    if (m_in_tree->GetEntry(i) < 0) {
      std::cout << "\nI/O error while reading event " << i << " in TTree '" << m_in_tree->GetName() << "'" << std::endl;
      continue;
    }
    if (fill_and_cut_variables()) {
      m_out_tree->Fill();
      ++count;
    }
    if (!(i%onepercent)) progress_boost();
  }

  return count;
}

void ParallelTreeLooper::worker(long long start_event, long long end_event, int worker_id, long long & count_out, std::string & worker_filename_out,
  std::mutex &cout_lock, std::mutex &clone_lock, long long update_every)
{
  ParallelTreeLooper *looper = nullptr;
  {
    std::lock_guard<std::mutex> lock(clone_lock);
    looper = this->clone();
  }
  {
    std::lock_guard<std::mutex> lock(cout_lock);
    std::cout << "Starting thread " << worker_id << " from event " << start_event << " to event " << end_event << '.' << std::endl;
  }

  auto const in_tree = looper->m_in_tree;
  auto const out_tree = looper->m_out_tree;

  const std::string worker_treename_out = out_tree->GetName();
  worker_filename_out = std::string(out_tree->GetCurrentFile()->GetName());
  worker_filename_out = worker_filename_out.substr(0, worker_filename_out.size() - 5);
  worker_filename_out += "_wrk" + std::to_string(worker_id) + "_" +
    std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count()) + //unique ID, to avoid conflicts
    ".root";

  // Read input filenames:
  auto in_tree_name = in_tree->GetName();
  std::vector<std::string> in_file_names;
  auto infile = dynamic_cast<TChain*>(in_tree);
  if (infile) {
    TObjArray *fileElements = infile->GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl = nullptr;
    while ((chEl = (TChainElement*)next())) {
      in_file_names.emplace_back(chEl->GetTitle());
    }
  }
  /*
  else {//assuming that it is a simple tree if it is no chain:
    in_file_names.emplace_back(in_tree->GetCurrentFile()->GetName());
  }
*/
  {
    std::lock_guard<std::mutex> lock(cout_lock);
    std::cout << "Starting to init looper " << worker_filename_out << " " << worker_treename_out << std::endl;
  }

  looper->m_out_file = nullptr; //that original out file is not closed
  looper->m_chain = nullptr; //that original in chain is not closed
  looper->m_out_tree = looper->openOutTree(worker_filename_out, worker_treename_out);
  looper->m_in_tree = looper->openChain(in_file_names, in_tree_name);
  looper->init();

  {
    std::lock_guard<std::mutex> lock(cout_lock);
    std::cout << "Temporary file for thread " << worker_id << ": " << looper->m_out_tree->GetCurrentFile()->GetName() << std::endl;
  }

  count_out = 0;
  for (long long i = start_event; i < end_event; ++i) {
    if (looper->m_in_tree->GetEntry(i) < 0) {
      std::lock_guard<std::mutex> lock(cout_lock);
      std::cout << "I/O error while reading event " << i << " in TTree '" << looper->m_in_tree->GetName() << "'" << std::endl;
      continue;
    }

    if (looper->fill_and_cut_variables()) {
      looper->m_out_tree->Fill();
      ++count_out;
    }

    if (!(i%update_every)){
      std::lock_guard<std::mutex> lock(cout_lock);
      this->progress_boost();
    }

  }
  looper->m_out_tree->Write("", TObject::kWriteDelete);
  looper->m_out_file->Write("", TObject::kWriteDelete);
  {
    std::lock_guard<std::mutex> lock(cout_lock);
    std::cout << "TreeLooper worker " << worker_id << " processed " << count_out << " events." << std::endl;
  }
  //looper->m_out_tree->Print();

  delete looper; //also closes files

}

void ParallelTreeLooper::progress_boost()
{
  static std::mutex progress_lock;
  static const int prog_w = 6;
  static bool firstrun = true;
  if (firstrun) firstrun = false, std::cout << std::string(prog_w,' ');

  {
    std::lock_guard<std::mutex> lock(progress_lock);
    progress += progress_step;
  }
  std::cout << std::string(prog_w, '\b') << std::setw(prog_w) << "\n" + std::to_string(progress) + " %" << (progress==100?"\n":"") << std::flush;
}

ParallelTreeLooper::~ParallelTreeLooper() {
  delete m_out_file;
  delete m_chain;
}