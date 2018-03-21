#ifndef PARALLELTREELOOPER_JN_H
#define PARALLELTREELOOPER_JN_H

/////////////////////////////////////////////////////////////////////////////////
// 
//  PARALLELTREELOOPER j.n. 2018
//  -----------------------------------------------------------------------------
//  A class for looping parallel over a tree and copy it.
//
//  Usage:
//  Branches of in and out tree has to be set in constructor of derived class
//  The function fill_variables() is called in every loop:
//    if there are custom variables, that not simple copy data 
//    from the in to the out tree fill them here
//
//  Some notes (important if you want to modify this base class):
//   - DO NOT clone() AFTER calling the init() function of a ParallelTreeLooper 
//     implementation, this would break the code, if not all pointers of the 
//     branches are set to nullptr after cloning
//   - If you use std::ref be sure to reserve enough space for your vectors,
//     because reallocating crashes the program, std::ref values
//   - Dont forget to call ROOT::EnableThreadSafety, 
//     this is needed if you want to open a .root file twice for reading
//   - Before merging cd() the original out_file, otherwise the current directory 
//     is possibly corrupted, and anyway the wrong one
//
/////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <mutex>
#include <vector>

class TTree;
class TFile;
class TChain;


class ParallelTreeLooper {

public:
  using cs = const std::string &;
  ParallelTreeLooper(const ParallelTreeLooper&) = default;
  ParallelTreeLooper &operator=(const ParallelTreeLooper&) = default;

  ParallelTreeLooper(const std::vector<std::string> &in_file_names, cs intreename, cs outfilename, cs outtreename);

  ParallelTreeLooper(TTree* in_tree, TTree* out_tree) :
    m_in_tree(in_tree),
    m_out_tree(out_tree)
  {
    errorCheck();
  }

  virtual ~ParallelTreeLooper();

  virtual void loop(int max_events = -1, int nThreads = 4);  //-1  means all events

  bool hasError() { return (!m_in_tree || !m_out_tree); }

  // a helper function to print all branches of a tree, so that you can easily copy them into your ParallelTreeLooper implementation
  static void printVariables(cs file, cs tree); 
  
  void worker(long long start_event, long long end_event, int worker_id, long long &count_out, 
    std::string &worker_filename_out, std::mutex &cout_lock, std::mutex &clone_lock, long long update_each, const std::string &uuid);

protected:
  ParallelTreeLooper() {};
  TTree* m_in_tree = nullptr;
  TTree* m_out_tree = nullptr;
  TFile* get_outfile() { return m_out_file; }

  TTree* openOutTree(cs filename, cs treename);
  TTree* openChain(const std::vector<std::string> & file_names, cs tree_name);

  void errorCheck();

  /////////////////////////////
  // HAVE TO BE IMPLEMENTED: //
  /////////////////////////////
  virtual bool fill_and_cut_variables() = 0; // If there are custom variables, that not simple copy data from one to the other tree fill them here
  virtual void init() = 0; // Set all in and out branches here, "virtual" constructor idiom: https://isocpp.org/wiki/faq/virtual-functions#virtual-ctors
  virtual ParallelTreeLooper* clone() const = 0; // Has to return a new instance of the derived ParallelTreeLooper class
  /////////////////////////////

private:
  TFile* m_out_file = nullptr;
  TChain* m_chain = nullptr;
  long long loop_multithreaded(long long nEvents, long long nThreads);
  long long loop_singlethreaded(long long nEvents);
  
  int progress = 0;
  int progress_step = 1;
  void progress_update();
  bool already_looped = false;
  
};
#endif //PARALLELTREELOOPER_JN_H