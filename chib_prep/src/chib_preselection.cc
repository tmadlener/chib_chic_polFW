#include "chib_preselection.h"

bool ChibPreselection::fill_and_cut_variables()
{
  return true;
}

void ChibPreselection::setup_new_branches()
{
}


int main() {
  ChibPreselection preselection({ "/afs/hephy.at/work/j/jnecker/data/full/chib_2016.root" }, "testdata_chib_preselection.root", "rootuple/chiTree", "data");
  preselection.AddBranchesToCopy({ "chi_p4" });
  preselection.process(100000, 4);
}