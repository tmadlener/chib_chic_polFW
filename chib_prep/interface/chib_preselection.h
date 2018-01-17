#include "TreeProcessor.h"

class ChibPreselection : public TreeProcessor
{

public:
  using TreeProcessor::TreeProcessor;

private:
  virtual bool fill_and_cut_variables() override;
  virtual ParallelTreeLooper* clone() const override { return new ChibPreselection(*this); }
  virtual void setup_new_branches() override;
};