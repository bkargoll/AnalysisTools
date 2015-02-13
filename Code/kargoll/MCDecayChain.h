#ifndef MCDecayChain_h
#define MCDecayChain_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class MCDecayChain : public Selection {

 public:
  MCDecayChain(TString Name_, TString id_);
  virtual ~MCDecayChain();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;

};
#endif
