#ifndef Parameters_h
#define Parameters_h

#include <vector>
#include "TString.h"

class Parameters {

 public:
  Parameters();
  Parameters(TString f, bool verbosity = false);
  virtual ~Parameters();

  void SetVerbosity(bool verbosity){verbose = verbosity;}


  void    SetFile(TString f);
  TString GetFile();

  void GetString(TString p, TString &v, TString dv="");
  void GetBool(TString p, bool &v, bool dv=false);
  void GetInt(TString p, int &v, int dv=0);
  void GetDouble(TString p, double &v, double dv=0.0);
  void GetVectorString(TString p, std::vector<TString> &v, TString dv="");
  void GetVectorStringDouble(TString p, std::vector<TString> &v1, std::vector<double> &v2);

 private:
  static TString file;
  static bool verbose;

  template<typename T> void GetParameter(TString p, T &v, T dv);

};
#endif
