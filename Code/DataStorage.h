#ifndef DataStorage_h
#define DataStorage_h

#include <vector>
#include "TString.h"

class DataStorage {
 public:
  DataStorage();
  virtual ~DataStorage();

  int  GetFile(TString key);
  void StoreFile(TString File, TString savedFile);

  TString timeStamp(bool Date = true, bool Time = true);

 private:
  TString mydir;

  // list of temporary files to be deleted
  std::vector<TString> filesToDelete;

 protected:
  TString inputFileName;
};
#endif
