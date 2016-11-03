#ifndef AlignmentAnalysis_h
#define AlignmentAnalysis_h

#include "BeamAnaBase.h"

#include "TString.h"
//#include "Utility.h"
#include "DataFormats.h"
#include "Histogrammer.h"

class TH1;
class AlignmentAnalysis : public BeamAnaBase {
 public:
  AlignmentAnalysis();
  ~AlignmentAnalysis();
  void beginJob();
  void eventLoop(); 
  std::pair<float, float>  GetOffsetVsZ(const char*, float **, float**);
  void bookHistograms();
  void clearEvent();
  void endJob();
  static double FuncStepGaus(Double_t * x, Double_t * par){
    double xx = x[0];
    double pitch = par[0];
    double sigma = par[1];
    double norm = par[2];
    double cte = par[3];
    double f =0;
    if (xx<0) f = norm*(1+TMath::Erf((xx+pitch/2.)/(sqrt(2)*sigma)));
    if (xx>0) f = norm*(1-TMath::Erf((xx-pitch/2.)/(sqrt(2)*sigma)));
    f += cte;
    return f;
  }
 private:
  std::string runNumber_;
  std::string outFile_;
  Histogrammer* hist_;
  unsigned long int nEntries_;
  bool isProduction_;
  std::string alignparFile_; 
  float zMin;
  float zStep;
  int zNsteps;
};
#endif
