//----------------------------------------------------------------------------
// utilities_macros - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef H_UTILITIES_MACROS
#define H_UTILITIES_MACROS

#include <vector>

#include "TString.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TChain.h"
#include "baby_full.hpp"
#include "bcut.hpp"

namespace ra2b{
  // Had to define the TColor objects in the cpp
  enum {
    c_t1t      = kGreen,   
    c_tchi     = 2,    // red
    c_tt_1l    = 2012, // ucsb_blue
    c_tt_2l    = 2011, // tar_heel_blue
    c_wjets    = 2004, // uo_green
    c_znunu    = 1, // penn_red
    c_wt       = 2002, // penn_red
    c_qcd      = 2001, // ucsb_gold
    c_other    = 2005  // tcu_purple
  };

}


struct pfeats{
  pfeats(const std::vector<int> &isamples, const TString &icut = "1", const TString &itagname="");

  std::vector<int> samples;
  TString cut, tagname;
};

class hfeats {
public:
  hfeats(TString ivarname, int inbins, float iminx, float imaxx, std::vector<int> isamples,
         TString ititle="", TString icuts="1", float icut=-1, TString itagname="",bool iskiplog=false, 
         std::vector<double> inevents= std::vector<double>(1,-1.));
  hfeats(TString ivarname, int inbins, float* ibinning, std::vector<int> isamples,
         TString ititle="", TString icuts="1", float icut=-1, TString itagname="",bool iskiplog=false, 
         std::vector<double> inevents= std::vector<double>(1,-1.));
  hfeats(TString ivarnamex, TString ivarnamey, int inbinsx, float iminx, float imaxx, int inbinsy, float iminy, float imaxy,  std::vector<int> isamples,
         TString ititlex, TString ititley, TString icuts, float icutx, float icuty, TString itagname);
  TString title, titlex, titley, varname, varnamex, varnamey, tag, cuts, unit;
  int nbins, nbinsx, nbinsy;
  float *binning;
  float minx, maxx, miny, maxy, cut, cutx, cuty,  maxYaxis, maxRatio;
  std::vector<int> samples;
  TString tagname;
  void format_tag();
  std::vector<double> nevents; //Added for track veto study. Useful to display number of events when hist is filled N times per event
  bool skiplog;
  TString whichPlots; // String that determines which of the [log_]lumi and [log_]shapes plots to make
  bool normalize; //normalizes isData to sum of histograms
  bool PU_reweight; 
  float moveRLegend; // Moves the right legend by this amount
  bool addOverflow; // If true, the overflow is added to the last bin
  
};

class sfeats {
public:
  sfeats(std::vector<TString> ifile, TString ilabel, int icolor=1, int istyle=1, TString icut="1",
	 TString samVariable="noPlot");
  std::vector<TString> file;
  TString label, cut, factor,tag;
  int color, style;
  bool isSig, doStack, isData, mcerr, doBand;
  TString samVariable; // Used to plot different variables in the same histogram
};

class tfeats {
public:
  tfeats(TString icuts, TString itag);
  void add(TString texname, TString tcut, TString option="");
  TString cuts, tag;
  std::vector<TString> texnames, tcuts, options;
};

class sysfeats {
public:
  sysfeats(TString iname, TString ititle);
  TString name;
  TString title;
  std::vector<TString> bincuts;
  std::vector<double> weights;
  void push_back(TString bincut, double weight);
  TString bincut(unsigned i);
  double weight(unsigned i);
  unsigned size();
};


void calc_chi2_diff(TH1D *histo1, TH1D *histo2, float &chi2, int &ndof, float &pvalue, float *average);
void calc_chi2(TH1D *histo, float &chi2, int &ndof, float &pvalue, float &average);
void dump_event(baby_full &tree, int entry);

std::vector<double> getYields(baby_full &baby, bcut baseline, std::vector<bcut> bincuts, std::vector<double> &yields, 
            std::vector<double> &w2, double lumi=1., bool do_trig=false);
std::vector<std::vector<double> > getYields(baby_full &baby, bcut baseline, std::vector<std::vector<bcut> > bincuts, 
              std::vector<std::vector<double> > &yields, 
              std::vector<std::vector<double> > &w2, double lumi=1., bool do_trig=false);

long getYieldErr(TChain& tree, TString cut, double& yield, double& uncertainty);

void plot_distributions(std::vector<sfeats> Samples, std::vector<hfeats> vars, TString luminosity="10", 
			TString filetype=".eps", TString namestyle="LargeLabels", TString dir = "1d", 
			bool doRatio=false, bool showcuts=false);
void plot_2D_distributions(std::vector<sfeats> Samples, std::vector<hfeats> vars, TString luminosity,
                           TString filetype, TString namestyle, TString dir);
TString cuts2title(TString title);
TString cuts2tex(TString cuts);
TString invertcut(TString cut);
TString format_tag(TString tag);
double gsl_ran_gamma (const double a, const double b, TRandom3 &rand);
double intGaus(double mean, double sigma, double minX, double maxX);
// yields[Nobs][Nsam] has the entries for each sample for each observable going into kappa
// weights[Nobs][Nsam] has the average weight of each observable for each sample
// powers[Nobs] defines kappa = Product_obs{ Sum_sam{yields[sam][obs]*weights[sam][obs]}^powers[obs] }
double calcKappa(std::vector<std::vector<float> > &entries, std::vector<std::vector<float> > &weights,
		 std::vector<float> &powers, float &mSigma, float &pSigma, bool do_data=false, 
		 bool verbose=false, bool do_plot=false, int nrep=100000);
float Efficiency(float den, float num, float &errup, float &errdown);

#endif
