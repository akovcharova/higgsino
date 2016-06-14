#include <iostream>
#include <vector>

#include "TString.h"

#include "styles.hpp"
#include "timer.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

namespace {
  TString luminosity="20";
  TString plot_type=".pdf";
  TString plot_style="RA4";
  //options: ttbar, other, all_but_qcd, all
  TString sample_set = "ttbar";
  TString nb_bins = "TTML";
  bool loose = false; //remove iso track vero, dr_max (skim still has 2.5), low_dphi
}

using namespace std;

int main(){ 
  time_t begtime, endtime;
  time(&begtime);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  
  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_higloose/");
  if (loose) folder.ReplaceAll("merged_higloose","merged_met100nb2nj4nl0");

  vector<TString> smpl;
  if (sample_set=="ttbar" || sample_set=="all_but_qcd" || sample_set=="all") {
    smpl.push_back(folder+"*_TTJets*Lept*");
    smpl.push_back(folder+"*_TTJets*HT*");
  }
  if (sample_set=="other" || sample_set=="all_but_qcd" || sample_set=="all"){
    smpl.push_back(folder+"*_WJetsToLNu*");
    smpl.push_back(folder+"*ZJetsToNuNu_HT*");
    smpl.push_back(folder+"*ST_tW*");
    smpl.push_back(folder+"*DYJetsToLL_M-50_HT*");
    smpl.push_back(folder+"*TTTT*");
    smpl.push_back(folder+"*_TTWJets*");
    smpl.push_back(folder+"*_WH_*");
    smpl.push_back(folder+"*_ttHJet*");
    smpl.push_back(folder+"*_ZH_*");
    smpl.push_back(folder+"*_ST_*channel*");
    smpl.push_back(folder+"*_TTGJets*");
    smpl.push_back(folder+"*_TTZTo*");
  }
  if (sample_set=="all") smpl.push_back(folder+"*_QCD_HT*");

  TString baseline("pass&&stitch&&njets>=4&&njets<=5&&!low_dphi&&ntks==0&&hig_drmax<2.2");  
  if (loose) baseline = "pass&&stitch&&njets>=4&&njets<=5&&hig_drmax<2.5";  
  // TString sigcut="hig_am>100&&hig_am<140&&hig_dm<40";
  // TString sbdcut="sbd";
  TString cut2b="nbt==2&&nbm==2", cut3b="nbt>=2&&nbm==3&&nbl==3", cut4b="nbt>=2&&nbm>=3&&nbl>=4";
  
  if(nb_bins=="TTTL"){
    cut2b = "nbt==2";
    cut3b = "nbt==3&&nbl==3";
    cut4b = "nbt==3&&nbl>=4";
  }
  if(nb_bins=="MMMM"){
    cut2b = "nbm==2";
    cut3b = "nbm==3";
    cut4b = "nbm>=4";
  }

  if (loose) cut3b.ReplaceAll("==3",">=3"); 

  vector<TString> metbins = {"met>100&&met<=150","met>150&&met<=250","met>250"};

  // Reading ntuples
  vector<sfeats> Samples; 
  TString label = "t#bar{t}";
  if (sample_set=="all") label = "Bkgd";
  if (sample_set=="all_but_qcd") label = "Bkgd (w/o QCD)"; 
  if (sample_set=="other") label = "Other (w/o QCD)";
  Samples.push_back(sfeats(smpl, label+", N_{b}=2", kGreen+1,1, "stitch&&pass&&"+cut2b));
  Samples.push_back(sfeats(smpl, label+", N_{b}=3", kAzure-3,1, "stitch&&pass&&"+cut3b)); 
  Samples.back().mcerr = true;
  Samples.push_back(sfeats(smpl, label+", N_{b}=4", kRed-4,1, "stitch&&pass&&"+cut4b)); 
  Samples.back().mcerr = true;
  Samples.back().style = 2;

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples

  vector<hfeats> vars;
  TString tag = sample_set+"_"+nb_bins;
  if (loose) tag += "_loose";
  for (size_t imet(0); imet < metbins.size(); imet++){
    float coarse = (imet==metbins.size()-1) ? 2:1;
    vars.push_back(hfeats("hig_dm",32/coarse,0,160,ra2b_sam, "#Deltam [GeV]",baseline+"&&"+metbins[imet],40, tag));
    vars.back().whichPlots = "34"; 
    vars.push_back(hfeats("hig_am",24/coarse,0,240,ra2b_sam, "<m> [GeV]",baseline+"&&"+metbins[imet], 100, tag));
    vars.back().whichPlots = "34"; vars.back().cut2 = 140;
    // vars.push_back(hfeats("hig_drmax",20,0,4,ra2b_sam, "#DeltaR_{max}",baseline+"&&"+metbins[imet],2.2, tag));
    // vars.back().whichPlots = "34"; 
    // vars.push_back(hfeats("Sum$(jets_hflavor==5)",7,0.5,7.5,ra2b_sam, "Number of true b-jets",baseline+"&&"+metbins[imet],-1, tag));
    // vars.back().whichPlots = "34"; 
  }


  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, tag,false,true);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
