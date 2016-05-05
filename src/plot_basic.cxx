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
}

using namespace std;

int main(){ 
  time_t begtime, endtime;
  time(&begtime);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder
  
  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/merged_met100nb2nj4nl0/");


  vector<TString> s_tchi;
  s_tchi.push_back(folder+"*TChiHH*mChi-400_*");
  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*Lept*");
  vector<TString> s_wjets;
  s_wjets.push_back(folder+"*_WJetsToLNu*");
  vector<TString> s_znunu;
  s_znunu.push_back(folder+"*ZJetsToNuNu_HT*");
  vector<TString> s_wt;
  s_wt.push_back(folder+"*ST_tW*");
  vector<TString> s_qcd;
  s_qcd.push_back(folder+"*_QCD_HT*");
  vector<TString> s_other;
  s_other.push_back(folder+"*DYJetsToLL_M-50_HT*");
  s_other.push_back(folder+"*TTTT*");
  s_other.push_back(folder+"*_TTWJets*");
  s_other.push_back(folder+"*_WH_*");
  s_other.push_back(folder+"*_ttHJet*");
  s_other.push_back(folder+"*_ZH_*");
  s_other.push_back(folder+"*_ST_*channel*");
  s_other.push_back(folder+"*_TTGJets*");
  s_other.push_back(folder+"*_TTZTo*");


  // Reading ntuples
  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_tchi, "TChiHH(400,1)",	 ra2b::c_tchi, 1)); Samples.back().isSig = true;
  Samples.push_back(sfeats(s_tt, "t#bar{t}",		 ra2b::c_tt_1l));
  Samples.push_back(sfeats(s_qcd, "QCD",		 ra2b::c_qcd));
  //Samples.push_back(sfeats(s_znunu, "Z#rightarrow#nu#nu",ra2b::c_znunu));
  Samples.push_back(sfeats(s_wt, "W+top",		 ra2b::c_wt));
  Samples.push_back(sfeats(s_wjets, "W+jets",		 ra2b::c_wjets));
  Samples.push_back(sfeats(s_other, "Other",		 ra2b::c_other)); 

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples



  vector<hfeats> vars;


  //// Baseline
  TString base = "nbl>=4&&nbm>=2&&njets<=5&&ntks==0&&hig_drmax<2.2&&!low_dphi", bin = "&&hig_am>100&&hig_am<140&&hig_dm<40";
  TString metcut = "&&met>250";

  vars.push_back(hfeats("jets_pt[0]",30,0,600,ra2b_sam, "Jet 1 p_{T} [GeV]",base+bin+metcut, 50));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[1]",35,0,350,ra2b_sam, "Jet 2 p_{T} [GeV]",base+bin+metcut, 50));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[2]",25,0,250,ra2b_sam, "Jet 3 p_{T} [GeV]",base+bin+metcut));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[3]",20,0,200,ra2b_sam, "Jet 4 p_{T} [GeV]",base+bin+metcut));
  vars.back().whichPlots = "12"; 

  vars.push_back(hfeats("3.14.16+hig_dphi",40,0,4,ra2b_sam, "#Delta#phi [GeV]",base+bin+metcut,0.5));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("hig_dm",40,0,200,ra2b_sam, "#Deltam [GeV]",base+"&&hig_am>100&&hig_am<140"+metcut,40));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("hig_am",25,0,250,ra2b_sam, "<m> [GeV]",base+"&&hig_dm<40"+metcut));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("hig_drmax",30,0,6,ra2b_sam, "#DeltaR_{max}","nbl>=4&&nbm>=2&&njets<=5&&ntks==0"+bin+metcut,2.2));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("met",20,100,1100,ra2b_sam, "E_{T}^{miss} [GeV]",base+bin,250));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("njets",10,3.5,13.5,ra2b_sam, "N_{jet}","nbl>=4&&nbm>=2&&ntks==0&&hig_drmax<2.2"+bin+metcut,5.5));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbl",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{L}","njets<=5&&ntks==0&&hig_drmax<2.2"+bin+metcut,3.5));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbm",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{M}","njets<=5&&ntks==0&&hig_drmax<2.2"+bin+metcut,2.5));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbt",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{T}","njets<=5&&ntks==0&&hig_drmax<2.2"+bin+metcut,1.5));
  vars.back().whichPlots = "12"; 

  // vars.push_back(hfeats("ht",34,0,3400,ra2b_sam, "H_{T} [GeV]",base+bin+metcut));
  // vars.back().whichPlots = "12"; 


  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, "basic5",false,true);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
