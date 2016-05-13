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
  TString folder_ns(bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/unskimmed/");


  vector<TString> s_tchi;
  s_tchi.push_back(folder+"*TChiHH*mChi-400_*");
  vector<TString> s_tchi_ns;
  s_tchi_ns.push_back(folder_ns+"*TChiHH*mChi-400_*");
  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*Lept*");
  s_tt.push_back(folder+"*_TTJets*HT*");
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
  Samples.push_back(sfeats(s_tchi, "TChiHH(400,1)",	 ra2b::c_tchi, 1, "pass")); Samples.back().isSig = true;
  Samples.push_back(sfeats(s_tt, "t#bar{t}",		 ra2b::c_tt_1l,1, "stitch&&pass"));
  //Samples.push_back(sfeats(s_tt, "t#bar{t}",		 ra2b::c_tt_1l, 1, "pass"));
  Samples.push_back(sfeats(s_qcd, "QCD",		 ra2b::c_qcd, 1, "pass"));
  //Samples.push_back(sfeats(s_znunu, "Z#rightarrow#nu#nu",ra2b::c_znunu, 1, "pass"));
  Samples.push_back(sfeats(s_wt, "W+top",		 ra2b::c_wt, 1, "pass"));
  Samples.push_back(sfeats(s_wjets, "W+jets",		 ra2b::c_wjets, 1, "pass"));
  Samples.push_back(sfeats(s_other, "Other",		 ra2b::c_other, 1, "pass")); 

  vector<int> ra2b_sam;
  unsigned nsam(Samples.size());
  for(unsigned sam(0); sam < nsam; sam++){
    ra2b_sam.push_back(sam);
  } // Loop over samples

  size_t ihig = Samples.size();
  Samples.push_back(sfeats(s_tchi_ns, "True H_{1}", ra2b::c_tchi, 1, "1", "Max$(mc_pt*(mc_id==25))")); 
  Samples.push_back(sfeats(s_tchi_ns, "True H_{2}", ra2b::c_tchi, 2, "1", "Min$(mc_pt*(mc_id==25)+99999*(mc_id!=25))")); 
  Samples.push_back(sfeats(s_tchi_ns, "Reco H_{1}", 4, 1, "1", "max(hig1_pt, hig2_pt)")); 
  Samples.push_back(sfeats(s_tchi_ns, "Reco H_{2}", 4, 2, "1", "min(hig1_pt, hig2_pt)")); 
  vector<int> hig_sam;
  for(unsigned sam(ihig); sam < Samples.size(); sam++){
    hig_sam.push_back(sam);
  } // Loop over samples


  vector<hfeats> vars;

  //// Higgs pT
  vars.push_back(hfeats("hig_pt",16,0,800,hig_sam, "Higgs p_{T} for TChiHH(400,1) [GeV]","njets>=4&&nbm>=4&&met<250",350));
  vars.back().whichPlots = "3";  
  vars.push_back(hfeats("hig_pt",16,0,800,hig_sam, "Higgs p_{T} for TChiHH(400,1) [GeV]","njets>=4&&nbm>=4&&met>250",350));
  vars.back().whichPlots = "3";  

  //// Sample to plot
  TString met_s = "himet", nb_s = "4b", tag = met_s+"_"+nb_s;

  //// Baseline
  TString bcut = "&&nbt>=2&&nbm>=3&&nbl>=4", jcut = "&&njets>=4&&njets<=5", drcut = "&&hig_drmax<2.2", dphicut = "&&!low_dphi";
  TString tkscut = "&&ntks==0", metcut = "&&met>250";
  TString bin = "&&hig_am>100&&hig_am<140&&hig_dm<40";

  if(nb_s=="23b") bcut = "&&nbt>=2&&nbl<=3";
  else bcut = "&&nbt>=2&&nbm>=3&&nbl>=4";
  if(met_s=="lowmet") metcut = "&&met>150&&met<=250";
  else metcut = "&&met>250";

  TString baseb = jcut+drcut+dphicut+tkscut+metcut, basephi = jcut+drcut+bcut+tkscut+metcut;
  TString basemet = bcut+jcut+drcut+dphicut+tkscut, basetks = bcut+jcut+drcut+dphicut+metcut;;
  TString basejets = bcut+drcut+dphicut+tkscut+metcut, basedr = bcut+jcut+dphicut+tkscut+metcut;
  TString base = bcut+jcut+drcut+dphicut+tkscut+metcut;


  vars.push_back(hfeats("met",20,100,600,ra2b_sam, "E_{T}^{miss} [GeV]",basemet+bin,150, tag));
  vars.back().whichPlots = "12";  vars.back().cut2 = 250;
  vars.push_back(hfeats("ht",28,0,1400,ra2b_sam, "H_{T} [GeV]",base+bin,-1, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("njets",10,3.5,13.5,ra2b_sam, "N_{jet}",basejets+bin,5.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbl",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{L}",baseb+bin,3.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbm",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{M}",baseb+bin,2.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("nbt",7,-0.5,6.5,ra2b_sam, "N_{b-jet}^{T}",baseb+bin,1.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("ntks",5,-0.5,4.5,ra2b_sam, "N_{track}",basetks+bin,0.5, tag));
  vars.back().whichPlots = "12"; 

  vars.push_back(hfeats("hig_dm",32,0,160,ra2b_sam, "#Deltam [GeV]",base+"&&hig_am>100&&hig_am<140",40, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("hig_am",25,0,250,ra2b_sam, "<m> [GeV]",base+"&&hig_dm<40", 100, tag));
  vars.back().whichPlots = "12"; vars.back().cut2 = 140;
  vars.push_back(hfeats("hig5_am",25,0,250,ra2b_sam, "Top 5 CSV <m> [GeV]",base+"&&hig_dm<40", 100, tag));
  vars.back().whichPlots = "12"; vars.back().cut2 = 140;
  vars.push_back(hfeats("hig_drmax",20,0,4,ra2b_sam, "#DeltaR_{max}",basedr+bin,2.2, tag));
  vars.back().whichPlots = "12"; 

  vars.push_back(hfeats("3.1416+hig_dphi",32,0,3.2,ra2b_sam, "Min #Delta#phi [GeV]",basephi+bin,0.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("dphi1",32,0,3.2,ra2b_sam, "#Delta#phi_{1} [GeV]",basephi+bin,0.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("dphi2",32,0,3.2,ra2b_sam, "#Delta#phi_{2} [GeV]",basephi+bin,0.5, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("dphi3",32,0,3.2,ra2b_sam, "#Delta#phi_{3} [GeV]",basephi+bin,0.3, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("dphi4",32,0,3.2,ra2b_sam, "#Delta#phi_{4} [GeV]",basephi+bin,0.3, tag));
  vars.back().whichPlots = "12"; 

  vars.push_back(hfeats("jets_pt[0]",30,0,600,ra2b_sam, "Jet 1 p_{T} [GeV]",base+bin, 50, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[1]",35,0,350,ra2b_sam, "Jet 2 p_{T} [GeV]",base+bin, 50, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[2]",25,0,250,ra2b_sam, "Jet 3 p_{T} [GeV]",base+bin, 30, tag));
  vars.back().whichPlots = "12"; 
  vars.push_back(hfeats("jets_pt[3]",20,0,200,ra2b_sam, "Jet 4 p_{T} [GeV]",base+bin, 30, tag));
  vars.back().whichPlots = "12"; 


  plot_distributions(Samples, vars, luminosity, plot_type, plot_style, tag,false,true);

  time(&endtime); 
  cout<<endl<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
