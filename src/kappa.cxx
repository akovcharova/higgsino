#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_full.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  double lumi(20.);
  TString nb_bins = "TML";
  bool loose = true; //remove iso track vero, dr_max (skim still has 2.5), low_dphi
}

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2016_04_29/mc/merged_higloose/";

  ////// Creating babies
  baby_full tt(folder+"*TTJets*Lept*");
  tt.Add(folder+"*TTJets*HT*");

  ////// Defining cuts

  TString baseline("pass&&stitch&&njets>=4&&njets<=5&&!low_dphi&&hig_drmax<2.2&&ntks==0");  
  if (loose) baseline = "pass&&stitch&&njets>=4&&njets<=5&&hig_drmax<2.5";  

  TString sigcut="hig_am>100&&hig_am<140&&hig_dm<40";
  TString sbdcut="sbd";
  TString cut2b="nbt==2&&nbm==2", cut3b="nbt>=2&&nbm==3&&nbl==3", cut4b="nbt>=2&&nbm>=3&&nbl>=4";
  
  if(nb_bins=="TTL"){
    cut2b = "nbt==2";
    cut3b = "nbt==3&&nbl==3";
    cut4b = "nbt==3&&nbl>=4";
  }
  if(nb_bins=="MMM"){
    cut2b = "nbm==2";
    cut3b = "nbm==3";
    cut4b = "nbm>=4";
  }
 
  vector<TString> metbins = {"met>100&&met<=150","met>150&&met<=250","met>250"};

  vector<vector<TString> > abcds; vector<TString> abcd_names;
  abcds.push_back(vector<TString>({sigcut+"&&"+cut2b, sbdcut+"&&"+cut2b, sigcut+"&&"+cut3b, sbdcut+"&&"+cut3b})); 
  abcd_names.push_back("3b");
  abcds.push_back(vector<TString>({sigcut+"&&"+cut2b, sbdcut+"&&"+cut2b, sigcut+"&&"+cut4b, sbdcut+"&&"+cut4b}));
  abcd_names.push_back("4b");

  ////// Combining cuts
  vector<bcut > bincuts;
  vector<bool> do_sbd;
  for(size_t ind(0); ind<metbins.size(); ind++){
    for (size_t ipl(0); ipl<abcds.size();ipl++){
      for(size_t obs(0); obs < abcds[ipl].size(); obs++){
        TString totcut(abcds[ipl][obs]+"&&"+metbins[ind]);
        if (totcut.Contains("sbd")) {
          do_sbd.push_back(true);
          totcut.ReplaceAll("&&sbd&&","&&").ReplaceAll("sbd&&","").ReplaceAll("&&sbd","").ReplaceAll("sbd","");
        } else do_sbd.push_back(false);
        bincuts.push_back(bcut(totcut));
      } // Loop over observables going into kappa
    } // Loop over signal bins
  }
    
  ////// Finding yields 
  vector<double> mcyield, mcw2;
  getYields(tt, baseline, bincuts, mcyield, mcw2, lumi, false, do_sbd);
  // for (size_t i=0; i< mcyield.size(); i++) 
  //   cout<<mcyield[i]<<endl;

  vector<float> pow_tot;
  pow_tot.push_back(1); 
  pow_tot.push_back(-1);
  pow_tot.push_back(-1);
  pow_tot.push_back(1); 

  float mSigma, pSigma, kappa;
  vector<vector<float> > kappas;
  for(size_t imet(0); imet<metbins.size(); imet++){
    for(size_t ipl(0); ipl<abcds.size(); ipl++){
      vector<TString> abcdcuts = abcds[ipl];
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      float k(1.);    
      for(size_t obs(0); obs < abcdcuts.size(); obs++){
        size_t index(abcds.size()*abcdcuts.size()*imet+ipl*abcdcuts.size()+obs);
        k *= pow(mcyield[index], pow_tot[obs]);
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(pow(mcyield[index],2)/mcw2[index]);
        weights.back().push_back(mcw2[index]/mcyield[index]);

      } // Loop over observables for MC
      kappa = calcKappa(entries, weights, pow_tot, mSigma, pSigma);
      if(mSigma<0) mSigma = 0;
      kappas.push_back(vector<float>({k, mSigma, pSigma}));

    } // Loop over signal bins
  }

  ///// Printing table
  TString outname = "txt/table_kappa_"+nb_bins+".tex";
  if (loose) outname.ReplaceAll(".tex","_loose.tex");
  ofstream out(outname);

  size_t digits(1);
  out << fixed << setprecision(digits);
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating,geometry}\n";
  out << "\\renewcommand{\\arraystretch}{1.3}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  out << "\\resizebox{\\textwidth}{!}{\n";

  out << "\n\\begin{tabular}[tbp!]{ l rrr }\\hline\\hline\n";
  out << "Bin & $\\kappa$ & $\\Delta\\kappa$ [\\%] & $\\Delta\\kappa$ [\\%]\\\\ \\hline"<<endl;
  for(size_t imet(0); imet<metbins.size(); imet++){
    out << " & \\multicolumn{3}{c}{"<<cuts2tex(metbins[imet])<<"}  \\\\ \\hline\n";
    for(size_t ipl(0); ipl<abcds.size(); ipl++){
      size_t ind(abcds.size()*imet+ipl);
      out<<abcd_names[ipl]<<" & "<<RoundNumber(kappas[ind][0],2) <<" & "<<kappas[ind][1] <<" & "<<kappas[ind][2];
        out<<" \\\\"<<endl;
      out<<"\\hline"<<endl;
    }
  }

  out<< "\\hline\\hline\n\\end{tabular}"<<endl<<endl;
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out.close();
  TString pdfname(outname); 
  cout<<" pdflatex "<<outname<<endl;


  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
