// table_yields: Generates a LaTeX file with a cutflow table for RA4

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <unistd.h> // getopt in Macs
#include <iomanip>  // setw

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"
#include "RooStats/NumberCountingUtils.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_full.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;
namespace {
  TString luminosity = "2.246";
}

void printTable(vector<sfeats> Samples, tfeats table, vector<vector<double> > yields, vector<vector<double> > w2, 
    vector<vector<double> > entries, size_t ini);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);


  //// Defining samples, i.e. columns in the table
  TString folder="/cms27r0/babymaker/2016_04_29/mc/merged_met100nb2nj4nl0/";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  folder = "/net/cms2"+folder;

  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*");
  // s_tt.push_back(folder+"*_TTJets_HT*");
  // vector<TString> s_t1t;
  // s_t1t.push_back(folder+"*T1tttt*1500_*");
  // vector<TString> s_t1tc;
  // s_t1tc.push_back(folder+"*T1tttt*1200_*");
  // vector<TString> s_other;
  // s_other.push_back(folder+"*DYJetsToLL*");
  // s_other.push_back(folder+"*_ZJet*");
  // s_other.push_back(folder+"*_WWTo*");
  // s_other.push_back(folder+"*ttHJetTobb*");
  // s_other.push_back(folder+"*_TTTT*");
  // s_other.push_back(folder+"*_WZ*.root");

  // vector<TString> s_qcd;
  // s_qcd.push_back(folder+"*_QCD_HT*");
  // s_qcd.push_back(folder+"*_TTJets_TuneCUET*");
  // vector<TString> s_wjets;
  // s_wjets.push_back(folder+"*_WJetsToLNu*");
  // vector<TString> s_ttv;
  // s_ttv.push_back(folder+"*_TTWJets*");
  // s_ttv.push_back(folder+"*_TTZTo*");
  // s_ttv.push_back(folder+"*_TTG*");
  // vector<TString> s_single;
  // s_single.push_back(folder+"*_ST_*");
 
  vector<sfeats> Samples; 
  // Samples.push_back(sfeats(s_other, "Other", 1001,1,"stitch"));
  // Samples.push_back(sfeats(s_qcd, "QCD", 1002, 1,"ntruleps==0"));
  // Samples.push_back(sfeats(s_ttv, "$t\\bar{t}V$", 1002));
  // Samples.push_back(sfeats(s_single, "Single $t$", 1005));
  // Samples.push_back(sfeats(s_wjets, "W+jets", 1004));
  Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ (1$\\ell$)", 1000,1, "ntruleps==1"));
  // Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ ($2\\ell$)", 1006,1,"ntruleps==2"));
  // Samples.push_back(sfeats(s_t1t, "T1tttt NC", 2));
  // Samples.push_back(sfeats(s_t1tc, "T1tttt C", 2,2));

  //// tables has a vector of the tables you want to print
  vector<tfeats> tables;
  TString baseline_s("stitch&&pass"); 

  //////////// Standard cutflow ////////////
  // Pushing first table and adding rows
  tables.push_back(tfeats("1", "an"));
  tables.back().add("No selection", "1");
  // tables.back().add("$1\\ell$", "nvleps==0");
  tables.back().add("$H_T>500$ GeV", "ht>500");
  // tables.back().add("MET$>200$ GeV", "met>200&&ht>500&&nleps==1");
  // tables.back().add("$N_{\\rm jets}\\geq6$", "njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("$N_{\\rm b}\\geq1$", "nbm>=1&&njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("$M_J>250$ GeV", "mj>250&&nbm>=1&&njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("$m_T>140$ GeV", "mt>140&&mj>250&&nbm>=1&&njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("$M_J>400$ GeV", "mt>140&&mj>400&&nbm>=1&&njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("$N_{\\rm b}\\geq2$", "mt>140&&mj>400&&mj>250&&nbm>=2&&njets>=6&&met>200&&ht>500&&nleps==1");
  // tables.back().add("MET$>400$ GeV", "mt>140&&mj>400&&nbm>=2&&njets>=6&&met>400&&ht>500&&nleps==1");
  // tables.back().add("$N_{\\rm jets}\\geq9$", "mt>140&&mj>400&&nbm>=2&&njets>=9&&met>400&&ht>500&&nleps==1");


  /////////////////////////////  No more changes needed down here to add tables ///////////////////////

  //// Concatenating cuts of all table rows in bincuts
  bcut baseline(baseline_s);
  vector<bcut> bincuts;
  for(size_t itab(0); itab < tables.size(); itab++){
    for(size_t icut(0); icut < tables[itab].tcuts.size(); icut++){
      bincuts.push_back(bcut(tables[itab].cuts+"&&"+tables[itab].tcuts[icut]));
    }
    tables[itab].cuts = "ht>500&&met>200&&"+baseline_s + "&&" + tables[itab].cuts;
  }
  //// Calculating yields per sample, all bins from all tables at a time
  vector<vector<double> > yields, w2, entries;
  vector<int> repeat_sam(Samples.size(), -1); // used when the same sample appears multiple times in Samples
  for(size_t sam(0); sam < Samples.size(); sam++){
    vector<bcut> samcuts;
    //// Adding specific sample cut to bin cuts
    for(size_t bin(0); bin < bincuts.size(); bin++)
      samcuts.push_back(bcut(bincuts[bin].cuts_+"&&"+Samples[sam].cut));
    for(size_t sam2(sam+1); sam2 < Samples.size(); sam2++){
      //// If 2 samples are the same, the bincuts are concatened so that the yields are found in one go
      if(Samples[sam].file == Samples[sam2].file) {
        repeat_sam[sam2] = sam;
        for(size_t bin(0); bin < bincuts.size(); bin++)
          samcuts.push_back(bcut(bincuts[bin].cuts_+"&&"+Samples[sam2].cut)); 
      }
    } // Loop over future samples
    if(repeat_sam[sam]==-1){
      //// Creating baby with all samples pushed to Samples[sam]
      baby_full baby(Samples[sam].file[0]);
      for(size_t file(1); file < Samples[sam].file.size(); file++) baby.Add(Samples[sam].file[file]);
      //// Finding yields for all bins (samcuts = bincuts&&Sample[sam].cut) for this sample
      yields.push_back(vector<double>());
      w2.push_back(vector<double>());
      time(&endtime); 
      cout<<setw(3)<<difftime(endtime, begtime)<<" seconds passed. Finding yields for \""<<Samples[sam].label
    <<"\" with "<<baby.GetEntries()<<" entries"<<endl;
      entries.push_back(getYields(baby, baseline, samcuts, yields.back(), w2.back(), luminosity.Atof()));
    } else {
      //// If the yields were calculated earlier, put them in the correct vector and erase them from the previous sample
      size_t ncuts(bincuts.size());
      entries.push_back(vector<double>(&entries[repeat_sam[sam]][ncuts], &entries[repeat_sam[sam]][2*ncuts+1]));
      entries[repeat_sam[sam]].erase(entries[repeat_sam[sam]].begin()+ncuts, entries[repeat_sam[sam]].begin()+ncuts*2);
      yields.push_back(vector<double>(&yields[repeat_sam[sam]][ncuts], &yields[repeat_sam[sam]][2*ncuts+1]));
      yields[repeat_sam[sam]].erase(yields[repeat_sam[sam]].begin()+ncuts, yields[repeat_sam[sam]].begin()+ncuts*2);
      w2.push_back(vector<double>(&w2[repeat_sam[sam]][ncuts], &w2[repeat_sam[sam]][2*ncuts+1]));
      w2[repeat_sam[sam]].erase(w2[repeat_sam[sam]].begin()+ncuts, w2[repeat_sam[sam]].begin()+ncuts*2);
    }
  } // Loop over samples

  //// Printing each table. ini keeps track of the index in yields where the specific table begins
  size_t ini(0);
  for(size_t itab(0); itab < tables.size(); itab++) {
    printTable(Samples, tables[itab], yields, w2, entries, ini);
    ini +=  tables[itab].tcuts.size();
  }
  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void printTable(vector<sfeats> Samples, tfeats table, vector<vector<double> > yields, vector<vector<double> > w2, 
    vector<vector<double> > entries, size_t ini) {
  int nsig(0), digits(1);
  for(unsigned sam(0); sam < Samples.size(); sam++) if(Samples[sam].isSig) nsig++;

  bool do_uncert(false);
  vector<double> av_w2(Samples.size()+1,0);

  TString outname = "txt/table_cutflow_"+table.tag+".tex";
  ofstream out(outname);

  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating}\n";
  out << "\\usepackage[landscape]{geometry}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\n\\begin{tabular}[tbp!]{ l | ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++) out << "r";
  out<<" | r ";
  for(int sam(0); sam < nsig; sam++) out<<"| r ";
  out<<"}\\hline\\hline\n";
  out << " \\multicolumn{1}{c|}{${\\cal L} = "<<luminosity<<"$ fb$^{-1}$} ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++)
    out << " & "<<Samples[sam].label;
  out<< " & SM bkg. ";
  for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++)
    out << " & "<<Samples[sam].label;
  out << "\\\\ \\hline \n ";
  for(size_t icut(0); icut < table.tcuts.size(); icut++){
    if(table.options[icut]=="=") do_uncert = true;
    if(icut==table.tcuts.size()-3) digits = 2;
    for(int ind(0); ind < table.options[icut].CountChar('='); ind++) out << " \\hline ";
    out<<table.texnames[icut];
    double bkg(0), ebkg(0);
    for(unsigned sam(0); sam < Samples.size()-nsig; sam++) {
      double val(yields[sam][ini+icut]), errval(sqrt(w2[sam][ini+icut]));
      if(w2[sam][ini+icut]>0) av_w2[sam] = w2[sam][ini+icut]/entries[sam][ini+icut];
      else errval = sqrt(av_w2[sam]);
      bkg += val;
      ebkg += pow(errval, 2);
      if(!do_uncert) out <<" & "<< RoundNumber(val,digits);
      else out <<" & $"<< RoundNumber(val,digits)<<"\\pm"<<RoundNumber(errval,digits)<<"$";
    } // Loop over background samples
    if(!do_uncert) out <<" & "<< RoundNumber(bkg,digits);
    else out <<" & $"<< RoundNumber(bkg,digits)<<"\\pm"<<RoundNumber(sqrt(ebkg),digits)<<"$";
    for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++){
      double val(yields[sam][ini+icut]), errval(sqrt(w2[sam][ini+icut]));
      if(!do_uncert) out <<" & "<< RoundNumber(val,digits);
      else out <<" & $"<< RoundNumber(val,digits)<<"\\pm"<<RoundNumber(errval,digits)<<"$";
    }
    out<<" \\\\ ";
    for(int ind(0); ind < table.options[icut].CountChar('-'); ind++) out << " \\hline";
    out<<endl;
  }// Loop over table cuts



  out << "\\hline\\multicolumn{1}{c|}{} ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++)
    out << " & "<<Samples[sam].label;
  out<< " & SM bkg. ";
  for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++)
    out << " & "<<Samples[sam].label;
  out << "\\\\ \n ";

  out<< "\\hline\\hline\n\\end{tabular}"<<endl<<endl;
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out.close();
  cout<<" pdflatex "<<outname<<endl;
}


