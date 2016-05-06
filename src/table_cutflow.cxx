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
  TString luminosity = "20.";
  TString nb_bins = "TML";
}

void printTable(vector<sfeats> Samples, tfeats table, vector<vector<double> > yields, vector<vector<double> > w2, 
    vector<vector<double> > entries, size_t ini);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);


  //// Defining samples, i.e. columns in the table
  TString folder="/cms27r0/babymaker/2016_04_29/mc/merged_met100nb2nj4nl0/";
  // TString folder="/cms27r0/babymaker/2016_04_29/mc/merged_higloose/";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  folder = "/net/cms2"+folder;

  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*Lept*");
  s_tt.push_back(folder+"*_TTJets_HT*");
  
  vector<TString> s_tchi400;
  s_tchi400.push_back(folder+"*-TChiHH_mChi-400*");
  // vector<TString> s_tchi200;
  // s_tchi200.push_back(folder+"*-TChiHH_mChi-200*");
  vector<TString> s_other;
  s_other.push_back(folder+"*DYJetsToLL*");
  s_other.push_back(folder+"*_WWTo*");
  s_other.push_back(folder+"*_TTTT*");
  s_other.push_back(folder+"*_WZ*.root");

  vector<TString> s_qcd;
  s_qcd.push_back(folder+"*_QCD_HT*");
  // s_qcd.push_back(folder+"*_TTJets_TuneCUET*");
  vector<TString> s_vjets;
  s_vjets.push_back(folder+"*_WJetsToLNu*");
  s_vjets.push_back(folder+"*_ZJetsToNuNu*");
  vector<TString> s_tx;
  s_tx.push_back(folder+"*_TTWJets*");
  s_tx.push_back(folder+"*_TTZTo*");
  s_tx.push_back(folder+"*_TTG*");
  s_tx.push_back(folder+"*ttHJetTobb*");
  s_tx.push_back(folder+"*_ST_*");

 
  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_other, "Other", 1001,1,"stitch"));
  Samples.push_back(sfeats(s_qcd, "QCD", 1002, 1,"ntruleps==0"));
  Samples.push_back(sfeats(s_tx, "$t\\bar{t}X$, single $t$", 1002));
  // Samples.push_back(sfeats(s_single, "Single $t$", 1005));
  Samples.push_back(sfeats(s_vjets, "W/Z+jets", 1004));
  // Samples.push_back(sfeats(s_zjets, "Z+jets", 1004));
  Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ (1$\\ell$)", 1000,1, "stitch&&ntruleps==1"));
  Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ ($2\\ell$)", 1006,1,"stitch&&ntruleps==2"));
  Samples.push_back(sfeats(s_tchi400, "TChiHH 400", 2));Samples.back().isSig = true;
  // Samples.push_back(sfeats(s_tchi200, "TChiHH 200", 2,2));

  //// tables has a vector of the tables you want to print
  vector<tfeats> tables;
  TString baseline_s("pass"); 

  TString cut2b="nbt==2&&nbm==2", cut3b="nbt>=2&&nbm==3&&nbl==3", cut4b="nbt>=2&&nbm>=3&&nbl>=4";
  if(nb_bins=="TTL"){
    cut2b = "nbt==2";
    cut3b = "nbt==3&&nbl==3";
    cut4b = "nbt==3&&nbl>=4";
  }
  if(nb_bins=="TMM"){
    cut2b = "nbt==2&&nbm==2";
    cut3b = "nbm==3";
    cut4b = "nbm>=4";
  }
  TString met0cut="met>250&&met<=300&&", met1cut="met>300&&";
  TString sigcut="hig_am>100&&hig_am<140&&hig_dm<40";

  //////////// Standard cutflow ////////////
  TString skim("met>150 && nbm>=2 && njets>=4 && nvleps==0&&"); 
  tables.push_back(tfeats("1", "an"));
  tables.back().add("MET $>150$, 2M b-tags, 4 jets, 0$\\ell$", skim);
  tables.back().add("Iso track veto", skim+" ntks==0");
  tables.back().add("$\\Delta\\phi_{min}$", skim+"ntks==0&&!low_dphi");
  tables.back().add("$N_{\\rm jets}\\leq 5$", skim+"ntks==0&&!low_dphi && njets<=5","-");

  TString preseln(skim+"ntks==0&&!low_dphi && njets<=5 &&");
  tables.back().add("$\\Delta m_{jj}<40$ ", preseln+"hig_dm<40");
  tables.back().add("$\\left< m_{jj} \\right> \\in(100,140)$ ", preseln+"hig_dm<40 && hig_am>100 && hig_am<140");
  tables.back().add("$\\Delta R_{max} < 2.2$ ", preseln+"hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2","-");
  TString higcuts("hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2 &&");
  tables.back().add("$MET \\geq 250$", preseln+ higcuts +"met>250");
  tables.back().add("$N_{b,T}\\geq2$", preseln+ higcuts +"met>250&&nbt>=2","-");

  tables.back().add("N-1, Iso track veto", skim + "met>250 && nbt>=2 &&" + "!low_dphi && njets<=5 && hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2","=");
  tables.back().add("N-1, $\\Delta\\phi_{min}$", skim + "met>250 && nbt>=2 &&" + "ntks==0 && njets<=5 && hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2");
  tables.back().add("N-1, $N_{\\rm jets}\\leq 5$", skim + "met>250 && nbt>=2 &&" + "ntks==0 && !low_dphi && hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2");
  tables.back().add("N-1, $\\Delta m_{jj}<40$", skim + "met>250 && nbt>=2 &&" + "ntks==0 && !low_dphi && njets<=5 && hig_am>100 && hig_am<140 && hig_drmax<2.2");
  tables.back().add("N-1, $\\left< m_{jj} \\right> \\in(100,140)$", skim + "met>250 && nbt>=2 &&" + "ntks==0 && !low_dphi && njets<=5 && hig_dm<40 && hig_drmax<2.2");
  tables.back().add("N-1, $\\Delta R_{max} < 2.2$", skim + "met>250 && nbt>=2 &&" + "ntks==0 && !low_dphi && njets<=5 && hig_dm<40 && hig_am>100 && hig_am<140","--");
  // tables.back().add("N-1, Iso track veto", skim + "met>250 && nbt>=2 &&" + "ntks==0 && !low_dphi && njets<=5 && hig_dm<40 && hig_am>100 && hig_am<140 && hig_drmax<2.2");

  //N-1 rows

  tables.back().add("$MET \\in (250,300), N_{b,T}\\geq2$", preseln+ higcuts +met0cut + "nbt>=2","-");
  tables.back().add("$MET \\in (250,300), N_{b,T}=2, N_{b,M}=2$", preseln+ higcuts + met0cut + cut2b);
  tables.back().add("$MET \\in (250,300), N_{b,T}\\geq2, N_{b,M}=3, N_{b,L}=3$", preseln+ higcuts + met0cut + cut3b);
  tables.back().add("$MET \\in (250,300), N_{b,T}\\geq2, N_{b,M}\\geq3, N_{b,L}\\geq4$", preseln+ higcuts + met0cut + cut4b,"--");
  tables.back().add("$MET>300, N_{b,T}\\geq2$", preseln+ higcuts +met1cut + "nbt>=2","-");
  tables.back().add("$MET > 300, N_{b,T}=2, N_{b,M}=2$", preseln+ higcuts + met1cut + cut2b);
  tables.back().add("$MET > 300, N_{b,T}\\geq2, N_{b,M}=3, N_{b,L}=3$", preseln+ higcuts + met1cut + cut3b);
  tables.back().add("$MET > 300, N_{b,T}\\geq2, N_{b,M}\\geq3, N_{b,L}\\geq4$", preseln+ higcuts + met1cut + cut4b);

  /////////////////////////////  No more changes needed down here to add tables ///////////////////////

  //// Concatenating cuts of all table rows in bincuts
  bcut baseline(baseline_s);
  vector<bcut> bincuts;
  for(size_t itab(0); itab < tables.size(); itab++){
    for(size_t icut(0); icut < tables[itab].tcuts.size(); icut++){
      bincuts.push_back(bcut(tables[itab].cuts+"&&"+tables[itab].tcuts[icut]));
    }
    //tables[itab].cuts = "ht>500&&met>200&&"+baseline_s + "&&" + tables[itab].cuts;
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
  int nsig(0), digits(2);
  for(unsigned sam(0); sam < Samples.size(); sam++) if(Samples[sam].isSig) nsig++;

  bool do_uncert(false);
  vector<double> av_w2(Samples.size()+1,0);

  TString outname = "txt/table_cutflow_"+table.tag+".tex";
  ofstream out(outname);

  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating,color}\n";
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
    if(table.options[icut].Contains("=")) do_uncert = true;
    // if(icut==table.tcuts.size()-3) digits = 2;
    for(int ind(0); ind < table.options[icut].CountChar('='); ind++) out << " \\hline ";
    TString color_on(""), color_off("");
    if (table.options[icut].Contains("b")) { color_on="\\textcolor{blue}{"; color_off="}"; }
    else if (table.options[icut].Contains("r")) { color_on="\\textcolor{red}{"; color_off="}"; }
    if (table.options[icut].Contains("e")) { color_on+="\\textbf{"; color_off+="}"; }
    out<<color_on<<table.texnames[icut]<<color_off;
    double bkg(0), ebkg(0);
    for(unsigned sam(0); sam < Samples.size()-nsig; sam++) {
      double val(yields[sam][ini+icut]), errval(sqrt(w2[sam][ini+icut]));
      if(w2[sam][ini+icut]>0) av_w2[sam] = w2[sam][ini+icut]/entries[sam][ini+icut];
      else errval = sqrt(av_w2[sam]);
      bkg += val;
      ebkg += pow(errval, 2);
      if(!do_uncert) out <<" & "+color_on<< RoundNumber(val,digits)<<color_off;
      else out <<" & "+color_on+"$"<< RoundNumber(val,digits)<<"\\pm"<<RoundNumber(errval,digits)<<"$"+color_off;
    } // Loop over background samples
    if(!do_uncert) out <<" & "+color_on<< RoundNumber(bkg,digits)<<color_off;
    else out <<" & "+color_on+"$"<< RoundNumber(bkg,digits)<<"\\pm"<<RoundNumber(sqrt(ebkg),digits)<<"$"+color_off;
    for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++){
      double val(yields[sam][ini+icut]), errval(sqrt(w2[sam][ini+icut]));
      if(!do_uncert) out <<" & "+color_on<< RoundNumber(val,digits)<<color_off;
      else out <<" & "+color_on+"$"<< RoundNumber(val,digits)<<"\\pm"<<RoundNumber(errval,digits)<<"$"+color_off;
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


