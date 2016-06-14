#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>


#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_full.hpp"
#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  double lumi(20.);
  TString nb_bins = "TTML";
  vector<TString> metlows = {"100","150","200","250"}; //low edge for each met bin
  //options: ttbar, all_but_qcd, all
  TString sample_set = "ttbar";
  bool loose = false; 
  TString hig_dm = "40";
  TString title_style("RA4"); //CMSPaper
}

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2016_04_29/mc/merged_higloose/";
  if (loose) folder.ReplaceAll("merged_higloose","merged_met100nb2nj4nl0");

  ////// Creating babies
  baby_full tt(folder+"*TTJets*Lept*");
  tt.Add(folder+"*TTJets*HT*");
  if (sample_set=="all_but_qcd" || sample_set=="all"){
    tt.Add(folder+"*_WJetsToLNu*");
    tt.Add(folder+"*ZJetsToNuNu_HT*");
    tt.Add(folder+"*ST_tW*");
    tt.Add(folder+"*DYJetsToLL_M-50_HT*");
    tt.Add(folder+"*TTTT*");
    tt.Add(folder+"*_TTWJets*");
    tt.Add(folder+"*_WH_*");
    tt.Add(folder+"*_ttHJet*");
    tt.Add(folder+"*_ZH_*");
    tt.Add(folder+"*_ST_*channel*");
    tt.Add(folder+"*_TTGJets*");
    tt.Add(folder+"*_TTZTo*");
  }
  if (sample_set=="all") tt.Add(folder+"*_QCD_HT*");
  ////// Defining cuts

  TString tag = "";
  TString baseline("pass&&stitch&&njets>=4&&njets<=5&&!low_dphi&&hig_drmax<2.2&&ntks==0");  
  if (loose) {
    baseline = "pass&&stitch&&njets>=4&&njets<=5&&hig_drmax<2.5";  
    tag +="_loose";    
  }

  if (hig_dm!="40") tag += "_higdm-"+hig_dm;
  TString sigcut="hig_am>100&&hig_am<140&&hig_dm<"+hig_dm;
  TString sbdcut="sbd";
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
 
  vector<TString> metbins, metnames;
  for (size_t ilow(0); ilow<metlows.size(); ilow++){
    if (ilow==metlows.size()-1) {
      metbins.push_back("met>"+metlows[ilow]);
      metnames.push_back("E_{T}^{miss} > "+metlows[ilow]);
    } else{
      metbins.push_back("met>"+metlows[ilow]+"&&met<="+metlows[ilow+1]);
      metnames.push_back(metlows[ilow]+" < E_{T}^{miss} #leq "+metlows[ilow+1]);
    }
  }

  // changing this is likely to break everything...
  vector<vector<TString> > abcdtypes; vector<TString> abcdnames;
  abcdtypes.push_back(vector<TString>({sigcut+"&&"+cut3b, sbdcut+"&&"+cut3b, sigcut+"&&"+cut2b, sbdcut+"&&"+cut2b})); 
  abcdnames.push_back("3b");
  abcdtypes.push_back(vector<TString>({sigcut+"&&"+cut4b, sbdcut+"&&"+cut4b, sigcut+"&&"+cut2b, sbdcut+"&&"+cut2b}));
  abcdnames.push_back("4b");
  vector<TString> nbnames = vector<TString>{"2b","3b","4b"};

  ////// Combining cuts
  vector<bcut > bincuts;
  vector<bool> do_sbd;
  for(size_t ind(0); ind<metbins.size(); ind++){
    for (size_t itype(0); itype<abcdtypes.size();itype++){
      for(size_t obs(0); obs < abcdtypes[itype].size(); obs++){
        TString totcut(abcdtypes[itype][obs]+"&&"+metbins[ind]);
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

  vector<float> pow_kappa = vector<float>({1,-1,-1,1});
  float mSigma, pSigma, kappa;
  // kappa indices: nabcdtypes*imet + iabcd_type, imet, valu+-err
  vector<vector<float> > kappas;
  // sig/sbd ratio indices: nb*imet +inb, valu+-err
  vector<vector<float> > rss;
  for(size_t imet(0); imet<metbins.size(); imet++){
    for(size_t itype(0); itype<abcdtypes.size(); itype++){
      vector<TString> abcdcuts = abcdtypes[itype];
      vector<vector<float> > entries;
      vector<vector<float> > weights;
      for(size_t obs(0); obs < abcdcuts.size(); obs++){
        size_t index(abcdtypes.size()*abcdcuts.size()*imet+itype*abcdcuts.size()+obs);
        entries.push_back(vector<float>());
        weights.push_back(vector<float>());
        entries.back().push_back(pow(mcyield[index],2)/mcw2[index]);
        weights.back().push_back(mcw2[index]/mcyield[index]);

      } // Loop over observables for MC
      kappa = calcKappa(entries, weights, pow_kappa, mSigma, pSigma); if(mSigma<0) mSigma = 0;
      kappas.push_back(vector<float>({kappa, mSigma<0 ? 0:mSigma, pSigma}));

      //piggy-back rss
      float rat(0); 
      vector<float> pow_inb;
      if (itype==0){ //calculate sig/sbd for 2b and 3b
        pow_inb = vector<float>({0.,0.,1.,-1.});
        rat = calcKappa(entries, weights, pow_inb, mSigma, pSigma);
        rss.push_back(vector<float>({rat, mSigma<0 ? 0:mSigma, pSigma}));
        pow_inb = vector<float>({1.,-1.,0.,0.});
        rat = calcKappa(entries, weights, pow_inb, mSigma, pSigma);
        rss.push_back(vector<float>({rat, mSigma<0 ? 0:mSigma, pSigma}));
      } else { // calculate sig/sbd for 4b
        pow_inb = vector<float>({1.,-1.,0.,0.});
        rat = calcKappa(entries, weights, pow_inb, mSigma, pSigma);
        rss.push_back(vector<float>({rat, mSigma<0 ? 0:mSigma, pSigma}));
      }
    } // Loop over signal bins
  }

  ///// Printing table
  TString outname = "txt/table_kappa_"+nb_bins+tag+".tex";
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
  cout<<"---------------------------------------------------"<<endl;
  cout << setw(20)<< "Bin"<<setw(10)<< "kappa"<< setw(10)<< "-"<< setw(10)<< "+"<<endl;
  cout<<"---------------------------------------------------"<<endl;
  for(size_t imet(0); imet<metbins.size(); imet++){
    out << " & \\multicolumn{3}{c}{"<<cuts2tex(metbins[imet])<<"}  \\\\ \\hline\n";
    for(size_t itype(0); itype<abcdtypes.size(); itype++){
      size_t ind(abcdtypes.size()*imet+itype);
      out<<abcdnames[itype]<<" & "<<RoundNumber(kappas[ind][0],2) <<" & "<<kappas[ind][1] <<" & "<<kappas[ind][2];
        out<<" \\\\"<<endl;
      out<<"\\hline"<<endl;
      cout<<setw(20)<<metnames[imet].Copy().ReplaceAll("#leq","<").ReplaceAll("E_{T}^{miss}","MET")+","+abcdnames[itype]
                    <<setw(10)<<RoundNumber(kappas[ind][0],2)
                    <<setw(10)<<RoundNumber(kappas[ind][1],2)
                    <<setw(10)<<RoundNumber(kappas[ind][2],2)<<endl;
    }
    cout<<"---------------------------------------------------"<<endl;
  }
  out<< "\\hline\\hline\n\\end{tabular}"<<endl<<endl;
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out.close();
  TString pdfname(outname); 
  cout<<" pdflatex "<<outname<<endl<<endl;

  //    Print also ratios
  //----------------------------
  cout<<"---------------------------------------------------"<<endl;
  cout << setw(20)<< "# b's"<<setw(10)<< "Sig/Sbd"<< setw(10)<< "-"<< setw(10)<< "+"<<endl;
  cout<<"---------------------------------------------------"<<endl;
  for(size_t imet(0); imet<metbins.size(); imet++){
    for(size_t inb(0); inb<nbnames.size(); inb++){
      size_t ind(nbnames.size()*imet+inb);
      cout<<setw(20)<<metnames[imet].Copy().ReplaceAll("#leq","<").ReplaceAll("E_{T}^{miss}","MET")+","+nbnames[inb]
                    <<setw(10)<<RoundNumber(rss[ind][0],2)
                    <<setw(10)<<RoundNumber(rss[ind][1],2)
                    <<setw(10)<<RoundNumber(rss[ind][2],2)<<endl;
    }
    cout<<"---------------------------------------------------"<<endl;
  }

  //       Plotting
  //------------------------
  float minh(0), maxh(metbins.size());
  TH1D histo("histo",nb_bins+", "+cuts2title(baseline).ReplaceAll("stitch, ","").ReplaceAll("pass, ",""),metbins.size(), minh, maxh);
  // if (title_style=="CMSPaper") histo.SetTitle("");
  for(unsigned imet(0); imet<metbins.size(); imet++){
    TString mettitle = cuts2title(metbins[imet]);
    mettitle.ReplaceAll("MET","E_{T}^{miss}");
    histo.GetXaxis()->SetBinLabel(1+imet, metnames[imet]);
  }
    
  TString stylename = "RA4";
  styles style(stylename);
  style.setDefaultStyle();
  float max_axis(1.5), min_axis(0.5);
  unsigned ntypes(abcdtypes.size());
  for(size_t imet(0); imet<metbins.size(); imet++){
    for(size_t itype(0); itype<ntypes; itype++){
      size_t ind(abcdtypes.size()*imet+itype);
      if(kappas[ind][0] > max_axis && kappas[ind][0]-kappas[ind][1] < max_axis) kappas[ind][0] = max_axis;      
    }
  }

  TCanvas can;
  TLatex cmslabel; 
  if (title_style=="CMSPaper"){
    cmslabel.SetNDC(kTRUE);
    cmslabel.SetTextAlign(11);
    cmslabel.DrawLatex(0.175,0.94,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");  
    cmslabel.SetTextAlign(31);
    cmslabel.DrawLatex(0.953,0.94,"13 TeV");  
  }

  TString ytitle("#left(#frac{SIG}{SBD}#right)_{nb} / #left(#frac{SIG}{SBD}#right)_{2b}"); 
  histo.GetYaxis()->CenterTitle(true);
  histo.GetYaxis()->SetTitleOffset(1.5);
  histo.GetYaxis()->SetTitleSize(0.06);
  histo.SetYTitle(ytitle);
  histo.SetMaximum(max_axis);
  histo.SetMinimum(min_axis);
  style.moveYAxisLabel(&histo, max_axis, false);
  histo.GetXaxis()->SetLabelOffset(0.007);
  double legX(style.PadLeftMargin+0.03), legY(0.895), legSingle = 0.052;
  double legW = 0.3, legH = legSingle*ntypes;
  legH = legSingle*((nbnames.size()+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(style.nFont); 
  leg.SetNColumns(2);

  //     Plot kappa
  //---------------------------
  histo.Draw();
  vector<TGraphAsymmErrors*> graph;
  int colors[] = {kGreen+1, kAzure-3, kRed-4}, styles[] = {20, 21, 22, 23};
  for(unsigned itype(0); itype<ntypes; itype++){
    vector<float> vy, veyl, veyh, vx, vexl, vexh;
    for(size_t imet(0); imet<metbins.size(); imet++){
      size_t ind(ntypes*imet+itype);
      vx.push_back(float(imet)+float(itype+1)/float(ntypes+1)); vexl.push_back(0.); vexh.push_back(0.);
      vy.push_back(kappas[ind][0]); veyl.push_back(kappas[ind][1]); veyh.push_back(kappas[ind][2]);
    }
    graph.push_back(new TGraphAsymmErrors(metbins.size(), &(vx[0]), &(vy[0]), &(vexl[0]), &(vexh[0]), &(veyl[0]), &(veyh[0])));
    graph.back()->SetMarkerStyle(styles[itype]); graph.back()->SetMarkerSize(1.65);  graph.back()->SetLineWidth(2);
    graph.back()->SetMarkerColor(colors[itype+1]); graph.back()->SetLineColor(colors[itype+1]);
    graph.back()->Draw("p same");   
    leg.AddEntry(graph.back(), abcdnames[itype], "p");
  }

  leg.Draw();
  // draw a line at 1
  TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
  line.DrawLine(minh, 1, maxh, 1);

  TString pname = "plots/kappa_"+nb_bins+"_"+sample_set+tag+".pdf";
  can.SaveAs(pname);
  cout<<endl<<" open "<<pname<<endl<<endl;

  //     Plot SIG/SBD
  //---------------------------
  unsigned nnb(nbnames.size());
  //recalculate maximum
  min_axis = 0.15; max_axis = 0.55; 
  for(size_t imet(0); imet<metbins.size(); imet++){
    for(size_t inb(0); inb<nnb; inb++){
      size_t ind(nnb*imet+inb);
      if(rss[ind][0] > max_axis && rss[ind][0]-rss[ind][1] < max_axis) rss[ind][0] = max_axis;
      if(rss[ind][0] < min_axis && rss[ind][0]+rss[ind][2] > min_axis) rss[ind][0] = min_axis;
    }
  }
  histo.SetMaximum(max_axis);
  histo.SetMinimum(min_axis);
  histo.SetYTitle("SIG/SBD ratio");
  histo.Draw(); 
  leg.Clear();

  for(unsigned inb(0); inb<nnb; inb++){
    vector<float> vy, veyl, veyh, vx, vexl, vexh;
    for(size_t imet(0); imet<metbins.size(); imet++){
      size_t ind(nnb*imet+inb);
      vx.push_back(float(imet)+float(inb+1)/float(nnb+1)); vexl.push_back(0.); vexh.push_back(0.);
      vy.push_back(rss[ind][0]); veyl.push_back(rss[ind][1]); veyh.push_back(rss[ind][2]);
    }
    graph.push_back(new TGraphAsymmErrors(metbins.size(), &(vx[0]), &(vy[0]), &(vexl[0]), &(vexh[0]), &(veyl[0]), &(veyh[0])));
    graph.back()->SetMarkerStyle(styles[inb]); graph.back()->SetMarkerSize(1.65);  graph.back()->SetLineWidth(2);
    graph.back()->SetMarkerColor(colors[inb]); graph.back()->SetLineColor(colors[inb]);
    graph.back()->Draw("p same");   
    leg.AddEntry(graph.back(), nbnames[inb], "p");
  }

  leg.Draw();
  line.DrawLine(minh, 1, maxh, 1);

  pname = "plots/rss_"+nb_bins+"_"+sample_set+tag+".pdf";
  can.SaveAs(pname);
  cout<<endl<<" open "<<pname<<endl<<endl;

  for (size_t igr(0); igr<graph.size(); igr++) delete graph[igr];

  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
