
#include "bcut.hpp"

using namespace std;

void bcut::parseWeights(TString weights){
  weights.ReplaceAll(" ", "");
  cutTypes_.clear();
  fWeights_.clear();
  fvWeights_.clear();
  indWeights_.clear();
  TString wgt;
  int pos;
  do{
    pos = weights.Index("*");
    wgt = weights;
    if(pos != -1) wgt.Remove(pos, wgt.Length());
    parseWeight(wgt);
    weights.Remove(0, pos+1);
  } while(pos != -1);
}

void bcut::parseWeight(TString wgt){
  cutTypes_.push_back(kFloat);
  fWeights_.push_back(NULL);
  fvWeights_.push_back(NULL);
  indWeights_.push_back(-1);
  constWeights_.push_back(1.);

  if(wgt=="weight")		fWeights_.back() = &baby_base::weight;
  else if(wgt=="w_lumi")	fWeights_.back() = &baby_base::w_lumi;
  else if(wgt=="w_pu")	fWeights_.back() = &baby_base::w_pu;
  else if(wgt=="w_lep")	fWeights_.back() = &baby_base::w_lep;
  else if(wgt=="w_fs_lep")	fWeights_.back() = &baby_base::w_fs_lep;
  else if(wgt=="w_toppt")	fWeights_.back() = &baby_base::w_toppt;
  else if(wgt=="w_btag")	fWeights_.back() = &baby_base::w_btag;
  else if(wgt=="eff_trig")	fWeights_.back() = &baby_base::eff_trig;
  else if(wgt.Contains("[")){ // if weight is a vector element
    TString index_s(wgt);
    wgt.Remove(wgt.Index("["), wgt.Length());
    index_s.Remove(0, index_s.Index("[")+1);
    index_s.Remove(index_s.Index("]"), index_s.Length());
    indWeights_.back() = index_s.Atoi();
    cutTypes_.back() = kvFloat;
    if(wgt=="w_pdf")        fvWeights_.back() = &baby_base::w_pdf;
    else if(wgt=="sys_pdf") fvWeights_.back() = &baby_base::sys_pdf;
    else if(wgt=="sys_isr") fvWeights_.back() = &baby_base::sys_isr;
    else if(wgt=="sys_mur") fvWeights_.back() = &baby_base::sys_mur;
    else if(wgt=="sys_muf") fvWeights_.back() = &baby_base::sys_muf;
    else if(wgt=="sys_murf") fvWeights_.back() = &baby_base::sys_murf;
    else if(wgt=="sys_trig") fvWeights_.back() = &baby_base::sys_trig;
    else if(wgt=="sys_lep") fvWeights_.back() = &baby_base::sys_lep;
    else if(wgt=="sys_fs_lep") fvWeights_.back() = &baby_base::sys_fs_lep;
    else if(wgt=="sys_bctag") fvWeights_.back() = &baby_base::sys_bctag;
    else if(wgt=="sys_fs_bctag") fvWeights_.back() = &baby_base::sys_fs_bctag;
    else if(wgt=="sys_udsgtag") fvWeights_.back() = &baby_base::sys_udsgtag;
    else if(wgt=="sys_fs_udsgtag") fvWeights_.back() = &baby_base::sys_fs_udsgtag;
    else {
      cout<<"Weight \""<<wgt<<" not defined. Add it to bcut::parseWeight in bcut.cpp"<<endl;
      exit(0);
    }
  }else if(wgt.Atof()>0) {
    constWeights_.back() = wgt.Atof();
    cutTypes_.back() = kConst;
  } else {
    cout<<"Weight \""<<wgt<<" not defined. Add it to bcut::parseWeight  in bcut.cpp"<<endl;
    exit(0);
  }   
}

float bcut::weight(baby_base *baby){
  float result(1.);
  for(size_t ind(0); ind < cutTypes_.size(); ind++){
    if(cutTypes_[ind] == kFloat) result *= (baby->*fWeights_[ind])();
    else if(cutTypes_[ind] == kvFloat) result *= (baby->*fvWeights_[ind])()[indWeights_[ind]];
    else result *= constWeights_[ind];
  } // Loop over individual weights

  return result;
}

void bcut::parseCuts(TString cuts){
  vcuts_.clear(); // Starting the parsing every time
  cuts.ReplaceAll(" ", "");
  TString cut;
  int pos;
  do{
    pos = cuts.Index("&&");
    cut = cuts;
    if(pos != -1) cut.Remove(pos, cut.Length());
    vcuts_.push_back(onecut(cut));
    cuts.Remove(0, pos+2);
  } while(pos != -1);
}

bool bcut::pass(baby_base *baby){
  for(size_t ind(0); ind < vcuts_.size(); ind++)
    if(!vcuts_[ind].pass(baby)) return false;
  return true;
}

void bcut::operator+=(TString &cut){
  if(cut!="") cuts_ += ("&&"+cut);
  parseCuts(cuts_);
}

bcut bcut::operator+(bcut &ibcut){
  return bcut(cuts_+"&&"+ibcut.cuts_);
}

bcut::bcut(TString cuts, TString weights):
  cuts_(cuts),
  weights_(weights){
  parseCuts(cuts);
  parseWeights(weights);
}

bcut::~bcut(){
}

onecut::onecut(TString cut):
  cut_(cut){
  parseCut(cut);
}

onecut::~onecut(){
}

bool onecut::pass(baby_base *baby){
  int vali; float valf;
  switch(cutType_){
  default:
  case kAlwaysTrue:
    return true;
  case kAlwaysFalse:
    return false;
  case kBool:
    return (baby->*bb_)();
  case kvBool:
    return (baby->*bvb_)()[ivector_];
  case kFloat:
  case kvFloat:
    if(cutType_ == kFloat) valf = (baby->*bf_)();
    if(cutType_ == kvFloat) valf = (baby->*bvf_)()[ivector_];
    switch(compType_){
    case kNotEqual:     return valf != cutf_;
    case kLess:         return valf <  cutf_;
    case kLessEqual:    return valf <= cutf_;
    case kEqual:        return valf == cutf_;
    case kGreaterEqual: return valf >= cutf_;
    default:
    case kGreater:      return valf >  cutf_;
    }
  case kInt:
  case kvInt:
    if(cutType_ == kInt) vali = (baby->*bi_)();
    if(cutType_ == kvInt) vali = (baby->*bvi_)()[ivector_];
    switch(compType_){
    case kNotEqual:     return vali != cuti_;
    case kLess:         return vali <  cuti_;
    case kLessEqual:    return vali <= cuti_;
    case kEqual:        return vali == cuti_;
    case kGreaterEqual: return vali >= cuti_;
    default:
    case kGreater:      return vali >  cuti_;
    }
  }    

  return false;
}

void onecut::parseCut(TString cut){
  if(cut=="" || cut=="1"){
    cutType_ = kAlwaysTrue;
    return;
  }
  if(cut.Contains("&") || cut.Contains("|") || cut.Contains("(") || cut.Contains(")")){
    cout<<"Cut "<<cut<<" not fundamental. Exiting"<<endl;
    cutType_ = kAlwaysFalse;
    return;
  }
  vector<TString> cTypes_s({"!=", "==", "<=", ">=", "<", ">"});
  vector<comparisonType> cTypes({kNotEqual, kEqual, kLessEqual, kGreaterEqual, kLess, kGreater});
  TString var(cut), val(cut);
  int pos(-1);
  for(size_t ind(0); ind < cTypes.size(); ind++){
    pos = cut.Index(cTypes_s[ind]);
    if(pos == -1) continue;

    var.Remove(pos, var.Length());
    val.Remove(0, pos+cTypes_s[ind].Length());
    assignBranch(var, val);
    compType_ = cTypes[ind];
    break;    
  } // Loop over comparison types
  if(pos == -1) assignBranch(var, val); // Assigning boolean branches
}

void onecut::assignBranch(TString var, TString val){
  if(var == "ht"){
    cutType_ = kFloat;
    bf_ = &baby_base::ht;
  }else if(var=="met"){
    cutType_ = kFloat;
    bf_ = &baby_base::met;
  }else if(var=="ht_ra2"){
    cutType_ = kFloat;
    bf_ = &baby_base::ht_ra2;
  }else if(var=="ht_clean"){
    cutType_ = kFloat;
    bf_ = &baby_base::ht_clean;
  }else if(var=="mht"){
    cutType_ = kFloat;
    bf_ = &baby_base::mht;
  }else if(var=="mt"){
    cutType_ = kFloat;
    bf_ = &baby_base::mt;
  }else if(var=="elelv_pt"){
    cutType_ = kFloat;
    bf_ = &baby_base::elelv_pt;
  }else if(var=="elel_pt"){
    cutType_ = kFloat;
    bf_ = &baby_base::elel_pt;
  }else if(var=="elelv_m"){
    cutType_ = kFloat;
    bf_ = &baby_base::elelv_m;
  }else if(var=="elel_m"){
    cutType_ = kFloat;
    bf_ = &baby_base::elel_m;
  }else if(var=="mumuv_pt"){
    cutType_ = kFloat;
    bf_ = &baby_base::mumuv_pt;
  }else if(var=="mumu_pt"){
    cutType_ = kFloat;
    bf_ = &baby_base::mumu_pt;
  }else if(var=="mumuv_m"){
    cutType_ = kFloat;
    bf_ = &baby_base::mumuv_m;
  }else if(var=="mumu_m"){
    cutType_ = kFloat;
    bf_ = &baby_base::mumu_m;
  }else if(var=="mj"){
    cutType_ = kFloat;
    bf_ = &baby_base::mj;
  }else if(var=="nleps"){
    cutType_ = kInt;
    bi_ = &baby_base::nleps;
  }else if(var=="nvels"){
    cutType_ = kInt;
    bi_ = &baby_base::nvels;
  }else if(var=="nels"){
    cutType_ = kInt;
    bi_ = &baby_base::nels;
  }else if(var=="nvmus"){
    cutType_ = kInt;
    bi_ = &baby_base::nvmus;
  }else if(var=="nmus"){
    cutType_ = kInt;
    bi_ = &baby_base::nmus;
  }else if(var=="ntruleps"){
    cutType_ = kInt;
    bi_ = &baby_base::ntruleps;
  }else if(var=="njets"){
    cutType_ = kInt;
    bi_ = &baby_base::njets;
  }else if(var=="njets_ra2"){
    cutType_ = kInt;
    bi_ = &baby_base::njets_ra2;
  }else if(var=="njets_clean"){
    cutType_ = kInt;
    bi_ = &baby_base::njets_clean;
  }else if(var=="run"){
    cutType_ = kInt;
    bi_ = &baby_base::run;
  }else if(var=="nbm"){
    cutType_ = kInt;
    bi_ = &baby_base::nbm;
  }else if(var=="pass"){
    cutType_ = kBool;
    bb_ = &baby_base::pass;
  }else if(var=="stitch"){
    cutType_ = kBool;
    bb_ = &baby_base::stitch;
  }else if(var=="pass_ra2"){
    cutType_ = kBool;
    bb_ = &baby_base::pass_ra2;
  }else if(var=="pass_jets"){
    cutType_ = kBool;
    bb_ = &baby_base::pass_jets;
  } else if(var.Contains("[")){ // if var is a vector element
    TString index_s(var);
    var.Remove(var.Index("["), var.Length());
    index_s.Remove(0, index_s.Index("[")+1);
    index_s.Remove(index_s.Index("]"), index_s.Length());
    ivector_ = index_s.Atoi();
    if(var=="trig"){
      cutType_ = kvBool;
      bvb_ = &baby_base::trig;
    }else if(var=="sys_pass"){
      cutType_ = kvBool;
      bvb_ = &baby_base::sys_pass;
    }else if(var=="sys_ht"){ 
      cutType_ = kvFloat;
      bvf_ = &baby_base::sys_ht;
    }else if(var=="sys_met"){ 
      cutType_ = kvFloat;
      bvf_ = &baby_base::sys_met;
    }else if(var=="leps_pt"){ 
      cutType_ = kvFloat;
      bvf_ = &baby_base::leps_pt;
    }else if(var=="sys_njets"){ 
      cutType_ = kvInt;
      bvi_ = &baby_base::sys_njets;
    }else if(var=="sys_nbm"){ 
      cutType_ = kvInt;
      bvi_ = &baby_base::sys_nbm;
    }else if(var=="sys_mj"){ 
      cutType_ = kvFloat;
      bvf_ = &baby_base::sys_mj;
    }else if(var=="sys_mt"){ 
      cutType_ = kvFloat;
      bvf_ = &baby_base::sys_mt;
    }else {
      cout<<"Branch \""<<var<<" not defined. Add it to onecut::assignBranch in bcut.cpp"<<endl;
      exit(0);
    }
  }else {
    cout<<"Branch \""<<var<<" not defined. Add it to onecut::assignBranch in bcut.cpp"<<endl;
    exit(0);
  } 

  if(cutType_ == kFloat || cutType_ == kvFloat) cutf_ = val.Atof();
  if(cutType_ == kInt || cutType_ == kvInt)     cuti_ = val.Atoi();
}
