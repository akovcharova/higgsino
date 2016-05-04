// Cut class

#ifndef H_BCUT
#define H_BCUT

#include <iostream>
#include <vector>

#include "TString.h"
#include "baby_full.hpp"

typedef float& (baby_base::*baby_float)();
typedef bool& (baby_base::*baby_bool)();
typedef int& (baby_base::*baby_int)();
typedef std::vector<bool>& (baby_base::*baby_vbool)();
typedef std::vector<int>& (baby_base::*baby_vint)();
typedef std::vector<float>& (baby_base::*baby_vfloat)();

class onecut {
public:
  // Allows selection of different cut comparisons
  enum comparisonType{
    kNotEqual = -3,
    kLess = -2,
    kLessEqual = -1,
    kEqual = 0,
    kGreater = 1,
    kGreaterEqual = 2
  };

  enum cutType{kFloat, kvFloat, kInt, kvInt, kBool, kvBool, kAlwaysTrue, kAlwaysFalse};
  cutType cutType_;
  comparisonType compType_;
  TString cut_;
  int ivector_;

  // Booleans
  baby_bool bb_;
  baby_vbool bvb_;
  // Floats
  baby_float bf_;
  float cutf_;
  baby_vfloat bvf_;
  // Ints
  baby_int bi_;
  baby_vint bvi_;
  int cuti_;

  bool pass(baby_base *baby);
  void parseCut(TString cut);
  void assignBranch(TString var, TString val);
  onecut(TString cut);
  ~onecut();
};

class bcut {
public:
  std::vector<onecut> vcuts_;
  TString cuts_, weights_;
  enum cutType{kFloat, kvFloat, kConst};
  std::vector<cutType> cutTypes_;
  std::vector<baby_float> fWeights_;
  std::vector<baby_vfloat> fvWeights_;
  std::vector<int> indWeights_;
  std::vector<float> constWeights_;

  void parseWeights(TString weights);
  void parseWeight(TString weight);
  void parseCuts(TString cuts);
  bool pass(baby_base *baby);
  float weight(baby_base *baby);
  void operator+=(TString &cut);
  bcut operator+(bcut &ibcut);
  bcut(TString cuts="", TString weights="weight");
  ~bcut();

};

#endif
