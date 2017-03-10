#ifndef HIGGSCSANDWIDTH_H
#define HIGGSCSANDWIDTH_H

#define PI 3.14159

#define  ID_ggToH  1
#define  ID_VBF    2
#define  ID_WH     3
#define  ID_ZH     4
#define  ID_ttH    5
#define  ID_Total  0 

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>

#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"


/**********************************************************/
/*            Class for Higgs Width and CS                */
/*                                                        */
/*  All numbers for CS and width are taken from official  */
/*  numbers on Higgs CS Twiki (Fall 2013)                 */
/*                                                        */
/*  Cross Sections are given in pb                        */
/*  Widths are given in GeV                               */
/*                                                        */
/*  These numbers are taken into memory and a spline      */
/*  interpolation is done       .                         */
/*                                                        */
/*  For any invalid process or mH out of range, -1 or 0   */
/*  willbe returned.                                      */
/*                                                        */
/*    Written by:                                         */
/*         Matt Snowball                                  */
/*         University of Florida                          */
/*         snowball@phys.ufl.edu                          */
/*                                                        */
/*       Last Update: April 14, 2014                      */
/*                                                        */
/**********************************************************/



class HiggsCSandWidth
{

 public:

  HiggsCSandWidth(std::string report = "YR3", std::string FileLoc = "include/txtFiles");
  ~HiggsCSandWidth();

  double HiggsCS(int ID, double mH, double sqrts);
  double HiggsCSErrPlus(int ID, double mH, double sqrts);
  double HiggsCSErrMinus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrPlus(int ID, double mH, double sqrts);
  double HiggsCSscaleErrMinus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrPlus(int ID, double mH, double sqrts);
  double HiggsCSpdfErrMinus(int ID, double mH, double sqrts);
  double HiggsWidth(int ID,double mH);
  double HiggsBR(int ID,double mH);
  double HiggsBRerrPlus(int ID,double mH);
  double HiggsBRerrMinus(int ID,double mH);
  double getInterpXS(int sqrts, int ID, double mH, int maxI, double mhArray[][6], double varArray[][6]);
  double getInterpBRWidth(bool width, int ID, double mH, int maxI, double mhArray[], double varArray[][26]);

 private:

  std::string fileName, fileLoc, Report;
  
  double max_MH_YR2, min_MH_YR2;
  double max_MH_YR3, min_MH_YR3;

  double max_MH_YR2_Assoc, max_MH_YR3_Assoc;

  double max_MH, min_MH;
  double max_MH_Assoc;

  double scratchMass;
  double mass_BR[500];
  double BR[500][26];

  double mass_BRerrPlus[500];
  double BRerrPlus[500][26];

  double mass_BRerrMinus[500];
  double BRerrMinus[500][26];

  double mass_XS_7tev[500][6];
  double CS_7tev[500][6];
  double CSerrPlus_7tev[500][6];
  double CSerrMinus_7tev[500][6];
  double CSscaleErrPlus_7tev[500][6];
  double CSscaleErrMinus_7tev[500][6];
  double CSpdfErrPlus_7tev[500][6];
  double CSpdfErrMinus_7tev[500][6];

  double mass_XS_8tev[500][6];
  double CS_8tev[500][6];
  double CSerrPlus_8tev[500][6];
  double CSerrMinus_8tev[500][6];
  double CSscaleErrPlus_8tev[500][6];
  double CSscaleErrMinus_8tev[500][6];
  double CSpdfErrPlus_8tev[500][6];
  double CSpdfErrMinus_8tev[500][6];

  double mass_XS_14tev[500][6];
  double CS_14tev[500][6];
  double CSerrPlus_14tev[500][6];
  double CSerrMinus_14tev[500][6];
  double CSscaleErrPlus_14tev[500][6];
  double CSscaleErrMinus_14tev[500][6];
  double CSpdfErrPlus_14tev[500][6];
  double CSpdfErrMinus_14tev[500][6];

  int N_BR;
  int N_CS_7tev[6];  
  int N_CS_8tev[6];  
  int N_CS_14tev[6];  


  int searchArray(int ID, int maxI, double mhArray[][6], double mh);
  int searchArray(int ID, int maxI, double mhArray[], double mh);

};

#endif
