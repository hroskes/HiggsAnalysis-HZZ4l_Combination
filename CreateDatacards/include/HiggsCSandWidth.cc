#ifndef HIGGSCSANDWIDTH_CC
#define HIGGSCSANDWIDTH_CC


#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <fstream>

#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"


#include "HiggsCSandWidth.h"

using namespace std;

HiggsCSandWidth::HiggsCSandWidth(std::string report, std::string FileLoc)
{

  max_MH_YR2 = 1000;
  min_MH_YR2 = 90;
  
  max_MH_YR3 = 1000; 
  min_MH_YR3 = 80;

  max_MH_YR2_Assoc = 300;
  max_MH_YR3_Assoc = 400;

  Report = report;
  fileLoc = FileLoc;

  if( report != "YR3" && report != "YR2")
    {
      cout << "Unknown report! Please choose YR3 or YR2." << endl;
      exit(1);
    }

  if( report == "YR3" )
    {
      max_MH = max_MH_YR3;
      min_MH = min_MH_YR3;
      max_MH_Assoc = max_MH_YR3_Assoc;
    }
  if( report == "YR2" )
    {
      max_MH = max_MH_YR2; 
      min_MH = min_MH_YR2;
      max_MH_Assoc = max_MH_YR2_Assoc;
    }


  ifstream file;

  int k = 0;
  fileLoc += "/"+report;
  // ---------------- Read BR into memory ------------------ //         
  fileName = fileLoc+"/HiggsBR_Official.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName << endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_BR[k] >> BR[k][1] >> BR[k][2] >> BR[k][3] >> BR[k][4] >> BR[k][5] >> BR[k][6] >> BR[k][7] >> BR[k][8] >> BR[k][9]
	   >> BR[k][10] >> BR[k][11] >> BR[k][12] >> BR[k][13] >> BR[k][14] >> BR[k][15] >> BR[k][16] >> BR[k][17] >> BR[k][18] >> BR[k][19] >> BR[k][20]
	   >> BR[k][21] >> BR[k][22] >> BR[k][23] >> BR[k][24] >> BR[k][25];
      k++;
      
    }
  N_BR = k;
  file.close();

  
  k = 0;
  fileName = fileLoc+"/HiggsBR_UpError_Official.txt";
  file.open(fileName.c_str());
  if( report == "YR3")
    {
      if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
      while (file.good())
	{
	  
	  file >> mass_BRerrPlus[k] >> BRerrPlus[k][1] >> BRerrPlus[k][2] >> BRerrPlus[k][3] >> BRerrPlus[k][4] 
	       >> BRerrPlus[k][5] >> BRerrPlus[k][6] >> BRerrPlus[k][7] >> BRerrPlus[k][8] >> BRerrPlus[k][9]
	       >> BRerrPlus[k][10] >> BRerrPlus[k][11] >> BRerrPlus[k][12] >> BRerrPlus[k][13] >> BRerrPlus[k][14] 
	       >> BRerrPlus[k][15] >> BRerrPlus[k][16] >> BRerrPlus[k][17] >> BRerrPlus[k][18] >> BRerrPlus[k][19] 
	       >> BRerrPlus[k][20] >> BRerrPlus[k][21] >> BRerrPlus[k][22] >> BRerrPlus[k][23] >> BRerrPlus[k][24] >> BRerrPlus[k][25];
	  k++;
	}
      if( k != N_BR ){cout << "N_BR does not match between HiggsBR_UpError_Official.txt and HiggsBR_Official.txt!" << endl; exit(1);}
      file.close();
    }
  else
    {
      for(int l = 0; l < N_BR; l++)
	{
	  mass_BRerrPlus[l] = mass_BR[k];
	  for(int m = 0; m < 25; m++)
	    {
	      BRerrPlus[l][m] = 0;
	      BRerrMinus[l][m] = 0;
	    }	
	}
    }

  k = 0;
  fileName = fileLoc+"/HiggsBR_DnError_Official.txt";
  file.open(fileName.c_str());
  if( report == "YR3" )
    {
      if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
      while (file.good())
	{
	  
	  file >> mass_BRerrMinus[k] >> BRerrMinus[k][1] >> BRerrMinus[k][2] >> BRerrMinus[k][3] >> BRerrMinus[k][4] 
	       >> BRerrMinus[k][5] >> BRerrMinus[k][6] >> BRerrMinus[k][7] >> BRerrMinus[k][8] >> BRerrMinus[k][9]
	       >> BRerrMinus[k][10] >> BRerrMinus[k][11] >> BRerrMinus[k][12] >> BRerrMinus[k][13] >> BRerrMinus[k][14] 
	       >> BRerrMinus[k][15] >> BRerrMinus[k][16] >> BRerrMinus[k][17] >> BRerrMinus[k][18] >> BRerrMinus[k][19] 
	       >> BRerrMinus[k][20] >> BRerrMinus[k][21] >> BRerrMinus[k][22] >> BRerrMinus[k][23] >> BRerrMinus[k][24] >> BRerrMinus[k][25];
	  k++;
	}
      if( k != N_BR ){cout << "N_BR does not match between HiggsBR_DnError_Official.txt and HiggsBR_Official.txt!" << endl;exit(1);}
      file.close();
    }

  k = 0;
  fileName = fileLoc+"/HiggsTotalWidth.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  if( report == "YR3" )
    {
      while (file.good())
	{
	  
	  file >> scratchMass >> BR[k][0] >> BRerrPlus[k][0] >> BRerrMinus[k][0];
	  k++;
	}
    }
  else{
    while (file.good())
      {
	
	file >> scratchMass >> BR[k][0];
	k++;
      }
  }
  if( k != N_BR ){cout << "N_BR does not match between HiggsTotalWidth.txt and HiggsBR_Official.txt!" << endl;exit(1);}
  file.close();
      
  // ---------------- Read 7 TeV CS into memory ------------------ //         
  k = 0;
  fileName = fileLoc+"/7TeV-ggH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_7tev[k][ID_ggToH] >> CS_7tev[k][ID_ggToH] >> CSerrPlus_7tev[k][ID_ggToH] >> CSerrMinus_7tev[k][ID_ggToH] 
	   >> CSscaleErrPlus_7tev[k][ID_ggToH] >> CSscaleErrMinus_7tev[k][ID_ggToH] >> CSpdfErrPlus_7tev[k][ID_ggToH] >> CSpdfErrMinus_7tev[k][ID_ggToH];
      k++;
    }
  N_CS_7tev[ID_ggToH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/7TeV-vbfH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_7tev[k][ID_VBF] >> CS_7tev[k][ID_VBF] >> CSerrPlus_7tev[k][ID_VBF] >> CSerrMinus_7tev[k][ID_VBF] >> CSscaleErrPlus_7tev[k][ID_VBF]
	   >> CSscaleErrMinus_7tev[k][ID_VBF] >> CSpdfErrPlus_7tev[k][ID_VBF] >> CSpdfErrMinus_7tev[k][ID_VBF];
      k++;
    }
  N_CS_7tev[ID_VBF] = k;
  file.close();
      
  k = 0;
  fileName = fileLoc+"/7TeV-ttH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_7tev[k][ID_ttH] >> CS_7tev[k][ID_ttH] >> CSerrPlus_7tev[k][ID_ttH] >> CSerrMinus_7tev[k][ID_ttH] >> CSscaleErrPlus_7tev[k][ID_ttH]
	   >> CSscaleErrMinus_7tev[k][ID_ttH] >> CSpdfErrPlus_7tev[k][ID_ttH] >> CSpdfErrMinus_7tev[k][ID_ttH];
      k++;
    }
  N_CS_7tev[ID_ttH] = k; 
  file.close();
      
  k = 0;
  fileName = fileLoc+"/7TeV-ZH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_7tev[k][ID_ZH] >> CS_7tev[k][ID_ZH] >> CSerrPlus_7tev[k][ID_ZH] >> CSerrMinus_7tev[k][ID_ZH] >> CSscaleErrPlus_7tev[k][ID_ZH]
	   >> CSscaleErrMinus_7tev[k][ID_ZH] >> CSpdfErrPlus_7tev[k][ID_ZH] >> CSpdfErrMinus_7tev[k][ID_ZH];
      k++;
    }
  N_CS_7tev[ID_ZH] = k;
  file.close();
     
  k = 0;
  fileName = fileLoc+"/7TeV-WH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_7tev[k][ID_WH] >> CS_7tev[k][ID_WH] >> CSerrPlus_7tev[k][ID_WH] >> CSerrMinus_7tev[k][ID_WH] >> CSscaleErrPlus_7tev[k][ID_WH]
	   >> CSscaleErrMinus_7tev[k][ID_WH] >> CSpdfErrPlus_7tev[k][ID_WH] >> CSpdfErrMinus_7tev[k][ID_WH];
      k++;
    }
  N_CS_7tev[ID_WH] = k;
  file.close();
      
  // ---------------- Read 8 TeV CS into memory ------------------ //         
  k = 0;
  fileName = fileLoc+"/8TeV-ggH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
	
      file >> mass_XS_8tev[k][ID_ggToH] >> CS_8tev[k][ID_ggToH] >> CSerrPlus_8tev[k][ID_ggToH] >> CSerrMinus_8tev[k][ID_ggToH] 
	   >> CSscaleErrPlus_8tev[k][ID_ggToH] >> CSscaleErrMinus_8tev[k][ID_ggToH] >> CSpdfErrPlus_8tev[k][ID_ggToH] >> CSpdfErrMinus_8tev[k][ID_ggToH];
      k++;
  }
  N_CS_8tev[ID_ggToH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/8TeV-vbfH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      file >> mass_XS_8tev[k][ID_VBF] >> CS_8tev[k][ID_VBF] >> CSerrPlus_8tev[k][ID_VBF] >> CSerrMinus_8tev[k][ID_VBF] >> CSscaleErrPlus_8tev[k][ID_VBF]
	   >> CSscaleErrMinus_8tev[k][ID_VBF] >> CSpdfErrPlus_8tev[k][ID_VBF] >> CSpdfErrMinus_8tev[k][ID_VBF];
      k++;
    }
  N_CS_8tev[ID_VBF] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/8TeV-ttH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_8tev[k][ID_ttH] >> CS_8tev[k][ID_ttH] >> CSerrPlus_8tev[k][ID_ttH] >> CSerrMinus_8tev[k][ID_ttH] >> CSscaleErrPlus_8tev[k][ID_ttH]
	   >> CSscaleErrMinus_8tev[k][ID_ttH] >> CSpdfErrPlus_8tev[k][ID_ttH] >> CSpdfErrMinus_8tev[k][ID_ttH];
      k++;
    }
  N_CS_8tev[ID_ttH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/8TeV-ZH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_8tev[k][ID_ZH] >> CS_8tev[k][ID_ZH] >> CSerrPlus_8tev[k][ID_ZH] >> CSerrMinus_8tev[k][ID_ZH] >> CSscaleErrPlus_8tev[k][ID_ZH]
	   >> CSscaleErrMinus_8tev[k][ID_ZH] >> CSpdfErrPlus_8tev[k][ID_ZH] >> CSpdfErrMinus_8tev[k][ID_ZH];
      k++;
    }
  N_CS_8tev[ID_ZH] = k;
  file.close();
      
  k = 0;
  fileName = fileLoc+"/8TeV-WH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_8tev[k][ID_WH] >> CS_8tev[k][ID_WH] >> CSerrPlus_8tev[k][ID_WH] >> CSerrMinus_8tev[k][ID_WH] >> CSscaleErrPlus_8tev[k][ID_WH]
	   >> CSscaleErrMinus_8tev[k][ID_WH] >> CSpdfErrPlus_8tev[k][ID_WH] >> CSpdfErrMinus_8tev[k][ID_WH];
      k++;
  }
  N_CS_8tev[ID_WH] = k;
  file.close();
      
      
  // ---------------- Read 14 TeV CS into memory ------------------ //         
  k = 0;
  fileName = fileLoc+"/14TeV-ggH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_14tev[k][ID_ggToH] >> CS_14tev[k][ID_ggToH] >> CSerrPlus_14tev[k][ID_ggToH] >> CSerrMinus_14tev[k][ID_ggToH] 
	   >> CSscaleErrPlus_14tev[k][ID_ggToH] >> CSscaleErrMinus_14tev[k][ID_ggToH] >> CSpdfErrPlus_14tev[k][ID_ggToH] >> CSpdfErrMinus_14tev[k][ID_ggToH];
      k++;
  }
  N_CS_14tev[ID_ggToH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/14TeV-vbfH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_14tev[k][ID_VBF] >> CS_14tev[k][ID_VBF] >> CSerrPlus_14tev[k][ID_VBF] >> CSerrMinus_14tev[k][ID_VBF] >> CSscaleErrPlus_14tev[k][ID_VBF]
	   >> CSscaleErrMinus_14tev[k][ID_VBF] >> CSpdfErrPlus_14tev[k][ID_VBF] >> CSpdfErrMinus_14tev[k][ID_VBF];
      k++;
  }
  N_CS_14tev[ID_VBF] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/14TeV-ttH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_14tev[k][ID_ttH] >> CS_14tev[k][ID_ttH] >> CSerrPlus_14tev[k][ID_ttH] >> CSerrMinus_14tev[k][ID_ttH] >> CSscaleErrPlus_14tev[k][ID_ttH]
	   >> CSscaleErrMinus_14tev[k][ID_ttH] >> CSpdfErrPlus_14tev[k][ID_ttH] >> CSpdfErrMinus_14tev[k][ID_ttH];
      k++;
  }
  N_CS_14tev[ID_ttH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/14TeV-ZH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_14tev[k][ID_ZH] >> CS_14tev[k][ID_ZH] >> CSerrPlus_14tev[k][ID_ZH] >> CSerrMinus_14tev[k][ID_ZH] >> CSscaleErrPlus_14tev[k][ID_ZH]
	   >> CSscaleErrMinus_14tev[k][ID_ZH] >> CSpdfErrPlus_14tev[k][ID_ZH] >> CSpdfErrMinus_14tev[k][ID_ZH];
      k++;
    }
  N_CS_14tev[ID_ZH] = k;
  file.close();

  k = 0;
  fileName = fileLoc+"/14TeV-WH.txt";
  file.open(fileName.c_str());
  if (!file.is_open()){ cout << "Could not find file "+fileName <<endl; exit(1);}
  while (file.good())
    {
      
      file >> mass_XS_14tev[k][ID_WH] >> CS_14tev[k][ID_WH] >> CSerrPlus_14tev[k][ID_WH] >> CSerrMinus_14tev[k][ID_WH] >> CSscaleErrPlus_14tev[k][ID_WH]
	   >> CSscaleErrMinus_14tev[k][ID_WH] >> CSpdfErrPlus_14tev[k][ID_WH] >> CSpdfErrMinus_14tev[k][ID_WH];
      k++;
  }
  N_CS_14tev[ID_WH] = k;
  file.close();



}


HiggsCSandWidth::~HiggsCSandWidth()
{
  //destructor

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV 
double HiggsCSandWidth::HiggsCS(int ID, double mH, double sqrts){

  /**********IDs*************/ 
  /*     ggToH = 1          */
  /*       VBF = 2          */ 
  /*        WH = 3          */ 
  /*        ZH = 4          */
  /*       ttH = 5          */
  /*     Total = 0          */
  /**************************/
 
  double val = 0;

  // If ID is unavailable return -1                                                                                                
  if(ID > ID_ttH || ID < ID_Total) return -1;
  // If Ecm is not 7 or 8 TeV return -1
  if(sqrts != 7 && sqrts != 8 && sqrts != 14) return -1;
  //Don't interpolate btw 0 and numbers for mH 400
  if(ID > ID_VBF && mH > max_MH_Assoc) return 0;
 
 // If mH is out of range return -1                                           
  // else find what array number to read         
  if( mH < min_MH || mH > max_MH){ return -1;}
  else{
    
    if(sqrts == 7)
      {
	if(ID == 0)
	  {
	    for(int i = 1; i <= 5; i++)
	      {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		val += getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CS_7tev);
	      }
	  }
	else val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CS_7tev);
      }
    else if (sqrts == 8)
      {
	if(ID == 0)
	  {
	    for(int i = 1; i <= 5; i++)
	      {
                if (mH > max_MH_Assoc && i > ID_VBF) continue;
		val += getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CS_8tev);
	      }
	  }
	else val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CS_8tev);
      }
    else if(sqrts == 14)
      {
        if(ID == 0)
          {
            for(int i =1; i <= 5; i++)
	      { 
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		val += getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CS_14tev);
	      } 
          }
	else val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CS_14tev);
      }
    else{cout << "HiggsCSandWidth::HiggsCS --- unknown sqrts! Choose 7,8, or 14." << endl; return -1;}
  }  

  return val;


}
  
  

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV 
double HiggsCSandWidth::HiggsCSErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;


  // If ID is unavailable return -1                                                                                    
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                                
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400                        
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}
 
  // If mH is out of range return -1                                                                        
  // else find what array number to read                                          
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    if(sqrts == 7)
      {
	if(ID == ID_Total)
	  {
	    double tmpVal = 0;
	    for(int i = 1; i <= 5; i++)
	      {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSerrPlus_7tev);
		tmpVal += interp*interp;
	      }
	    val = sqrt(tmpVal);
	  }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSerrPlus_7tev);}
      }
    else if (sqrts == 8)
      {
	if(ID == ID_Total)
          {
            double tmpVal = 0;
	    for(int i = 1; i <= 5; i++)
	      { 
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSerrPlus_8tev);
		tmpVal += interp*interp;
	      }
            val= sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSerrPlus_8tev);}
      }
    else if(sqrts == 14)
      {
	if(ID == ID_Total)
          {
            double tmpVal = 0;
	    for(int i = 1; i <= 5; i++)
	      { 
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSerrPlus_14tev);
		tmpVal += interp*interp;
	      }
            val= sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSerrPlus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}
  }

  val *= .01;
  
  return val;
  
}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                                                                                       
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                           
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400                                                        
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}

  // If mH is out of range return -1                                                                           
  // else find what array number to read                                                                 
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    if(sqrts == 7)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
	    for(int i = 1; i <= 5; i++)
	      {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSerrMinus_7tev);
		tmpVal += interp*interp;
	      }
            val = -1*sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSerrMinus_7tev);}
      }
    else if (sqrts == 8)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSerrMinus_8tev);
                tmpVal += interp*interp;
              }
            val = -1*sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSerrMinus_8tev);}
      }
    else if(sqrts == 14)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSerrMinus_14tev);
                tmpVal += interp*interp;
              }
            val = -1*sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSerrMinus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }

  val *= .01;

  return val;

}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSscaleErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/


  double val = 0;

  // If ID is unavailable return -1                                                         
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400                                           
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}

  // If mH is out of range return -1                                                         
  // else find what array number to read                                                      
  if( mH < min_MH || mH > max_MH){return -1;}
  else{
    
    
    if(sqrts == 7)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSscaleErrPlus_7tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSscaleErrPlus_7tev);}
      }
    else if (sqrts == 8)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSscaleErrPlus_8tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSscaleErrPlus_8tev);}
      }
    else if(sqrts == 14)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSscaleErrPlus_14tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSscaleErrPlus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSscaleErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  
  val *= .01; //Account for percentage  
  
  return val;
  
}

//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSscaleErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                     
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                               
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400                                   
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}

  // If mH is out of range return -1                        
  // else find what array number to read                              
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    if(sqrts == 7)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSscaleErrMinus_7tev);
                tmpVal += interp*interp;
              }
            val = -1*sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSscaleErrMinus_7tev);}
      }
    else if (sqrts == 8)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSscaleErrMinus_8tev);
                tmpVal += interp*interp;
              }
            val = -1*sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSscaleErrMinus_8tev);}
      }
    else if(sqrts == 14)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSscaleErrMinus_14tev);
                tmpVal += interp*interp;
              }
            val = -1*sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSscaleErrMinus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSscaleErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}
  }
  
  val *= .01; //Account for percentage  
  
  return val;

}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSpdfErrPlus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                                                                           
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                                                         
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400                                                  
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}


  // If mH is out of range return -1                                                                                  
  // else find what array number to read                                                              
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    if(sqrts == 7)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSpdfErrPlus_7tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
	else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSpdfErrPlus_7tev);}
      }
    else if (sqrts == 8)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSpdfErrPlus_8tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSpdfErrPlus_8tev);}
      }
    else if(sqrts == 14)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSpdfErrPlus_14tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSpdfErrPlus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSpdfErrPlus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  val *= .01; //Account for percentage  
  
  return val;


}


//Higgs CS takes process ID, higgs mass mH, and COM energy sqrts in TeV
double HiggsCSandWidth::HiggsCSpdfErrMinus(int ID, double mH, double sqrts){

  /**********IDs*************/
  /*     ggToH = 1          */
  /*       VBF = 2          */
  /*        WH = 3          */
  /*        ZH = 4          */
  /*       ttH = 5          */
  /**************************/

  double val = 0;

  // If ID is unavailable return -1                           
  if(ID > ID_ttH || ID < ID_Total){return -1;}
  // If Ecm is not 7 or 8 TeV return -1                                                                 
  if(sqrts != 7 && sqrts != 8 && sqrts != 14){return -1;}
  //Don't interpolate btw 0 and numbers for mH 400             
  if(ID > ID_VBF && mH > max_MH_Assoc){return 0;}


  // If mH is out of range return -1                                                              
  // else find what array number to read                            
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    if(sqrts == 7)
      {
	if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_7tev[i],mass_XS_7tev,CSpdfErrMinus_7tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_7tev[ID],mass_XS_7tev,CSpdfErrMinus_7tev);}
      }
    else if (sqrts == 8)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_8tev[i],mass_XS_8tev,CSpdfErrMinus_8tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_8tev[ID],mass_XS_8tev,CSpdfErrMinus_8tev);}
      }
    else if(sqrts == 14)
      {
        if(ID == ID_Total)
          {
            double tmpVal = 0;
            for(int i = 1; i <= 5; i++)
              {
		if (mH > max_MH_Assoc && i > ID_VBF) continue;
		double interp = getInterpXS(sqrts,i,mH,N_CS_14tev[i],mass_XS_14tev,CSpdfErrMinus_14tev);
                tmpVal += interp*interp;
              }
            val = sqrt(tmpVal);
          }
        else{val = getInterpXS(sqrts,ID,mH,N_CS_14tev[ID],mass_XS_14tev,CSpdfErrMinus_14tev);}
      }
    else{cout << "HiggsCSandWidth::HiggsCSpdfErrMinus --- unknown sqrts! Choose 7 or 8." << endl; return -1;}

  }
  
  val *= .01; //Account for percentage  
  
  return val;

}



// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsWidth(int ID, double mH){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/

  double Width = 0;

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 0){return -1;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    Width = getInterpBRWidth(true, ID, mH, N_BR, mass_BR, BR);

  }

  return Width;

} 


// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsBR(int ID, double mH){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/


  double BranchRatio = 0;

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 1){return -1;}


  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    BranchRatio = getInterpBRWidth(false, ID, mH, N_BR, mass_BR, BR);

  }

  return BranchRatio;

} 



// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsBRerrPlus(int ID, double mH){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/

  double BRunc = 0;

  if( Report == "YR2" ){ return 0;}

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 1){return -1;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    BRunc = getInterpBRWidth(false, ID, mH, N_BR, mass_BRerrPlus, BRerrPlus);

  }

  BRunc *=.01; //Account for percentage

  return BRunc;

} 

// HiggsWidth takes process ID and higgs mass mH
double HiggsCSandWidth::HiggsBRerrMinus(int ID, double mH){


  /***********************IDs************************/
  /*                       Total = 0                */
  /*                       H->bb = 1                */
  /*                   H->tautau = 2                */
  /*                     H->mumu = 3                */
  /*                       H->ss = 4                */
  /*                       H->cc = 5                */
  /*                       H->tt = 6                */
  /*                       H->gg = 7                */
  /*                   H->gamgam = 8                */
  /*                     H->gamZ = 9                */
  /*                       H->WW = 10               */
  /*                       H->ZZ = 11               */
  /*                       H->4e = 12               */
  /*                    H->2e2mu = 13               */
  /*              H->4lep (e,mu) = 14               */
  /*          H->4lep (e,mu,tau) = 15               */
  /*                H->e+nu e-nu = 16               */
  /*               H->e+nu mu-nu = 17               */
  /*    H->2l2nu(l=e,mu)(nu=any) = 18               */
  /* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
  /*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
  /* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
  /* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
  /*  H->2nu2q (nu=any)(q=udcsb) = 23               */
  /*            H->4q (q=udcsb) = 24                */
  /*      H->4f (f=any fermion) = 25                */
  /**************************************************/

  double BRunc = 0;

  if( Report == "YR2" ){return 0;}

  // If ID is unavailable return -1                                           
  if(ID > 25 || ID < 1){return -1;}

  // If mH is out of range return -1                                            
  // else find what array number to read                                        
  if( mH < min_MH || mH > max_MH){return -1;}
  else{

    BRunc = getInterpBRWidth(false, ID, mH, N_BR, mass_BRerrMinus, BRerrMinus);

  }

  BRunc *=.01; //Account for percentage

  return BRunc;

} 




double HiggsCSandWidth::getInterpXS(int sqrts, int ID, double mH, int maxI, double mhArray[][6], double varArray[][6])
{

  using namespace std;

  int i = 0;
  double reqCS = 0;
  const int index = 4;
  double xmh[index], sig[index];
  
  i = searchArray(ID, maxI, mhArray, mH);

  //Do the interpolation
  if(i < 1){i = 1;}
  if(i+2 >= maxI){i = maxI - 3;}
  xmh[0]=mhArray[i-1][ID];  xmh[1]=mhArray[i][ID];  xmh[2]=mhArray[i+1][ID];  xmh[3]=mhArray[i+2][ID];
  sig[0]=varArray[i-1][ID]; sig[1]=varArray[i][ID]; sig[2]=varArray[i+1][ID]; sig[3]=varArray[i+2][ID];
  
  TGraph *graph = new TGraph(index, xmh, sig);
  TSpline3 *gs = new TSpline3("gs",graph);
  gs->Draw();
  reqCS = gs->Eval(mH);
  delete gs;
  delete graph;
  
  return reqCS;  
  
}



double HiggsCSandWidth::getInterpBRWidth(bool width, int ID, double mH, int maxI, double mhArray[], double varArray[][26])
{

  using namespace std;

  int i = 0;
  double reqBR = 0;

  i = searchArray(ID, maxI, mhArray, mH);

  if(i < 1){i = 1;}
  if(i+2 >= maxI){i = maxI - 3;}
  const int indexW = 4;
  double xmhW[indexW], sigW[indexW];

  xmhW[0]=mhArray[i-1];xmhW[1]=mhArray[i];xmhW[2]=mhArray[i+1];xmhW[3]=mhArray[i+2];
  sigW[0]=varArray[i-1][ID]; sigW[1]=varArray[i][ID]; sigW[2]=varArray[i+1][ID]; sigW[3]=varArray[i+2][ID];
  
  TGraph *graphW = new TGraph(indexW, xmhW, sigW);
  TSpline3 *gsW = new TSpline3("gsW",graphW);
  gsW->Draw();
  reqBR = gsW->Eval(mH);
  delete gsW;
  delete graphW;
  
  if(width && ID > 0)
    {
      const int indexPW = 4;
      double xmhPW[indexPW], sigPW[indexPW];
      xmhPW[0]=mhArray[i-1];xmhPW[1]=mhArray[i];xmhPW[2]=mhArray[i+1];xmhPW[3]=mhArray[i+2];
      sigPW[0]=varArray[i-1][0]; sigPW[1]=varArray[i][0]; sigPW[2]=varArray[i+1][0]; sigPW[3]=varArray[i+2][0];
      
      TGraph *graphPW = new TGraph(indexPW, xmhPW, sigPW);
      TSpline3 *gsPW = new TSpline3("gsPW",graphPW);
      gsPW->Draw();
      reqBR *= gsPW->Eval(mH);
      delete gsPW;
      delete graphPW;
    }

  
  return reqBR;  

}



int HiggsCSandWidth::searchArray(int ID, int maxI, double mhArray[][6], double mh)
{
  
  int min_dist = 10000; 
  int min_index = 0; 

  for( int i = 0; i < maxI; i++ )
    { 
      if( mhArray[i][ID] == mh ) { return i;} 
      else {
	if ( abs(mhArray[i][ID] - mh) < min_dist ) 
	  {
	    min_dist = abs(mhArray[i][ID] - mh);
	    min_index = i;
	    int dist = mhArray[i][ID] - mh;
	    if(dist > 0){ return min_index;}
	  }
      }
    }
  return min_index;

}

int HiggsCSandWidth::searchArray(int ID, int maxI, double mhArray[], double mh)
{
  
  int min_dist = 10000; 
  int min_index = 0; 

  for( int i = 0; i < maxI; i++ )
    { 
      if( mhArray[i] == mh ) { return i;} 
      else {
	if ( abs(mhArray[i] - mh) < min_dist ) 
	  {
	    min_dist = abs(mhArray[i] - mh);
	    min_index = i;
	    int dist = mhArray[i] - mh;
	    if(dist > 0){ return min_index;}
	  }
      }
    }
  return min_index;

}


#endif
