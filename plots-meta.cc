// ./X

// nutotal vs offaxis ( mHNL = 0.1, ... , 1.9 )

#include <iostream> 
#include <fstream> 
#include <vector>

#include <stdlib.h>
#include <gsl/gsl_histogram.h>

#include <TH1.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

int main(){
	
	// ******************************************************
	// Verificar ENERGÍA 	
	TFile f1("/home/sane/mainhnlchain/dirac/120/meson2hnl/data/120-431.root","open"); // INPUT

	
	TTree* nui;

	f1.GetObject("nu",nui);
	
	double rebeam, rmhnl, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
  rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex;
  int auxpdg, rdet_id, rmotherid;

	nui->SetBranchAddress("ebeam",&rebeam);
  nui->SetBranchAddress("mhnl",&rmhnl);
  nui->SetBranchAddress("offaxis",&roffaxis);
  nui->SetBranchAddress("det_id",&rdet_id);
  nui->SetBranchAddress("id",&rid);
  nui->SetBranchAddress("tProd",&rtProd);
  nui->SetBranchAddress("xProd",&rxProd);
  nui->SetBranchAddress("yProd",&ryProd);
  nui->SetBranchAddress("zProd",&rzProd);
  nui->SetBranchAddress("tDec",&rtDec);
  nui->SetBranchAddress("xDec",&rxDec);
  nui->SetBranchAddress("yDec",&ryDec);
  nui->SetBranchAddress("zDec",&rzDec);
  nui->SetBranchAddress("e",&re);
  nui->SetBranchAddress("px",&rpx);
  nui->SetBranchAddress("py",&rpy);
  nui->SetBranchAddress("pz",&rpz);
  nui->SetBranchAddress("pT",&rpT);
  nui->SetBranchAddress("theta",&rtheta);
  nui->SetBranchAddress("phi",&rphi);
  nui->SetBranchAddress("y",&ry);
  nui->SetBranchAddress("eta",&reta);
  //nui->SetBranchAddress("motherid",&rmotherid); 
	
	int jEvents = nui->GetEntries();
	
	cout << "jEvents = " << jEvents << endl;
	
	// 0-mHNL y 41 offaxis [0, ..., 40]
  int D_mhnl = 15; // debería ser 15
  double v_meta1[D_mhnl][41][3];
  double v_meta2[D_mhnl][41][3];


  for( int i=0; i < D_mhnl; ++i){
    for (int j=0; j < 41; ++j){
      for(int k=0; k < 3; ++k){
        v_meta1[i][j][k]=0;
        v_meta2[i][j][k]=0;
      }
    }
  }

  int k = 0;
  int l = 0;

  int nEvents;
  nEvents = jEvents;

  cout << "nEvents = " << nEvents << endl;
  cout << "Realizando análisis..." << endl;
  // *********************  MAIN LOOP  **************************

  for (int i=0; i<nEvents; ++i){

		nui->GetEntry(i);
    
    if( (rdet_id==0    /*|| rdet_id==1*/   || rdet_id==2) &&
        rmhnl >= 0.5 ){
      
      k = int((rmhnl-0.5)/0.1+0.0001);
      //printf("%.2f %d \n", rmhnl, k);

      l = int(roffaxis+0.000001);

      switch(int(rid+0.0001)){
        case 12:
          v_meta1[k][l][0]++; 
          break;
        case 14:
          v_meta1[k][l][1]++; 
          break;
        case 16:
          v_meta1[k][l][2]++; 
          break;
        default:
          break;
      }

      switch(int(rid-0.0001)){
        case -12:
          v_meta2[k][l][0]++; 
          break;
        case -14:
          v_meta2[k][l][1]++; 
          break;
        case -16:
          v_meta2[k][l][2]++; 
          break;
        default:
          break;
      }     
            
    }
		
    if( (i+1) % 10000 == 0){	printf("%.2f% \r",1.*i/nEvents*100);	}
	}

  ofstream data1("meta-1.dat"); 
  ofstream data2("meta-2.dat"); 

  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        data1 << v_meta1[k][j][m] << " ";
      }
      data1 << endl;
    }
  }

  for (int k = 0; k < D_mhnl; ++k){
    for(int m = 0; m < 3; ++m){
      for (int j = 0; j < 41; ++j){
        data2 << v_meta2[k][j][m] << " ";
      }
      data2 << endl;
    }
  }
  
  cout << endl << "SUCCESS!" << endl;

  return 0;
	
}

