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
  rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rpmother;
  int mother1, auxpdg, rdet_id;

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
  nui->SetBranchAddress("pmother",&rpmother);
	
	double treesize = nui->GetEntries();
	
	cout << treesize << endl;
	
	size_t nbins=20;
  
  gsl_histogram * h1_e00 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e01 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e02 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e03 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e04 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e10 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e11 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e12 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e13 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e14 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e20 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e21 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e22 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e23 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e24 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e30 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e31 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e32 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e33 = gsl_histogram_alloc (nbins);
  gsl_histogram * h1_e34 = gsl_histogram_alloc (nbins);

  gsl_histogram * h2_e00 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e01 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e02 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e03 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e04 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e10 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e11 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e12 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e13 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e14 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e20 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e21 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e22 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e23 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e24 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e30 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e31 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e32 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e33 = gsl_histogram_alloc (nbins);
  gsl_histogram * h2_e34 = gsl_histogram_alloc (nbins);  

  vector < gsl_histogram * > row;
  vector < vector < gsl_histogram * > > v_h1;
  vector < vector < gsl_histogram * > > v_h2;

  // nu histogram vector

  row.push_back(h1_e00);
  row.push_back(h1_e01);
  row.push_back(h1_e02);
  row.push_back(h1_e03);
  row.push_back(h1_e04);
  v_h1.push_back(row);
  row.clear();

  row.push_back(h1_e10);
  row.push_back(h1_e11);
  row.push_back(h1_e12);
  row.push_back(h1_e13);
  row.push_back(h1_e14);
  v_h1.push_back(row);
  row.clear();

  row.push_back(h1_e20);
  row.push_back(h1_e21);
  row.push_back(h1_e22);
  row.push_back(h1_e23);
  row.push_back(h1_e24);
  v_h1.push_back(row);
  row.clear();

  row.push_back(h1_e30);
  row.push_back(h1_e31);
  row.push_back(h1_e32);
  row.push_back(h1_e33);
  row.push_back(h1_e34);
  v_h1.push_back(row);
  row.clear();

  // nubar histogram vector

  row.push_back(h2_e00);
  row.push_back(h2_e01);
  row.push_back(h2_e02);
  row.push_back(h2_e03);
  row.push_back(h2_e04);
  v_h2.push_back(row);
  row.clear();

  row.push_back(h2_e10);
  row.push_back(h2_e11);
  row.push_back(h2_e12);
  row.push_back(h2_e13);
  row.push_back(h2_e14);
  v_h2.push_back(row);
  row.clear();

  row.push_back(h2_e20);
  row.push_back(h2_e21);
  row.push_back(h2_e22);
  row.push_back(h2_e23);
  row.push_back(h2_e24);
  v_h2.push_back(row);
  row.clear();

  row.push_back(h2_e30);
  row.push_back(h2_e31);
  row.push_back(h2_e32);
  row.push_back(h2_e33);
  row.push_back(h2_e34);
  v_h2.push_back(row);
  row.clear();


  gsl_histogram * histo;
  histo = gsl_histogram_alloc (nbins);
  
  double xmin=0., xmax=50.;
  double xrange[nbins+1];
  for(int i=0; i < nbins+1; ++i){
    xrange[i] =  xmin + (xmax-xmin)/nbins*i;
  }

  //gsl_histogram_set_ranges_uniform (histo, XRANGE, XRANGE_Dim);
  for (int i=0; i<4; ++i){
    for(int j=0; j<5; ++j){
      gsl_histogram_set_ranges (v_h1[i][j], xrange, nbins+1);
      gsl_histogram_set_ranges (v_h2[i][j], xrange, nbins+1);
    }
  }

  double v_meta1[4][5][2];
  double v_meta2[4][5][2];

  for (int i=0; i < 4; ++i ){
    for(int j = 0; j < 5; ++j){
      for (int p = 0; p < 2; ++p){
        v_meta1[i][j][p]=0;
        v_meta2[i][j][p]=0;
      }
    }
  }

  int k = 0;
  int l = 0;

  int nevents;
  nevents = treesize;
  nevents = 1e8;


  // *********************  MAIN LOOP  **************************

  for (int i=0; i<nevents; ++i){

		nui->GetEntry(i);

    if( ( rdet_id==0    || rdet_id==1   || rdet_id==2 ) &&
        ( rmhnl==0.5    || rmhnl==1.0   || rmhnl==1.5 || rmhnl==1.9) &&
        ( roffaxis==0   || roffaxis==10 || 
          roffaxis==20  || roffaxis==30 || roffaxis==40)
        ){
      
      if(rmhnl == 1.9){
        k = 3;
      }
      else {
        k = (rmhnl-0.5)/0.5;
      }

      l = roffaxis/10;

      //cout << k << " " << l << endl;

      if( ( rid == 12 || rid == 14 || rid == 16) ){
        gsl_histogram_increment(v_h1[k][l], re);
        v_meta1[k][l][0] = rmhnl;
        v_meta1[k][l][1] = roffaxis;
      }
      if( ( rid ==-12 || rid == -14 || rid == -16) ){
        gsl_histogram_increment(v_h2[k][l], re);
        v_meta2[k][l][0] = rmhnl;
        v_meta2[k][l][1] = roffaxis;
      }
      
    }
		
    if( (i+1) % 10000 == 0){	cout << 1.*(i+1)/nevents * 100. << "%" <<"\r";	}
	}

  ofstream data1("data1.dat"); 
  ofstream data2("data2.dat"); 

  for (int k=0; k < 4; ++k){
    for(int l=0; l < 5; ++l){
      data1 << v_meta1[k][l][0] << " " << v_meta1[k][l][1]<< " ";
      for (int i=0; i < nbins; ++i){
        data1 << gsl_histogram_get( v_h1[k][l], i)<< " ";
      }
      data1 << endl;
    }
  }

  for ( k=0; k < 4; ++k){
    for( l=0; l < 5; ++l){
      data2 << v_meta2[k][l][0] << " " << v_meta2[k][l][1]<< " ";
      for (int i=0; i < nbins; ++i){
        data2 << gsl_histogram_get( v_h2[k][l], i)<< " ";
      }
      data2 << endl;
    }
  }

  /*
  for(int i=0; i < nbins; ++i){
    odata << xrange[i] << " " << xrange[i+1] << " " << gsl_histogram_get(histo , i) << endl;
  }
  */

  return 0;
	
}

