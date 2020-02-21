// ./root-merge ebeam nrooti

// Añade la data del archivo nrooti.root al archivo ebeam.root

// Agregar data de trees de un root "IN" a trees de un root "OUT".
// El root "OUT" ya debe contener trees a los cuales se agregará la data.
// Se debe conocer los nombres de los trees de interés en ambos root.

#include <iostream>
#include <fstream>


#include <TTree.h> 
#include <TFile.h>
#include <TH1.h>

#include <string>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]){
	
	// *******************************************************
	// Check ENERGY and SEED
	// OUTPUT no puede estar vacío!! (Debe ser inicializado con mergeini)

	stringstream ifile_ss;
	ifile_ss << argv[1] << "-" << argv[2] <<".root";
	string ifile_s = ifile_ss.str();
	char ifile_c[ifile_s.size()+1];
	strcpy(ifile_c, ifile_s.c_str());
	
	ifstream ifiletest;
	ifiletest.open(ifile_s);
	//cout << ifiletest.good();
	
	if (ifiletest.good()==0){
		cout << endl << ifile_s <<" NOT FOUND"<<endl;
		exit(EXIT_FAILURE);
	}
	
	TFile f1(ifile_c, "open");

	stringstream ofile_ss;
	ofile_ss << argv[1]  <<".root";
	string ofile_s = ofile_ss.str();
	char ofile_c[ofile_s.size()+1];
	strcpy(ofile_c, ofile_s.c_str());
	
	TFile f2(ofile_c,"update"); // ¡Solo Para actualizar!
	
	// *******************************************************
	
	double var_e, var_theta, var_pT;

	double rebeam, rmhl, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re, rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rpmother, rnevents, 
		rnuall, rnudet;

	
	// DEFINE TREES FROM INPUT FILE:
	
	f1.cd();
	
	TTree* nui;
		
	f1.GetObject("nu",nui);

	nui->SetBranchAddress("ebeam",&rebeam);
	nui->SetBranchAddress("mhnl",&rmhl);
	nui->SetBranchAddress("offaxis",&roffaxis);
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
	nui->SetBranchAddress("pindex",&rpindex);
	nui->SetBranchAddress("pmother",&rpmother);

		//DEFINE TREES FROM OUTPUT FILE
	
	f2.cd();
	
	TTree* nuo;
		
	f2.GetObject("nu",nuo);

	nuo->SetBranchAddress("ebeam",&rebeam);
	nuo->SetBranchAddress("mhnl",&rmhl);
	nuo->SetBranchAddress("offaxis",&roffaxis);
	nuo->SetBranchAddress("id",&rid);
	nuo->SetBranchAddress("tProd",&rtProd);
	nuo->SetBranchAddress("xProd",&rxProd);
	nuo->SetBranchAddress("yProd",&ryProd);
	nuo->SetBranchAddress("zProd",&rzProd);
	nuo->SetBranchAddress("tDec",&rtDec);
	nuo->SetBranchAddress("xDec",&rxDec);
	nuo->SetBranchAddress("yDec",&ryDec);
	nuo->SetBranchAddress("zDec",&rzDec);
	nuo->SetBranchAddress("e",&re);
	nuo->SetBranchAddress("px",&rpx);
	nuo->SetBranchAddress("py",&rpy);
	nuo->SetBranchAddress("pz",&rpz);
	nuo->SetBranchAddress("pT",&rpT);
	nuo->SetBranchAddress("theta",&rtheta);
	nuo->SetBranchAddress("phi",&rphi);
	nuo->SetBranchAddress("y",&ry);
	nuo->SetBranchAddress("eta",&reta);
	nuo->SetBranchAddress("pindex",&rpindex);
	nuo->SetBranchAddress("pmother",&rpmother);

	
			
	// UPDATE TREES
	
	for (int i=0; i<nui->GetEntries();++i){
		nui->GetEntry(i); //almacenar e(i) en var1
		nuo->Fill(); // almacenar var1 en tDp
	}

	
		
	//Write trees (removing history).
	f2.cd();
	
	nuo->Write("",2);

	return 0;


	
}
