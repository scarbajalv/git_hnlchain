// ./multigunvX ebeam mHNL.dat maxoffaxis deltaoffaxis seed nroot index

// Basado en main21.cc
// Single particle gun a partir de un idata file que contiene
// vectores del tipo (e,theta,phi) de una partícula idGun.
// Imprime los parámetros de las partículas finales luego del decay de idGun.

#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h> 
#include <sys/stat.h>
#include <stdlib.h> 
#include <vector>

#include "TH1.h"
#include "THStack.h"
#include "TTree.h"
#include "TFile.h"

#include <sstream>
#include <string>
using namespace Pythia8;

#include "detect.h"

// ---------------- COUNT LINES -----------------
/// Function to count lines of file.
int countlines(char* filename){
	
	int count=0;
	string line;
	ifstream ifile;
	ifile.open(filename);	
	while (getline(ifile,line)){ /// While it is possible to store lines
		///cout<<line<<endl; /// Print Line
		count++;		/// count lines
	}
	ifile.clear(); ifile.seekg(0);
	return count;
}

// ---------------- SINGLE PARTICLE GUN FUNCTION -----------------
/// Input: flavour, energy, direction (theta, phi).
/// If theta < 0 then random choice over solid angle.
/// Optional final argument to put particle at rest => E = m.
void fillParticle(int id, double ee, double thetaIn, double phiIn,
  Event& event, ParticleData& pdt, Rndm& rndm, bool atRest = false,
  bool hasLifetime = false) {

  /// Reset event record to allow for new event.
  event.reset();
  /// Select particle mass; where relevant according to Breit-Wigner.
  double mm = pdt.mSel(id);
  double pp = sqrtpos(ee*ee - mm*mm);
  /// Special case when particle is supposed to be at rest.
  if (atRest) {
    ee = mm;
    pp = 0.;
  }
  /// Angles as input or uniform in solid angle.
  double cThe, sThe, phi;
  if (thetaIn >= 0.) {
    cThe = cos(thetaIn);
    sThe = sin(thetaIn);
    phi  = phiIn;
  } else {
    cThe = 2. * rndm.flat() - 1.;
    sThe = sqrtpos(1. - cThe * cThe);
    phi = 2. * M_PI * rndm.flat();
  }
  /// Store the particle in the event record.
  int iNew = event.append( id, 1, 0, 0, pp * sThe * cos(phi),
    pp * sThe * sin(phi), pp * cThe, ee, mm);
  /// Generate lifetime, to give decay away from primary vertex.
  if (hasLifetime) event[iNew].tau( event[iNew].tau0() * rndm.exp() );
}



//============================= PYTHIA8 ================================

int main(int argc, char *argv[]) {
	
	Pythia pythia; 
	
	if(argc<7){
		cout<<"FATAL ERROR: Se deben ingresar: ebeam mHNL.dat maxoffaxis deltaoffaxis seed nroot index"<<endl;
		exit (EXIT_FAILURE);
	}
	
	// *********************** MAIN PARAMETERS ************************
	
	stringstream idata_ss;
	idata_ss << "./mainconfig/d431data-120-" << argv[7] << ".dat";
	string idata_s = idata_ss.str();
	char idata_c[idata_s.size()+1];
	strcpy(idata_c, idata_s.c_str());

	fstream idata(idata_s, std::ios_base::in);
	
	int jEvents = countlines(idata_c);
	
	/*
	fstream idata("./mainconfig/d431data-120-1.dat", std::ios_base::in);
	
	int jEvents = countlines("./mainconfig/d431data-120-1.dat");
	*/
	
	int idGun  = 431;

		
	int nEvent = jEvents;  /// less or equal than jEvents
	
	// ****************************************************************
  
  int nList = 1;
  bool   hasLifetime =  true; ///if false, decays at origin.

  int 	 idhnl  = 2000000001;
  bool   atRest = false; /// if true, ignores energy and sets eeGun=m
  

  
 	/// pythia config files
  
 	pythia.readFile("./mainconfig/pythiaconfig.ini");
	pythia.readFile(argv[2]); /// Leer el otro config file como argumento.
	
	// Read ebeam
	stringstream ebeam_ss;
	ebeam_ss << argv[1];
	string ebeam_s = ebeam_ss.str();
	double ebeam = stod(ebeam_s);

	/// Read HNL mass from mainconfig/mHNL.dat (its first line is #mHNL)
	ifstream getmass(argv[2]);
	string massstring;
	getmass >> massstring;
	massstring=massstring.substr(1,massstring.size());
	double hnlmass=stod(massstring);

	 string maxoffaxiss=argv[3];
  	int maxoffaxis = stoi(maxoffaxiss);

  	string deltaoffaxiss=argv[4];
	int  deltaoffaxis = stoi(deltaoffaxiss) ;

	// Set random seed as 3th argument - eg 001
	string seed = argv[5];
	stringstream seedconfigss;
	seedconfigss<<"Random:seed = "<<seed;
	string seedconfig = seedconfigss.str();
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedconfig);

	stringstream nroot_ss;
	nroot_ss << ebeam_s << "-"<<argv[6]<<".root";
	string nroot_s = nroot_ss.str();
	char nroot_c[nroot_s.size()+1];
	strcpy(nroot_c, nroot_s.c_str());
	TFile outFile(nroot_c, "recreate");
  
  // Generator; shorthand for event and particleData.

  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;
  
  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Switch off ProcessLevel. Enable decays, update masses.
  pythia.readString("ProcessLevel:all = off"); 
  
  // Initialize.
  pythia.init();
  
  cout<<endl;
   	cout<<"index = "<<argv[7]<<endl;
   	cout<<"Ebeam = "<<ebeam<<endl;
  	cout<<"mHNL = "<<hnlmass<<endl;
  	//cout<<"Off-Axis = "<<offaxis<<endl;
  	cout<<"seed = "<<seed<<endl;
	cout<<"jEvents = "<<jEvents<<endl;
  	cout<<"nEvents = "<<nEvent<<endl<<endl;

  // Import Ds+ data (e,theta,phi) from idata to vector v
   // Number of columns (e, theta, phi)
  cout<<"Importando dsdata..."<<endl;
   vector <double> row;
   vector < vector <double> > v;
   double a;
   int jcount = 0, icount = 0;

	while (idata>>a){
		jcount++;
		row.push_back(a);
		if (jcount % 3 == 0){
			v.push_back(row);
			row.clear();
			icount++;
			cout<<"\r"<<setprecision(2)<<(icount)/1000000.<<" M";

		}
	}
	cout<<endl;

  
  // Variables auxiliares
  double xprod, yprod, zprod, x, y, z, px, py, pz, param;
  int pdg, mother;
  
  // Vectores para almacenar data
  vector < vector <double> > hnlallvector;
  vector < vector <double> > hnldetvector;
  vector < vector <double> > allvector;
  vector < vector <double> > nuallvector;
  vector < vector <double> > nudetvector;
  
  // ****************** EVENT LOOP ********************
 
  cout<<"Realizando cadena de decays para mHNL = "<<massstring<<" GeV"<<endl;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

  	/// Set up single particle (id, energy, theta, phi, ...)
    fillParticle
    (idGun, v[iEvent][0], v[iEvent][1], v[iEvent][2], event, pdt, pythia.rndm, atRest, hasLifetime);

    /// Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << "WFT! Event generation aborted prematurely, owing to error!\n";
      break;
    }

    /// List first few events.
    //if (iEvent < nList) {
    //  event.list();
    //}
    //

    // Loop over all particles (analysis).
    for (int i = 0; i < event.size(); ++i) {
			
			if (pythia.event[i].id()==12
				||pythia.event[i].id()==14
				||pythia.event[i].id()==16){
				
				row.clear();
				row.push_back(pythia.event[i].id());
				row.push_back(pythia.event[i].tProd());
				row.push_back(pythia.event[i].xProd());
				row.push_back(pythia.event[i].yProd());
				row.push_back(pythia.event[i].zProd());
				row.push_back(pythia.event[i].tDec());
				row.push_back(pythia.event[i].xDec());
				row.push_back(pythia.event[i].yDec());
				row.push_back(pythia.event[i].zDec());
				row.push_back(pythia.event[i].e());
				row.push_back(pythia.event[i].px());
				row.push_back(pythia.event[i].py());
				row.push_back(pythia.event[i].pz());
				row.push_back(pythia.event[i].pT());
				row.push_back(pythia.event[i].theta());
				row.push_back(pythia.event[i].phi());
				row.push_back(pythia.event[i].y());
				row.push_back(pythia.event[i].eta());
				row.push_back(pythia.event[i].index());
				row.push_back(pythia.event[i].mother1());
					
				nuallvector.push_back(row);
					
			} // End of nuall

			if (pythia.event[i].id()==idhnl||pythia.event[i].id()==-idhnl){
				
				row.clear();
				row.push_back(pythia.event[i].id());
				row.push_back(pythia.event[i].tProd());
				row.push_back(pythia.event[i].xProd());
				row.push_back(pythia.event[i].yProd());
				row.push_back(pythia.event[i].zProd());
				row.push_back(pythia.event[i].tDec());
				row.push_back(pythia.event[i].xDec());
				row.push_back(pythia.event[i].yDec());
				row.push_back(pythia.event[i].zDec());
				row.push_back(pythia.event[i].e());
				row.push_back(pythia.event[i].px());
				row.push_back(pythia.event[i].py());
				row.push_back(pythia.event[i].pz());
				row.push_back(pythia.event[i].pT());
				row.push_back(pythia.event[i].theta());
				row.push_back(pythia.event[i].phi());
				row.push_back(pythia.event[i].y());
				row.push_back(pythia.event[i].eta());
				row.push_back(pythia.event[i].index());
				row.push_back(pythia.event[i].mother1());
					
				hnlallvector.push_back(row);
					
			} // End of nuall
								
	} // End of analysis

	cout<<"\r"<<setprecision(2)<<(iEvent+1)/1000000.<<" M";
	
  } // End of event loop.  
  cout<<endl<<"Decays finalizados."<<endl<<endl;


  // ********************************************************************
  // ************************ DATA ANALYSIS *****************************
  // ********************************************************************



	int noffaxis = maxoffaxis/deltaoffaxis + 1; 	
	//ofstream alldata[noffaxis];
	//ofstream nudata[noffaxis];
	ofstream nudet[noffaxis];

  	if (maxoffaxis % deltaoffaxis != 0){
 			cout<<"FATAL ERROR: maxoffaxis debe ser múltiplo entero de deltaoffaxis"<<endl;
 			exit (EXIT_FAILURE);
 		} 	


  	//cout<<"Total Particles: "<<allvector.size()<<endl;
  	cout<<"Total Neutrinos: "<<nuallvector.size()<<endl;
  	cout<<"Total HNL: "<<hnlallvector.size()<<endl<<endl;
  	cout<<endl<<"Iniciando análisis..."<<endl<<endl;


 	 	stringstream metadatass;

 	// Create ROOT variables
 	double rebeam, rmhnl, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
	rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rpmother;
  	int mother1, mother2, auxpdg, rdet_id;

	TTree nu("nu","nu");
	nu.Branch("ebeam",&rebeam,"rebeam/D");
	nu.Branch("mhnl",&rmhnl,"rmhnl/D");
	nu.Branch("offaxis",&roffaxis,"roffaxis/D");
	nu.Branch("det_id",&rdet_id,"rdet_id/I");
	nu.Branch("id",&rid,"rid/D");
	nu.Branch("tProd",&rtProd,"rtProd/D");
	nu.Branch("xProd",&rxProd,"rxProd/D");
	nu.Branch("yProd",&ryProd,"ryProd/D");
	nu.Branch("zProd",&rzProd,"rzProd/D");
	nu.Branch("tDec",&rtDec,"rtDec/D");
	nu.Branch("xDec",&rxDec,"rxDec/D");
	nu.Branch("yDec",&ryDec,"ryDec/D");
	nu.Branch("zDec",&rzDec,"rzDec/D");
	nu.Branch("e",&re,"re/D");
	nu.Branch("px",&rpx,"rpx/D");
	nu.Branch("py",&rpy,"rpy/D");
	nu.Branch("pz",&rpz,"rpz/D");
	nu.Branch("pT",&rpT,"rpT/D");
	nu.Branch("theta",&rtheta,"rtheta/D");
	nu.Branch("phi",&rphi,"rphi/D");
	nu.Branch("y",&ry,"ry/D");
	nu.Branch("eta",&reta,"reta/D");
	nu.Branch("pindex",&rpindex,"rpindex/D");
	nu.Branch("pmother",&rpmother,"rpmother/D");

	// Define detector geometry
	double w = 7.; // width X (m)
	double h = 3.; // height Y (m)
	double l = 5.; // length Z (m)
	double z0det = 574; // z0 (m)
	double unitfactor = 1000; // 1000 for mm

	double	xidet, xfdet, yidet , yfdet, zidet,	zfdet;
	double	ximpd, xfmpd, rmpd , czmpd, cympd;
	
// ******************************************************************
// ******************** BEGIN DATA ANALYSIS *************************	

	int ioffaxis = 0;

 	while (ioffaxis < noffaxis){

 		row.clear();
 		nudetvector.clear();
 		hnldetvector.clear();

	 	int  offaxis = deltaoffaxis*ioffaxis;

		xidet = (-w/2+offaxis)*unitfactor, 
		xfdet = (w/2+offaxis)*unitfactor, 
		yidet = -h/2*unitfactor, 
		yfdet = h/2*unitfactor,
		zidet = z0det*unitfactor,
		zfdet = (z0det+l)*unitfactor;

		ximpd = (-w/2+offaxis)*unitfactor;
		xfmpd = (w/2+offaxis)*unitfactor;
		rmpd = 2.5*unitfactor;
		czmpd = (z0det+l+rmpd)*unitfactor;
		cympd = 0;


		// Seleccionar neutrinos que ingresan al detector
		for(int i=0; i<nuallvector.size(); ++i){
			double xx[3]={nuallvector[i][2],nuallvector[i][3],nuallvector[i][4]};
			double pp[3]={nuallvector[i][10],nuallvector[i][11],nuallvector[i][12]};				
			// Atraviesa el LArTPC y el MPD
			if (detect(xx,pp,offaxis)&&detectmpd(xx,pp,offaxis)){				
				row.clear();
				rdet_id = 2;
				for (int j = 0; j < 20; ++j){
					row.push_back(nuallvector[i][j]);
				}
				row.push_back(offaxis); // j=20
				row.push_back(rdet_id); // j=21
				nudetvector.push_back(row);
			} // en of nudet
			// Atraviesa solo el LArTPC
			if (detect(xx,pp,offaxis)&&!detectmpd(xx,pp,offaxis)){				
				row.clear();
				rdet_id = 0;
				for (int j = 0; j < 20; ++j){
					row.push_back(nuallvector[i][j]);
				}
				row.push_back(offaxis); // j=20
				row.push_back(rdet_id); // j=21
				nudetvector.push_back(row);
			} // en of nudet
			// Atraviesa solo el MPD
			if (detectmpd(xx,pp,offaxis)&!detect(xx,pp,offaxis)){
				row.clear();
				rdet_id = 1;
				for (int j = 0; j < 20; ++j){
					row.push_back(nuallvector[i][j]);
				}
				row.push_back(offaxis); // j=20
				row.push_back(rdet_id); // j=21
				nudetvector.push_back(row);
			} // en of nudet
		}

		// HNL decay dentro de LarTPC
		for(int i=0; i<hnlallvector.size(); ++i){			
			if (hnlallvector[i][6]>xidet && hnlallvector[i][6]<xfdet &&
				hnlallvector[i][7]>yidet && hnlallvector[i][7]<yfdet &&
				hnlallvector[i][8]>zidet && hnlallvector[i][8]<zfdet){
				row.clear();
				rdet_id = 0;
				for (int j = 0; j < 20; ++j){
					row.push_back(hnlallvector[i][j]);
				}
				row.push_back(offaxis);
				row.push_back(rdet_id); // j=21
				hnldetvector.push_back(row);
			} // en of nudet
		}

		// HNL decay dentro de MPD
		for(int i=0; i<hnlallvector.size(); ++i){
			if (hnlallvector[i][6]>ximpd && hnlallvector[i][6]<xfmpd &&
				(pow(hnlallvector[i][7]-cympd,2)+pow(hnlallvector[i][8]-czmpd,2) 
					< pow(rmpd,2)) ){
				row.clear();
				rdet_id = 1;
				for (int j = 0; j < 20; ++j){
					row.push_back(hnlallvector[i][j]);
				}
				row.push_back(offaxis);
				row.push_back(rdet_id); // j=21
				hnldetvector.push_back(row);
			} // en of nudet
		}

		// Fill TTree
	  	
		for (int i=0; i<nudetvector.size(); ++i){			
			rebeam = ebeam;
			rmhnl = hnlmass;
			roffaxis = nudetvector[i][20];
			rdet_id = nudetvector[i][21];
			rid = nudetvector[i][0];			
			rtProd = nudetvector[i][1];
			rxProd = nudetvector[i][2]; 
			ryProd = nudetvector[i][3]; 
			rzProd = nudetvector[i][4];
			rtDec = nudetvector[i][5]; 
			rxDec = nudetvector[i][6];
			ryDec = nudetvector[i][7];
			rzDec = nudetvector[i][8];
			re = nudetvector[i][9];
			rpx = nudetvector[i][10];
			rpy = nudetvector[i][11];
			rpz = nudetvector[i][12];
			rpT = nudetvector[i][13];
			rtheta = nudetvector[i][14];
			rphi = nudetvector[i][15];
			ry = nudetvector[i][16];
			reta = nudetvector[i][17];
			rpindex = nudetvector[i][18];
			rpmother = nudetvector[i][19];
			nu.Fill();
		}

		for (int i=0; i<hnldetvector.size(); ++i){			
			rebeam = ebeam;
			rmhnl = hnlmass;
			roffaxis = hnldetvector[i][20];
			rdet_id = nudetvector[i][21];
			rid = hnldetvector[i][0];
			rtProd = hnldetvector[i][1];
			rxProd = hnldetvector[i][2]; 
			ryProd = hnldetvector[i][3]; 
			rzProd = hnldetvector[i][4];
			rtDec = hnldetvector[i][5]; 
			rxDec = hnldetvector[i][6];
			ryDec = hnldetvector[i][7];
			rzDec = hnldetvector[i][8];
			re = hnldetvector[i][9];
			rpx = hnldetvector[i][10];
			rpy = hnldetvector[i][11];
			rpz = hnldetvector[i][12];
			rpT = hnldetvector[i][13];
			rtheta = hnldetvector[i][14];
			rphi = hnldetvector[i][15];
			ry = hnldetvector[i][16];
			reta = hnldetvector[i][17];
			rpindex = hnldetvector[i][18];
			rpmother = hnldetvector[i][19];
			nu.Fill();
		}

		//cout<<nuallvector[j][0]<<endl;
		cout	<<"off-axis = "<<offaxis<<" => "
					<<"Detected Neutrinos: "<<nudetvector.size()<<endl
					<<"HNL decay in detector: "<<hnldetvector.size()<<endl<<endl;

		ioffaxis = ioffaxis + 1;

	} // End of offaxis loop
	

	// Write Root File
  	cout<<"Exporting root file..."<<endl;
  	outFile.cd();
  	nu.Write("",2);
   
	cout<<endl<<"SUCCESS!"<<endl;

  // Done.
  return 0;
  

}
