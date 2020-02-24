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
	
	// *********************** MAIN PARAMETERS ************************
	
	TFile outFile("001.root", "recreate");

	fstream idata("./mainconfig/dsdata-120-01.dat", std::ios_base::in);
	
	int jEvents = countlines("./mainconfig/dsdata-120-01.dat");
	
	int idGun  = 431;
	
	
	
	// Set random seed as 1sy argument - eg 001
	string seed = argv[1];
	stringstream seedconfigss;
	seedconfigss<<"Random:seed = "<<seed;
	string seedconfig = seedconfigss.str();
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedconfig);
	

	int nEvent = jEvents;  /// less or equal than jEvents
	
	// ****************************************************************
  
  int 	nList = 1;
  bool   hasLifetime =  true; ///if false, decays at origin.

  /*
  int 	 idhnl  = 2000000001;
  */
  bool   atRest = false; /// if true, ignores energy and sets eeGun=m
  

  
 	/// pythia config files
  
 	pythia.readFile("./mainconfig/pythiaconfig.ini");
  /*
	pythia.readFile(argv[1]); /// Leer el otro config file como argumento.
	
	/// Read HNL mass
	ifstream getmass(argv[1]);
	string massstring;
	getmass >> massstring;
	massstring=massstring.substr(1,massstring.size());
	double hnlmass=stod(massstring);
	*/

  
  // Generator; shorthand for event and particleData.

  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;
  
  // Switch off automatic event listing in favour of manual.
  ///pythia.readString("Random:seed = 003");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Switch off ProcessLevel. Enable decays, update masses.
  pythia.readString("ProcessLevel:all = off"); 
  
  // Initialize.
  pythia.init();
  
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
  vector < vector<double> > allvector;
  vector < vector<double> > nuallvector;
  vector < vector<double> > nudetvector;
  
  // ****************** EVENT LOOP ********************
  
  cout<<"Realizando cadena de decays de BG..."<<endl;

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
    /*if (iEvent < nList) {
      event.list();
    }
    */

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
				row.push_back(pythia.event[pythia.event[i].mother1()].id());
					
				nuallvector.push_back(row);
					
			} // End of nuall
								
	} // End of Loop over all particles (analysis)

	cout<<"\r"<<setprecision(2)<<(iEvent+1)/1000000.<<" M";
	
  } // End of event loop.  
  cout<<endl;


  // ************************************************************************
  // ************************ PRE DATA ANALYSIS *****************************

  	string maxoffaxiss=argv[2];
  	int maxoffaxis = stoi(maxoffaxiss);

  	string deltaoffaxiss=argv[3];
	int  deltaoffaxis = stoi(deltaoffaxiss) ; /// offaxis (metros)

	int noffaxis = maxoffaxis/deltaoffaxis + 1; 	
	//ofstream alldata[noffaxis];
	//ofstream nudata[noffaxis];
	//ofstream nudet[noffaxis];

  	if (maxoffaxis % deltaoffaxis != 0){
 			cout<<"FATAL ERROR: maxoffaxis debe ser múltiplo entero de deltaoffaxis"<<endl;
 			exit (EXIT_FAILURE);
 		} 	

 	cout<<endl<<"Exporting data..."<<endl;
 	cout<<endl;
  	//cout<<"Off-Axis = "<<offaxis<<endl;
  	cout<<"seed = "<<seed<<endl;
  	cout<<"nevents = "<<nEvent<<endl;
  	cout<<"Total Particles: "<<allvector.size()<<endl;
  	cout<<"Total Neutrinos: "<<nuallvector.size()<<endl<<endl;


 	int ioffaxis = 0;
 	//stringstream metadatass;
 	
 	// Create folder ./seed to store data
 	/*
 	stringstream seedfolderss;
 	seedfolderss<<"./"<<seed;
 	std::string seedfolders=seedfolderss.str();
 	char seedfolderchar[seedfolders.length()+1];
 	strcpy(seedfolderchar, seedfolders.c_str());  
 	mkdir(seedfolderchar, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	*/
 	// Create ROOT variables
 	double roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
	rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rpmother;
  	int mother1, mother2, auxpdg;

	TTree nu("nu","nu");
	nu.Branch("offaxis",&roffaxis,"roffaxis/D");
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

	double rnevents = nEvent, rnuall, rnudet;
	TTree meta("meta","meta");
	meta.Branch("offaxis",&roffaxis,"roffaxis/D");
	meta.Branch("nevents",&rnevents,"rnevents/D");
	meta.Branch("nuall",&rnuall,"rnuall/D");
	meta.Branch("nudet",&rnudet,"rnudet/D");


	// ******************************************************************
	// ******************** BEGIN DATA ANALYSIS *************************


 	while (ioffaxis < noffaxis){

 		row.clear();
 		nudetvector.clear();

	 	int  offaxis = deltaoffaxis*ioffaxis;

	 	// Detector geometry (fiducial): 3m(W) x 2m(H)
	 	double 	xidet=-1500+offaxis*1000., 
				xfdet=1500+offaxis*1000., 
				yidet=-1000, 
				yfdet=1000, 
				zdet=574000;

		for(int i=0; i<nuallvector.size(); ++i){
			// Seleccionar neutrinos que ingresan al detector
			xprod=nuallvector[i][2];
			yprod=nuallvector[i][3];
			zprod=nuallvector[i][4];
			px=nuallvector[i][10];
			py=nuallvector[i][11];
			pz=nuallvector[i][12];
				
			param=(zdet-zprod)/pz;
			
			x=xprod+param*px;
			y=yprod+param*py;
			z=zprod+param*pz;
				
			if (x>xidet && x<xfdet && y>yidet && y<yfdet && param>0){
				
				row.clear();
				for (int j = 0; j < 20; ++j){
					row.push_back(nuallvector[i][j]);
				}

				nudetvector.push_back(row);
			} // en of nudet
			
			//cout<<nuallvector[j][0]<<endl;

		}

		//  DEFINE OUTPUT *.dat FILENAMES 
		/*
		string offaxiss = std::to_string(offaxis);
	
		std::stringstream nudetss;			
	 	nudetss<<"./"<<seed<<"/datanudet-"<<offaxiss<<".dat";
		std::string nudets = nudetss.str();
	  	nudet[ioffaxis].open(nudets);
		*/
		//std::stringstream alldatass;
		//std::stringstream nudatass;
		//alldatass<<"dataall-"<<massstring<<"-"<<offaxiss<<".dat";
		//nudatass<<"datanuall-"<<massstring<<"-"<<offaxiss<<".dat";
		//std::string alldatas = alldatass.str();
		//std::string nudatas = nudatass.str();
	  	//alldata[ioffaxis].open(alldatas);
	  	//nudata[ioffaxis].open(nudatas);


	  // EXPORT DATA TO *.dat and ROOT FILEs
	  
	  	/**
	  	alldata[ioffaxis]<<hnlmass<<" "<<offaxis<<" "<<allvector.size()<<" "<<nuallvector.size()<<" "<<nudetvector.size()<<endl;
	  	for (int i=0; i<allvector.size(); ++i){
			for (int j=0; j<20; ++j){
			alldata[ioffaxis]<<setprecision(12)<<allvector[i][j]<<" ";
			}
			alldata[ioffaxis]<<endl;
		}		
		nudata[ioffaxis]<<hnlmass<<" "<<offaxis<<" "<<allvector.size()<<" "<<nuallvector.size()<<" "<<nudetvector.size()<<endl;
		for (int i=0; i<nuallvector.size(); ++i){
			for (int j=0; j<20; ++j){
			nudata[ioffaxis]<<setprecision(12)<<nuallvector[i][j]<<" ";
			}
			nudata[ioffaxis]<<endl;
		}
		**/		
		//nudet[ioffaxis]<<" "<<offaxis<<" "<<allvector.size()<<" "<<nuallvector.size()<<" "<<nudetvector.size()<<endl;
		
		// Fill meta TTree
		roffaxis = offaxis;
		rnuall = nuallvector.size();
		rnudet = nudetvector.size();
		meta.Fill();

		for (int i=0; i<nudetvector.size(); ++i){
			// Fill *.dat file
			/*
			for (int j=0; j<20; ++j){			
				nudet[ioffaxis]<<setprecision(12)<<nudetvector[i][j]<<" ";
			}
			nudet[ioffaxis]<<endl;
			*/

			// Fill nu TTree
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

		cout<<"off-axis = "<<offaxis<<" => "
		<<"Detected Neutrinos: "<<nudetvector.size()<<endl<<endl;

		// add info to stringstream to be exported to metadata-seed.dat
		//stringstream metadatass;
		/*
		metadatass<<offaxiss<<".dat"
		<<" "<<offaxiss		
		<<" "<<nEvent
		<<" "<<nuallvector.size()
		<<" "<<nudetvector.size()	
		<<endl;
		*/

		ioffaxis = ioffaxis + 1;
 
	}
	// ************************ END DATA ANALYSIS *************************
	// ********************************************************************


  // EXPORT METADATA
  /*
  char data[]="metadata.dat"; /// data file to update
	char auxdata[]="auxdata.dat"; /// aux file
	ifstream ifile(data);
  ofstream ofile(auxdata); 
      
  /// If file exists, update it.
  if (ifile){
		ofile<<ifile.rdbuf();
		
		ofile<<massstring<<"-"<<offaxiss<<".dat"
		<<" "<<massstring
		<<" "<<offaxiss		
		<<" "<<nEvent
		<<" "<<nuallvector.size()
		<<" "<<nudetvector.size()	
		<<endl;
		
		remove(data);
		rename(auxdata,data);
	}
	/// If not, create it.
	else{

		ofile<<massstring<<"-"<<offaxiss<<".dat"
		<<" "<<massstring
		<<" "<<offaxiss
		<<" "<<nEvent
		<<" "<<nuallvector.size()
		<<" "<<nudetvector.size()
		<<endl;
		
		rename(auxdata,data);
	}

	*/
  	

//  ********************** EXPORT METADATA *************************
 
	// create string to be exported to metadata-seed.dat
	/*
	std::string metadatas = metadatass.str();
	cout<<metadatas<<endl;

	// Create char with content "metadata-seed.dat"
	std:: stringstream metadatafnamess;
	metadatafnamess<<"./"<<seed<<"/metadata-"<<seed<<".dat";
	std::string metadatafname=metadatafnamess.str();
	char data[metadatafname.length()+1];
	strcpy(data, metadatafname.c_str()); 

	char auxdata[]="auxdata.dat"; /// aux file
	ifstream ifile(data);
	ofstream ofile(auxdata); 
	

  
      
  /// If file exists, update it.
  if (ifile){
		ofile<<ifile.rdbuf();
		
		ofile<<metadatas;
		
		remove(data);
		rename(auxdata,data);
	}
	/// If not, create it.
	else{

		ofile<<metadatas;
		
		rename(auxdata,data);
	}
  	*/

	// Write Root File
  	cout<<"Exporting root file..."<<endl;
  	outFile.cd();
  	nu.Write("",2);
  	meta.Write("",2);


	cout<<endl<<"SUCCESS!"<<endl;

  // Done.
  return 0;


}
