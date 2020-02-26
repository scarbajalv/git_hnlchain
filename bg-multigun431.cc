  
// ./X ebeam maxoffaxis deltaoffaxis seed nroot index

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

#include "detect.h"

#include <sstream>
#include <string>

#include "cpp-getEnv.h"

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

	if(argc<6){
		cout<<"FATAL ERROR: ./X ebeam maxoffaxis deltaoffaxis seed nroot index"<<endl;
		exit (EXIT_FAILURE);
	}
	
	// *********************** MAIN PARAMETERS ************************
	
	int 	idGun  = 431;
  int   nList = 1;
  bool  hasLifetime =  true; ///if false, decays at origin.
  int   idhnl  = 2000000001;
  bool  atRest = false; /// if true, ignores energy and sets eeGun=m  
	
	// ****************************************************************
  	
	stringstream ebeam_ss;
	ebeam_ss << argv[1];
	string ebeam_s = ebeam_ss.str();
	double ebeam = stod(ebeam_s);

	string maxoffaxiss=argv[2];
	int maxoffaxis = stoi(maxoffaxiss);

	string deltaoffaxiss=argv[3];
	int  deltaoffaxis = stoi(deltaoffaxiss) ; /// offaxis (metros)
	
	string seed = argv[4];
	stringstream seedconfigss;
	seedconfigss<<"Random:seed = "<<seed;
	string seedconfig = seedconfigss.str();
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedconfig);

	stringstream nroot_ss;
	nroot_ss << ebeam_s << "-"<<argv[5]<<".root";
	string nroot_s = nroot_ss.str();
	char nroot_c[nroot_s.size()+1];
	strcpy(nroot_c, nroot_s.c_str());
	TFile outFile(nroot_c, "recreate");

	stringstream idata_ss;
	idata_ss << getEnv("mainconfig") <<"d431data-120-" << argv[6] << ".dat"; 
	string idata_s = idata_ss.str();
	char idata_c[idata_s.size()+1];
	strcpy(idata_c, idata_s.c_str());
	fstream idata(idata_s, std::ios_base::in);	
	int jEvents = countlines(idata_c);
	int nEvent = jEvents;  /// less or equal than jEvents
	
	// ****************************************************************
    
 	/// pythia config files
  
 	pythia.readFile("./mainconfig/pythiaconfig-bg.ini");
    
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
  cout<<"Importando "<< idata_s <<endl;
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
  idata.close();

  
  // Variables auxiliares
  double xprod, yprod, zprod, x, y, z, px, py, pz, param;
  int pdg, mother;
  
  int noffaxis = maxoffaxis/deltaoffaxis + 1;   
  
    if (maxoffaxis % deltaoffaxis != 0){
      cout<<"FATAL ERROR: maxoffaxis debe ser múltiplo entero de deltaoffaxis"<<endl;
      exit (EXIT_FAILURE);
    }   

 	cout<<endl<<"Exporting data..."<<endl;
 	cout<<endl;
  	//cout<<"Off-Axis = "<<offaxis<<endl;
  cout<<"seed = "<<seed<<endl;
  cout<<"nevents = "<<nEvent<<endl;

	double rebeam, roffaxis, rid, rtProd, rxProd, ryProd, rzProd, rtDec, rxDec, ryDec, rzDec, re,
	rpx, rpy, rpz, rpT, rtheta, rphi, ry, reta, rpindex, rpmother;
  int mother1, mother2, auxpdg, rdet_id;

	TTree nu("nu","nu");
	nu.Branch("ebeam",&rebeam,"rebeam/D");
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

	int counter0=0, counter1=0, counter2=0;

  // *******************************************************************
  // ************************ MAIN LOOP ********************************
  
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

    
    for (int i = 0; i < event.size(); ++i) {// Loop over all particles (analysis).
			
			if (pythia.event[i].id()==12||pythia.event[i].id()==-12
				||pythia.event[i].id()==14||pythia.event[i].id()==-14
				||pythia.event[i].id()==16||pythia.event[i].id()==-16){ // if nu

          double xx[3]={pythia.event[i].xProd(),pythia.event[i].yProd(),pythia.event[i].zProd()};
          double pp[3]={pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz()};  

        int ioffaxis = 0;

        while (ioffaxis < noffaxis){  // being off-axis loop

          int  offaxis = deltaoffaxis*ioffaxis;
            
            bool saveintree = false; // asumimos que no ingresa al detector

            // Atraviesa el LArTPC y el MPD
            if (detect(xx,pp,offaxis)&&detectmpd(xx,pp,offaxis)){       
              rdet_id = 2;
              saveintree = true;          
            } // en of nudet
            // Atraviesa solo el LArTPC
            if (detect(xx,pp,offaxis)&&!detectmpd(xx,pp,offaxis)){        
              rdet_id = 0;
              saveintree = true;            
            } // en of nudet
            // Atraviesa solo el MPD
            if (detectmpd(xx,pp,offaxis)&&!detect(xx,pp,offaxis)){
              rdet_id = 1;
              saveintree = true;        
            } // en of nudet

            if (saveintree == true){
              rebeam = ebeam;
              roffaxis = offaxis;
              rid = pythia.event[i].id();
              rtProd = pythia.event[i].tProd();
              rxProd = pythia.event[i].xProd(); 
              ryProd = pythia.event[i].yProd(); 
              rzProd = pythia.event[i].zProd();
              rtDec = pythia.event[i].tDec(); 
              rxDec = pythia.event[i].xDec();
              ryDec = pythia.event[i].yDec();
              rzDec = pythia.event[i].zDec();
              re = pythia.event[i].e();
              rpx = pythia.event[i].px();
              rpy = pythia.event[i].py();
              rpz = pythia.event[i].pz();
              rpT = pythia.event[i].pT();
              rtheta = pythia.event[i].theta();
              rphi = pythia.event[i].phi();
              ry = pythia.event[i].y();
              reta = pythia.event[i].eta();
              rpmother = pythia.event[pythia.event[i].mother1()].id();
              nu.Fill();
            }

          ioffaxis = ioffaxis + 1;
       
        } // End of offaxis loop
				
		
			} // End of if nu.
								
    } // End of Loop over all particles (analysis)

	cout<<"\r"<<setprecision(2)<<(iEvent+1)/1000000.<<" M";
	
  } // End of event loop.  
  cout<<endl;
 
	// Write Root File
  	cout<<"Exporting root file..."<<endl;
  	outFile.cd();
  	nu.Write("",2);


	cout<<endl<<"SUCCESS!"<<endl;

  // Done.
  return 0;


}
