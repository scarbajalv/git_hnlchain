// ./X ebeam mHNL.dat maxoffaxis deltaoffaxis seed nroot index
     
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
#include <iomanip>

#include "TH1.h"
#include "THStack.h"
#include "TTree.h"
#include "TFile.h"

#include <sstream>
#include <string>

#include "detect.h"
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
	
	if(argc<7){
		cout<<"FATAL ERROR: ./X ebeam mHNL.dat maxoffaxis deltaoffaxis seed nroot index"<<endl;
		exit (EXIT_FAILURE);
	}
	
	// *********************** MAIN PARAMETERS ************************
	
	int 	idGun	=	431;
	int 	nList = 1;
	bool	hasLifetime = true; ///if false, decays at origin.
	int 	idhnl  = 2000000001;
	bool 	atRest = false; /// if true, ignores energy and sets eeGun=m	
	
	
	// ****************************************************************
   
	// Read ebeam
	stringstream ebeam_ss;
	ebeam_ss << argv[1];
	string ebeam_s = ebeam_ss.str();
	double ebeam = stod(ebeam_s);

	/// Read HNL mass from $mainconfig/mHNL.dat (first line is #mHNL)
	ifstream getmass(argv[2]);
	string massstring;
	getmass >> massstring;
	massstring=massstring.substr(1,massstring.size());
	double hnlmass=stod(massstring);

	string maxoffaxiss=argv[3];
  int maxoffaxis = stoi(maxoffaxiss);

  string deltaoffaxiss=argv[4];
	int deltaoffaxis = stoi(deltaoffaxiss) ;

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

	stringstream idata_ss;
	idata_ss << getEnv("mainconfig") <<"d431data-120-" << argv[7] << ".dat";	
  string idata_s = idata_ss.str();
  char idata_c[idata_s.size()+1];
	strcpy(idata_c, idata_s.c_str());
	fstream idata(idata_s, std::ios_base::in);	
	int jEvents = countlines(idata_c);
	int nEvent = jEvents;  /// less or equal than jEvents
  
  // pythia config files
  
 	pythia.readFile("./mainconfig/pythiaconfig.ini");
	pythia.readFile(argv[2]); /// Leer el otro config file como argumento.

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
  cout<<"Ejecutando "<< argv[0] <<endl;
  cout<<"idGun = "<< idGun <<endl;
   	cout<<"index = "<<argv[7]<<endl;
   	cout<<"Ebeam = "<<ebeam<<endl;
  	cout<<"mHNL = "<<hnlmass<<endl;
  	//cout<<"Off-Axis = "<<offaxis<<endl;
  	cout<<"seed = "<<seed<<endl;
	cout<<"jEvents = "<<jEvents<<endl;
  	cout<<"nEvents = "<<nEvent<<endl<<endl;

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
      if((icount+1) % 100000 == 0){
      cout<<"\r"<<fixed<<setprecision(1)<<(icount+1)/1000000.<<" M"<<flush;
      }

		}
	}
	cout<<endl;

  
  // Variables auxiliares
  double xprod, yprod, zprod, x, y, z, px, py, pz, param;
  int pdg, mother;
  
  
  
  int noffaxis = maxoffaxis/deltaoffaxis + 1;   
  //ofstream alldata[noffaxis];
  //ofstream nudata[noffaxis];
  ofstream nudet[noffaxis];

    if (maxoffaxis % deltaoffaxis != 0){
      cout<<"FATAL ERROR: maxoffaxis debe ser múltiplo entero de deltaoffaxis"<<endl;
      exit (EXIT_FAILURE);
    }   

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
  nu.Branch("pmother",&rpmother,"rpmother/D");
  
  double unitfactor = 1000; // 1000 for mm

  // Define LarTPC detector geometry
  double w = 7.; // width X (m)
  double h = 3.; // height Y (m)
  double l = 5.; // length Z (m)
  double z0det = 574; // z0 (m)

  // Define LarTPC detector geometry
  double w_mpd = 5.; // largo del cilindro (m)
  double r_mpd = 0.5*5; // radio (m)


  double  xidet, xfdet, yidet , yfdet, zidet, zfdet;
  double  ximpd, xfmpd, rmpd , czmpd, cympd;
  int counter0=0, counter1=0, counter2=0;
  int counterhnl0=0, counterhnl1=0;



// *******************************************************************
// ************************ MAIN LOOP ********************************
   
  cout<<"Realizando cadena de decays para mHNL = "<<massstring<<" GeV"<<endl;

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) { // Main event Loop

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
    for (int i = 0; i < event.size(); ++i) { // loop over all particles
			
			if (pythia.event[i].id()==12||pythia.event[i].id()==-12
        ||pythia.event[i].id()==14||pythia.event[i].id()==-14
        ||pythia.event[i].id()==16||pythia.event[i].id()==-16){ // begin if nu
				
        double xx[3]={pythia.event[i].xProd(),pythia.event[i].yProd(),pythia.event[i].zProd()};
        double pp[3]={pythia.event[i].px(),pythia.event[i].py(),pythia.event[i].pz()};  

        int ioffaxis = 0;

        while (ioffaxis < noffaxis){  // being nu off-axis loop

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
              rmhnl = hnlmass;
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
       
        } // End nu offaxis loop			

      } // End if nu


			if (pythia.event[i].id()==idhnl||pythia.event[i].id()==-idhnl){ // if HNL

        int ioffaxis = 0;

        while (ioffaxis < noffaxis){ // begin HNL offaxis loop

          bool saveintree = false; // asumimos que no ingresa al detector

          int  offaxis = deltaoffaxis*ioffaxis; // ¡está en metros!

          xidet = (-w/2+offaxis)*unitfactor, 
          xfdet = (w/2+offaxis)*unitfactor, 
          yidet = -h/2*unitfactor, 
          yfdet = h/2*unitfactor,
          zidet = z0det*unitfactor,
          zfdet = (z0det+l)*unitfactor;

          ximpd = (-w_mpd/2+offaxis)*unitfactor;
          xfmpd = (w_mpd/2+offaxis)*unitfactor;
          czmpd = (z0det+l+rmpd)*unitfactor;
          cympd = 0;

          // HNL decay dentro de LarTPC
              
          if (pythia.event[i].xDec()>xidet && pythia.event[i].xDec()<xfdet &&
            pythia.event[i].yDec()>yidet && pythia.event[i].yDec()<yfdet &&
            pythia.event[i].zDec()>zidet && pythia.event[i].zDec()<zfdet){
              rdet_id = 0;
              saveintree = true;            
          } // en of nudet

          if (pythia.event[i].xDec()>ximpd && pythia.event[i].xDec()<xfmpd &&
          (pow(pythia.event[i].yDec()-cympd,2)+pow(pythia.event[i].zDec()-czmpd,2) 
            < pow(rmpd,2)) ){
              rdet_id = 1;
              saveintree = true;
          } // en of nudet
          
          if (saveintree == true){
              rebeam = ebeam;
              rmhnl = hnlmass;
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

        } // end HNL offaxis loop
					
			} // End if HNL
								
	} // End loop over all particles

  if((iEvent+1) % 10000 == 0){
	 cout<<"\r"<<fixed<<setprecision(2)<<(iEvent+1)/1000000.<<" M"<<flush;
	}
  
  } // End of // Main event Loop

  
  cout<<endl<<"Decays finalizados."<<endl<<endl;


  // Write Root File
  cout<<"Exporting root file..."<<endl;
  outFile.cd();
  nu.Write("",2);

	cout<<endl<<"SUCCESS!"<<endl;

  // Done.
  return 0;
  

}
