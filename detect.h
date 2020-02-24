#include <vector>
#include <cmath> 

using namespace std; 

int detect(double x[3], double p[3], double offaxis){
	
	double w = 7.; // width X (m)
	double h = 3.; // height Y (m)
	double l = 5.; // length Z (m)
	double z0det = 574; // z0 (m)
	double unitfactor = 1000; // 1000 for mm

	double	xidet = (-w/2+offaxis)*unitfactor, 
					xfdet = (w/2+offaxis)*unitfactor, 
					yidet = -h/2*unitfactor, 
					yfdet = h/2*unitfactor,
					zidet = z0det*unitfactor,
					zfdet = (z0det+l)*unitfactor;
	
	double omega[6][3];
	double eta[6][3];

	for (int i=0; i<3; ++i){
		omega[i][0] = xidet;
		omega[i][1] = yidet;
		omega[i][2] = zidet;
	}
	for (int i=3; i<6; ++i){
		omega[i][0] = xfdet;
		omega[i][1] = yfdet;
		omega[i][2] = zfdet;
	}

	eta[0][0]=0; 	eta[0][1]=0; 	eta[0][2]=-1;
	eta[1][0]=-1; eta[1][1]=0;	eta[1][2]=0;
	eta[2][0]=0;	eta[2][1]=1;	eta[2][2]=0;
	eta[3][0]=1;	eta[3][1]=0;	eta[3][2]=0;
	eta[4][0]=0;	eta[4][1]=0;	eta[4][2]=1;
	eta[5][0]=0;	eta[5][1]=-1;	eta[5][2]=0;

	bool detected=false; 
	vector  <vector <double>> points;
	vector <double> row;

	double xx[3];

	for (int i=0; i<6; ++i){
		for (int j=0; j<3; ++j){
			xx[j]	=	x[j]+p[j]
							*((omega[i][0]-x[0])*eta[i][0]
								+(omega[i][1]-x[1])*eta[i][1]
								+(omega[i][2]-x[2])*eta[i][2])/
							(p[0]*eta[i][0]+p[1]*eta[i][1]+p[2]*eta[i][2]);
		}
		//cout <<"point: "<< xx[0]<< " " << xx[1] <<" " << xx[2]<<endl;
		if(i==0||i==4){
			if(xx[0]>xidet && xx[0]<xfdet && xx[1]>yidet && xx[1]<yfdet){
				detected = true;
				row.clear();
				row.push_back(xx[0]);
				row.push_back(xx[1]);
				row.push_back(xx[2]);
				points.push_back(row);
			}
		}
		if(i==1||i==3){
			if(xx[1]>yidet && xx[1]<yfdet && xx[2]>zidet && xx[2]<zfdet ){
				detected = true;
				row.clear();
				row.push_back(xx[0]);
				row.push_back(xx[1]);
				row.push_back(xx[2]);
				points.push_back(row);
			}
		}
		if(i==2||i==5){
			if(xx[0]>xidet && xx[0]<xfdet && xx[2]>zidet && xx[2]<zfdet){
				detected = true;
				row.clear();
				row.push_back(xx[0]);
				row.push_back(xx[1]);
				row.push_back(xx[2]);
				points.push_back(row);
			}
		}
	}

	/*
	if(detected==true){
		for(int i=0; i<points.size(); ++i){
			cout 	<< points[i][0] << " "
						<< points[i][1] << " "
						<< points[i][2] <<endl;
		}
	}
	*/
	

	return detected;

}

int detectmpd(double x[3], double p[3], double offaxis){

	bool detected = false;

	double unitfactor = 1000; // default Pythia8: 1000

	double z0Argoncube = 574; // default DUNE: 574

	double l = (5 + z0Argoncube) * unitfactor; // largo del cilindro
	double r = 2.5 * unitfactor;
	double lmpd = 5 * unitfactor;


	double xsol[2];
	double ysol[2];
	double zsol[2];

	double ysol2a;
	double zsol2a;
	double ysol2b;
	double zsol2b;

	// LADO DEL CILINDRO

	// Discriminante
	double discr = -4*(
										pow(r,2)*(-1+pow(p[1],2)-pow(p[2],2))
										+2*r*p[1]*(l*p[1]+p[2]*x[1]-p[1]*x[2])
										+pow(l*p[1]+p[2]*x[1]-p[1]*x[2],2)
										);

	if (discr >= 0){

		xsol[0]	=	x[0]-(1/(pow(p[1],2)+pow(p[2],2)))*p[0]*
							(
								p[1]*x[1]-p[2]*(l+r-x[2])
								+sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);
		xsol[1]	=	x[0]+(1/(pow(p[1],2)+pow(p[2],2)))*p[0]*
							(
								-p[1]*x[1]+p[2]*(l+r-x[2])
								+sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);

		ysol[0] = (1/(pow(p[1],2)+pow(p[2],2)))*
							(
								l*p[1]*p[2] + p[2]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
								-p[1]*sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);
		ysol[1] = (1/(pow(p[1],2)+pow(p[2],2)))*
							(
								l*p[1]*p[2] + p[2]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
								+p[1]*sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);

		zsol[0]	=	x[2]-(1/(pow(p[1],2)+pow(p[2],2)))*p[2]*
							(
								p[1]*x[1]-p[2]*(l+r-x[2])
								+sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);
		zsol[1]	=	x[2]+(1/(pow(p[1],2)+pow(p[2],2)))*p[2]*
							(
								-p[1]*x[1]+p[2]*(l+r-x[2])
								+sqrt(
									-pow(l,2)*pow(p[1],2)+pow(r,2)*pow(p[2],2)-pow(p[2]*x[1]-p[1]*x[2],2)
									-2*l*p[1]*(r*p[1]+p[2]*x[1]-p[1]*x[2])
									+2*r*p[1]*(-p[2]*x[1]+p[1]*x[2])
									)
							);

		if ((-0.5*lmpd + offaxis*unitfactor <= xsol[0]) && (xsol[0] <= 0.5*lmpd + offaxis*unitfactor)){
			detected = true;
			//cout << "cilindro: "<< xsol[0] <<", "<< ysol[0] <<", " << zsol [0] << endl;
		}
		if ((-0.5*lmpd + offaxis*unitfactor <= xsol[1]) && (xsol[1] <= 0.5*lmpd + offaxis*unitfactor)){
			detected = true;
			//cout << "cilindro: "<< xsol[1] <<", " << ysol[1] <<", " << zsol [1] << endl;
		}	

		// Tapas del cilindro

		if (p[0]!=0){
			ysol2a = x[1] + p[1]*(0.5*lmpd + offaxis*unitfactor - x[0])/p[0];
			zsol2a = x[2] + p[2]*(0.5*lmpd + offaxis*unitfactor - x[0])/p[0];
			ysol2b = x[1] + p[1]*(-0.5*lmpd + offaxis*unitfactor - x[0])/p[0];
			zsol2b = x[2] + p[2]*(-0.5*lmpd + offaxis*unitfactor - x[0])/p[0];
			if (pow(ysol2a,2) + pow(zsol2a - l - r, 2) <= pow(r,2)){
				detected = true;
				//cout << "tapas: " << 0.5*lmpd + offaxis << ", " << ysol2a << ", " << zsol2a << endl;
			}
			if (pow(ysol2b,2) + pow(zsol2b - l - r, 2) <= pow(r,2)){
				detected = true;
				//cout << "tapas: " << -0.5*lmpd + offaxis << ", " << ysol2b << ", " << zsol2b << endl;
			}
		}
	}

	return detected;

}