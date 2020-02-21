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