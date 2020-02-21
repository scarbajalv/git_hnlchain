#!/usr/bin/env wolframscript
(* ::Package:: *)

#!/usr/bin/env wolframscript

Print["Majorana HNL production & decay"];

(* HNL DECAYS TOO MESONS/LEPTONS*)
(* ESTE SCRIPT GENERA LOS ARCHIVOS PCONFIG QUE UTILIZA PYTHIA8. *)
(* NOTAR QUE SIEMPRE TOMA COMO ARGUMENTO LA MASA DEL HNL *)
e=.5*10^-3; mu=0.1056583; tau=1.776; leptons={e,mu,tau};
tfactor=1.5193*10^24; (* Para pasar de s a GeV^-1 *)
mmc=3*10^11; (* para pasar de s a mm/c *)
ckm=({
 {0.9742, 0.2243, 0.0365},
 {0.218, 0.997, 0.0422},
 {0.0081, 0.0394, 1.019}
});
gf=1.166378*10^-5;
xw=0.2223; 
gll=-(1/2)+xw;
glr=xw;

(* tau meson *)
tau15 = {290.3*10^-15*tfactor, 1.77686, Null, Null, 15};

(* charged pseudoscalar *)
d431={504*10^-15*tfactor, 1.96834, 249.0*10^-3, ckm[[2,2]], 431, 0.0548};
d431bar={504*10^-15*tfactor, 1.96834,249.0*10^-3, ckm[[2,2]], -431, 0.0548}; 
d411={1040*10^-15*tfactor, 1.86958, 211.9*10^-3, ckm[[2,1]], 411, 0.0012};
d411bar={1040*10^-15*tfactor, 1.86958, 211.9*10^-3, ckm[[2,1]], -411, 0.0012};
k321={1.238*10^-8*tfactor,493.677*10^-3,155.6*10^-3,ckm[[1,2]]}; 
pi211={2.6*10^-8*tfactor,139.57061*10^-3,130.2*10^-3,ckm[[1,1]]}; 

(* Charged vector *)
k323={Null,895.55*10^-3,0.1827,ckm[[1,2]]}; (* K*(892) *)
rho770c={Null,775.26*10^-3,0.162,ckm[[1,1]]};

(* neutral pseudoscalar *)
pi111={Null,134.977*10^-3,130.2*10^-3,Null}; 
eta221={Null,547.862*10^-3,81.7*10^-3,Null}; (* eta *)
eta331={Null,957.78*10^-3,-94.7*10^-3,Null}; (* eta'(958) *)

(* Neutral vector *)
rho770n={Null,775.26*10^-3,0.162,1-2xw};
omega782={Null,782.65*10^-3,0.153,4/3 xw}; 
phi1020={Null, 1.019461,0.234,4/3 xw-1};

lambda[a_,b_,c_]:=a^2+b^2+c^2-2a*b-2b*c-2c*a;
i1[x_,y_]:=Sqrt[lambda[1,x,y]]((1-x)^2-y(1+x));
i2[x_,y_]:=Sqrt[lambda[1,x,y]]((1+x-y)(1+x+2y)-4x);

i11[x_,y_,z_]:=12NIntegrate[1/s (s-x-y)(1+z-s)Sqrt[lambda[s,x,y]] Sqrt[lambda[1,s,z]],{s,(Sqrt[x]+Sqrt[y])^2,(1-Sqrt[z])^2}]; 
i22[x_,y_,z_]:=24Sqrt[y*z]NIntegrate[1/s (1+x-s)Sqrt[lambda[s,y,z]] Sqrt[lambda[1,s,x]],{s,(Sqrt[y]+Sqrt[z])^2,(1-Sqrt[x])^2}];

datae=Import["./mainconfig/mixe.dat"];
datamu=Join[Table[{i,1},{i,0.01,0.1,0.01}],Import["./mainconfig/mixmu.dat"]];
datatau=Import["./mainconfig/mixtau.dat"];
fmixe=Interpolation[datae,InterpolationOrder->1];
fmixmu=Interpolation[datamu,InterpolationOrder->1];
fmixtau=Interpolation[datatau,InterpolationOrder->1];

width01[n_,l_,h_,u_]:=h[[4]]^2 u^2 (gf^2 h[[3]]^2 n^3)/(16Pi) i1[(l/n)^2,(h[[2]]/n)^2];
width02[n_,h_,u_]:=u^2*(gf^2 h[[3]]^2 n^3)/(16*Pi) i1[0,(h[[2]]/n)^2];
width03[n_,l_,h_,u_]:=(u^2*h[[4]]^2 gf^2 h[[3]]^2 n^3)/(16Pi*h[[2]]^2) i2[(l/n)^2,(h[[2]]/n)^2];
width04[n_,h_,u_]:=u^2*(gf^2 h[[3]]^2 h[[4]]^2 n^3)/(16Pi*h[[2]]^2) i2[0,(h[[2]]/n)^2];
width05[n_,l1_,l2_,u_]:=(gf^2 n^5)/(192Pi^3) u^2 i11[0,(l1/n)^2,(l2/n)^2];
width06a[n_,l_,u_]:=(gf^2 n^5)/(96Pi^3) u^2 ((gll*glr+glr)i22[0,(l/n)^2,(l/n)^2]+(gll^2+glr^2+1+2gll)i11[0,(l/n)^2,(l/n)^2]);
width06b[ n_,l_,u_]:=(gf^2 n^5)/(96Pi^3) u^2 ((gll*glr+0)i22[0,(l/n)^2,(l/n)^2]+(gll^2+glr^2+0)i11[0,(l/n)^2,(l/n)^2]);
(* nu_i nu_j nubar_j  Notar que se ha dividido entre 3 para tomar cada término de la sumatoria de Pascoli2009 Eq. C.8 *)
width07[n_,u_]:=u^2 (gf^2 n^5)/(3*96Pi^3); 

(* (1) : N > l- P+ *)
ePi[n_,u_]:=Re[width01[n,e,pi211,u]*UnitStep[n-pi211[[2]]-e]];
eK[n_,u_]:=Re[width01[n,e,k321,u]*UnitStep[n-k321[[2]]-e]];
eD[n_,u_]:=Re[width01[n,e,d411,u]*UnitStep[n-d411[[2]]-e]];
muPi[n_,u_]:=Re[width01[n,mu,pi211,u]*UnitStep[n-pi211[[2]]-mu]];
muK[n_,u_]:=Re[width01[n,mu,k321,u]*UnitStep[n-k321[[2]]-mu]];
muD[n_,u_]:=Re[width01[n,mu,d411,u]*UnitStep[n-d411[[2]]-mu]];
tauPi[n_,u_]:=Re[width01[n,tau,pi211,u]*UnitStep[n-pi211[[2]]-tau]];
tauK[n_,u_]:=Re[width01[n,tau,k321,u]*UnitStep[n-k321[[2]]-tau]];
tauD[n_,u_]:=Re[width01[n,tau,d411,u]*UnitStep[n-d411[[2]]-tau]];

(* (2) : N > l- V+ *)
eKstar[n_,u_]:=Re[width03[n,e,k323,u]*UnitStep[n-k323[[2]]-e]];
erhoc[n_,u_]:=Re[width03[n,e,rho770c,u]*UnitStep[n-rho770c[[2]]-e]];
muKstar[n_,u_]:=Re[width03[n,mu,k323,u]*UnitStep[n-k323[[2]]-mu]];
murhoc[n_,u_]:=Re[width03[n,mu,rho770c,u]*UnitStep[n-rho770c[[2]]-mu]];
tauKstar[n_,u_]:=Re[width03[n,tau,k323,u]*UnitStep[n-k323[[2]]-tau]];
taurhoc[n_,u_]:=Re[width03[n,tau,rho770c,u]*UnitStep[n-rho770c[[2]]-tau]];

(* (3) : N > nu P0 *)
nupi0[n_,u_]:=Re[width02[n,pi111,u]*UnitStep[n-pi111[[2]]]];
nueta[n_,u_]:=Re[width02[n,eta221,u]*UnitStep[n-eta221[[2]]]];
nuetaprime[n_,u_]:=Re[width02[n,eta331,u]*UnitStep[n-eta331[[2]]]];

(* (4) : N > nu V0 *)
nurho[n_,u_]:=Re[width04[n,rho770n,u]*UnitStep[n-rho770n[[2]]]];
nuomega[n_,u_]:=Re[width04[n,omega782,u]*UnitStep[n-omega782[[2]]]];
nuphi[n_,u_]:=Re[width04[n,phi1020,u]*UnitStep[n-phi1020[[2]]]];

(* (5) (par de leptones diferentes) N > l1- W+ (W+ > l2+ nu2) / l1!=l2  *)
enumumu[n_,u_]:=Re[width05[n,e,mu,u]*UnitStep[n-e-mu]];
munuee[n_,u_]:=Re[width05[n,mu,e,u]*UnitStep[n-mu-e]];
enutautau[n_,u_]:=Re[width05[n,e,tau,u]*UnitStep[n-e-tau]];
taunuee[n_,u_]:=Re[width05[n,tau,e,u]*UnitStep[n-tau-e]];
munutautau[n_,u_]:=Re[width05[n,mu,tau,u]*UnitStep[n-mu-tau]];
taunumumu[n_,u_]:=Re[width05[n,tau,mu,u]*UnitStep[n-tau-mu]];

(* (6a) (par de leptones y neutrino iguales) N > nu1 l1- l1+ *)
nueee[n_,u_]:=Re[width06a[n,e,u]*UnitStep[n-2e]];
numumumu[n_,u_]:=Re[width06a[n,mu,u]*UnitStep[n-2mu]];
nutautautau[n_,u_]:=Re[width06a[n,tau,u]*UnitStep[n-2tau]];

(* (6b) (par de leptones iguales, neutrino diferente) N > nu1 l2- l2+ *)
nuemumu[n_,u_]:=Re[width06b[n,mu,u]*UnitStep[n-2mu]];
nuetautau[n_,u_]:=Re[width06b[n,tau,u]*UnitStep[n-2tau]];
numuee[n_,u_]:=Re[width06b[n,e,u]*UnitStep[n-2e]];
numutautau[n_,u_]:=Re[width06b[n,tau,u]*UnitStep[n-2tau]];
nutauee[n_,u_]:=Re[width06b[n,e,u]*UnitStep[n-2e]];
nutaumumu[n_,u_]:=Re[width06b[n,mu,u]*UnitStep[n-2mu]];

(* (7) : nui nuj nujbar *)
nununu[n_,u_]:=Re[width07[n,u]];

(* Notar que los factores "2" solo estan presentes para Majorana. Para Dirac solo cambiarlos a "1" *)
totalw[n_,mixe_,mixmu_,mixtau_]:=
(* 65 canales *)
2*(ePi[n,mixe]+muPi[n,mixmu]+tauPi[n,mixtau]+eK[n,mixe]+muK[n,mixmu]+eD[n,mixe])+ (*1*)
2*(eKstar[n,mixe]+muKstar[n,mixmu]+erhoc[n,mixe]+murhoc[n,mixmu])+ (*2*)
1*(nupi0[n,mixe]+nupi0[n,mixmu]+nupi0[n,mixtau]+nueta[n,mixe]+nueta[n,mixmu]+nueta[n,mixtau]+nuetaprime[n,mixe]+nuetaprime[n,mixmu]+nuetaprime[n,mixtau])+ (*3*)
1*(nurho[n,mixe]+nurho[n,mixmu]+nurho[n,mixtau]+nuomega[n,mixe]+nuomega[n,mixmu]+nuomega[n,mixtau]+nuphi[n,mixe]+nuphi[n,mixmu]+nuphi[n,mixtau])+ (*4*)
2*(enumumu[n,mixe]+enutautau[n,mixe]+munuee[n,mixmu]+munutautau[n,mixmu]+taunuee[n,mixtau]+taunumumu[n,mixtau])+ (*5*)
1*(nueee[n,mixe]+numumumu[n,mixmu])+ (*6a*)
1*(numuee[n,mixmu]+nutauee[n,mixtau]+nuemumu[n,mixe]+nutaumumu[n,mixtau])+ (*6b*)
3*(nununu[n,mixe]+nununu[n,mixmu]+nununu[n,mixtau]); (*7*)


(* CC: Dirac & Majorana / 16 canales *)
expbrcc1[n_,mixe_,mixmu_,mixtau_]:={
	(* (1) - 6 canales *)
	"2000000001:oneChannel = 1 "<>ToString[AccountingForm[eD[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 411",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eK[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muK[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[ePi[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muPi[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[tauPi[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 211",	
	(* (2) - 4 canales*)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eKstar[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muKstar[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[erhoc[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 213",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[murhoc[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 213",
	(* (5) - 6 canales*)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[enumumu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 14 -13",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[enutautau[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 11 16 -15",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munuee[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 12 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munutautau[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 13 16 -15",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunuee[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 12 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 15 14 -13"};

	(* NC: Majorana % Dirac / 33 canales *)
	expbrnc[n_,mixe_,mixmu_,mixtau_]:={
	(* (3) - 9 canales *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 111",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 111",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nupi0[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 111",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 221",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 221",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueta[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 221",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 331",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 331",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuetaprime[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 331",
	(* (4) - 9 canales *)	
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 113",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 113",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nurho[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 113",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 223",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 223",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuomega[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 223",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 333",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 333",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuphi[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 333",
	(* (6a) - 2 canales [¿memode  = 22/0?] *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nueee[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 11 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[numumumu[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 13 -13",
	(* (6b) - 4 canales [¿memode  = 22/0?] *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[numuee[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 11 -11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nutauee[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 11 -11",	
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nuemumu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 13 -13",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nutaumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 13 -13",
	(* (7) - 9 canales *)
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 12 -12",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 14 -14",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 12 16 -16",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 12 -12",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 14 -14",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 14 16 -16",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 12 -12",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 14 -14",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[nununu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 16 16 -16"};

(* CC: Solo Majorana / 16 canales *)
expbrcc2[n_,mixe_,mixmu_,mixtau_]:={
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[ePi[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muPi[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[tauPi[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -15 -211",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eK[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muK[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -321",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eD[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -411",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[eKstar[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[muKstar[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -323",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[erhoc[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -213",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[murhoc[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -213",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[enumumu[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -14 13",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[enutautau[n,mixe]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -11 -16 15",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munuee[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -12 11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[munutautau[n,mixmu]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -13 -16 15",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunuee[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -15 -12 11",
	"2000000001:addChannel = 1 "<>ToString[AccountingForm[taunumumu[n,mixtau]/totalw[n,mixe,mixmu,mixtau],16]]<>" 0 -15 -14 13"};




(* Final string with all the HNL Brs *)
(* exbrcc2 solo se agrega para Majorana (asegurarse también de haber cambiado los factores en totalw) *)
expbr[n_,mixe_,mixmu_,mixtau_]:=Join[expbrcc1[n,mixe,mixmu,mixtau],expbrnc[n,mixe,mixmu,mixtau],expbrcc2[n,mixe,mixmu,mixtau]];



(*************** MESON DECAY ****************** *)
(* meson decay BRs function *)
width2body[meson_, l_, n_, u_]:=
	Re[(gf^2*meson[[3]]^2*meson[[2]]*n^2)/(8*Pi)*meson[[4]]^2*u^2
	*(1-n^2/meson[[2]]^2+2l^2/meson[[2]]^2+l^2/n^2(1-l^2/meson[[2]]^2))	
	*Sqrt[(1 + n^2/meson[[2]]^2 - l^2/meson[[2]]^2)^2-4*n^2/meson[[2]]^2]
	*UnitStep[meson[[2]] - l]];
newwidth[meson_, n_, mixe_, mixmu_, mixtau_] :=  
	Re[1/meson[[1]] 
  + width2body[meson, e, n, mixe] 
  + width2body[meson, mu, n, mixmu] 
  + width2body[meson, tau, n, mixtau]];
newwidth2hnl[meson_, n_, mixe_, mixmu_, mixtau_] :=  
	Re[width2body[meson, e, n, mixe] 
  + width2body[meson, mu, n, mixmu] 
  + width2body[meson, tau, n, mixtau]];

mesonparam[meson_,n_,mixe_,mixmu_,mixtau_]:={
	ToString[meson[[5]]]<>":mayDecay = true",
	ToString[meson[[5]]]<>":tau0 = "<>ToString[AccountingForm[Re[meson[[1]]*1/tfactor*mmc],16]],
	ToString[meson[[5]]]<>":m0 = "<>ToString[AccountingForm[meson[[2]],16]]};
(* Meson decay BRs*)
mesondecay[meson_,n_,mixe_,mixmu_,mixtau_]:={
	ToString[meson[[5]]]<>":oneChannel = 1 1 0 -15 16"};
mesondecaybar[meson_,n_,mixe_,mixmu_,mixtau_]:={
	ToString[meson[[5]]]<>":oneChannel = 1 1 0 15 -16"};


(******************************************************************)
(************************* TAU DECAY ******************************)
(* tau > HNL P  - Shaposhnikov *)
(* u siempre es utau *)
tauwidth01[meson_,n_,u_]:=
Re[u^2/(16Pi) gf^2*meson[[4]]^2*meson[[3]]^2*tau15[[2]]^2*((1-n^2/tau15[[2]]^2)^2
	-meson[[2]]^2/tau15[[2]]^2 (1+n^2/tau15[[2]]^2))*Sqrt[(1-(meson[[2]]^2-n^2)/tau15[[2]]^2)
	(1-(meson[[2]]^2+n^2)/tau15[[2]]^2)]*UnitStep[tau15[[2]]-n-meson[[2]]]];
(* tau > HNL V - Shaposhnikov *)
(* u siempre es utau *)
tauwidth02[meson_,n_,u_]:=
Re[u^2/(16Pi) gf^2*meson[[4]]^2*meson[[3]]^2*tau15[[2]]^2*((1-n^2/tau15[[2]]^2)^2
	-meson[[2]]^2/tau15[[2]]^2 (1+n^2/tau15[[2]]^2-2*meson[[2]]^2/tau15[[2]]^2))
	*Sqrt[(1-(meson[[2]]^2-n^2)/tau15[[2]]^2)(1-(meson[[2]]^2+n^2)/tau15[[2]]^2)]*UnitStep[tau15[[2]]-n-meson[[2]]]];
(* tau > HNL l1 nubar1 - Bondarenko 2.8 *)
(* u SIEMPRE es utau *)
tauwidth03[n_,l_,u_]:=
Re[(gf^2 tau15[[2]]^5)/(96Pi^3) u^2*NIntegrate[1/z^3 (z-(l/tau15[[2]])^2)^2 
*Sqrt[lambda[1,z,(n/tau15[[2]])^2]]*((z+2(l/tau15[[2]])^2)(1-(n/tau15[[2]])^2)^2+z(z-(l/tau15[[2]])^2)(1+(n/tau15[[2]])^2
	-(l/tau15[[2]])^2)-z (l/tau15[[2]])^4-2z^3),{z,(l/tau15[[2]])^2,(1-(n/tau15[[2]]))^2}]*UnitStep[tau15[[2]]-n-l]];
(* tau > nutau l1 HNL - Bondarenko 2.8 *)
(* u puede ser de cualquier sabor *)
tauwidth04[l_,n_,u_]:=
Re[(gf^2 tau15[[2]]^5)/(96Pi^3) u^2*NIntegrate[1/z^3 (1-z)^2 Sqrt[lambda[z,(n/tau15[[2]])^2,(l/tau15[[2]])^2]]*(2z^3+z-z(1-z)
	(1-(n/tau15[[2]])^2-(l/tau15[[2]])^2)-(2+z)((n/tau15[[2]])^2-(l/tau15[[2]])^2)^2)
,{z,(l/tau15[[2]]+n/tau15[[2]])^2,1}]*UnitStep[tau15[[2]]-l-n]];

(* u siempre es utau *)
tau2PiN[n_,mixtau_]:=tauwidth01[pi211,n,mixtau];
tau2KN[n_,mixtau_]:=tauwidth01[k321,n,mixtau];
tau2KstarN[n_,mixtau_]:=tauwidth02[k323,n,mixtau];
tau2rhocN[n_,mixtau_]:=tauwidth02[rho770c,n,mixtau];
tau2Nenue[n_,mixtau_]:=tauwidth03[n,e,mixtau];
tau2Nmunumu[n_,mixtau_]:=tauwidth03[n,mu,mixtau];
(* u puede ser cualquiera *)
tau2nutaueN[n_,mixe_]:=tauwidth04[e,n,mixe];
tau2nutaumuN[n_,mixmu_]:=tauwidth04[mu,n,mixmu];
tautotalw[n_,mixe_,mixmu_,mixtau_]:=
Re[1/tau15[[1]]
+tau2PiN[n,mixtau]+tau2KN[n,mixtau]
+tau2KstarN[n,mixtau]+tau2rhocN[n,mixtau]
+tau2Nenue[n,mixtau]+tau2Nmunumu[n,mixtau]
+tau2nutaueN[n,mixe]+tau2nutaumuN[n,mixmu]];
tau2hnlw[n_,mixe_,mixmu_,mixtau_]:=
Re[tau2PiN[n,mixtau]+tau2KN[n,mixtau]
+tau2KstarN[n,mixtau]+tau2rhocN[n,mixtau]
+tau2Nenue[n,mixtau]+tau2Nmunumu[n,mixtau]
+tau2nutaueN[n,mixe]+tau2nutaumuN[n,mixmu]];

(* tau param *)
tauparam[n_,mixe_,mixmu_,mixtau_]:=
{"15:mayDecay = true",
"15:tau0 = "<>ToString[AccountingForm[Re[1/tautotalw[n,mixe,mixmu,mixtau]]*1/tfactor*mmc,16]],
"15:m0 = "<>ToString[AccountingForm[tau15[[2]],16]]};

(* BRs a HNL. Escalados, suman 1. *)
expbrtau[n_,mixe_,mixmu_,mixtau_]:={
(*(1)-6 canales*)
"15:oneChannel = 1 "<>ToString[AccountingForm[tau2PiN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 2000000001 -211",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2KN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 2000000001 -321",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2KstarN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 2000000001 -323",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2rhocN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 2000000001 -213",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2Nenue[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 2000000001 11 -12",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2Nmunumu[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 2000000001 13 -14",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2nutaueN[n,mixe]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 16 11 2000000001",
"15:addChannel = 1 "<>ToString[AccountingForm[tau2nutaumuN[n,mixmu]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 16 13 2000000001"};

(* taubar BRs: NO ES NECESARIO *)
(*
expbrtaubar[n_,mixe_,mixmu_,mixtau_]:={
(*(1)-6 canales*)
"-15:oneChannel = 1 "<>ToString[AccountingForm[tau2PiN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 -2000000001 211",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2KN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 -2000000001 321",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2KstarN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 -2000000001 323",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2rhocN[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1521 -2000000001 213",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2Nenue[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 -2000000001 -11 12",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2Nmunumu[n,mixtau]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 -2000000001 -13 14",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2nutaueN[n,mixe]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 -16 -11 2000000001",
"-15:addChannel = 1 "<>ToString[AccountingForm[tau2nutaumuN[n,mixmu]/tau2hnlw[n,mixe,mixmu,mixtau],16]]<>" 1531 -16 -13 2000000001"};
*)

(* FACTOR BR *)
expbrtaufactor[meson_,n_,mixe_,mixmu_,mixtau_]:={
	"# BR factor "<>ToString[meson[[5]]]<>": "<>ToString[AccountingForm[(meson[[6]])*(meson[[1]]^(-1)/(meson[[1]]^(-1)
		+newwidth2hnl[meson,n,mixe,mixmu,mixtau]))*(tau2hnlw[n,mixe,mixmu,mixtau]/tautotalw[n,mixe,mixmu,mixtau]),16]]}

(*  HNL CONFIG  *)
hnlparams = {"2000000001:new",
"2000000001:isResonance = false",
"2000000001:name = HNL",
"#2000000001:antiName = HNLbar", (* Esto lo define como Majorana *)
"2000000001:spinType = 2",
"2000000001:chargeType = 0",
"2000000001:colType = 0"}

hnlmass[n_]:={"2000000001:m0 = "<>ToString[AccountingForm[n,16]]}
(* lifetime  *)
lifetime[n_,mixe_,mixmu_,mixtau_]:=Quiet[Re[1/totalw[n,mixe,mixmu,mixtau]]*1/tfactor*mmc];
(* PYTHIA8 string *)
explifetime[n_,mixe_,mixmu_,mixtau_]:=Quiet[{
	"2000000001:tau0 = "<>ToString[AccountingForm[Re[1/totalw[n,mixe,mixmu,mixtau]]*1/tfactor*mmc,16]]}]; 

(* Join ALL config for HNL (particleData + BRs) *)
hnldecay[n_,mixe_,mixmu_,mixtau_]:=Join[hnlparams,hnlmass[n],explifetime[n,mixe,mixmu,mixtau],expbr[n,mixe,mixmu,mixtau]];

(* EXPORTA PYTHIA CONFIG FILES PARA "nfiles" MASAS CON <0,1.9>.  *)

Print["Todos los par\[AAcute]metros y funciones definidos correctamente."]

paso=ToExpression[$ScriptCommandLine[[2]]];
names={};
Do[
name=ToString[n]<>".dat";
Print["Generando "<>name<>" ..."];
AppendTo[names,name];
Quiet[Export["./mathbrs/"<>name,
	Join[
		{"#"<>ToString[n]},
		expbrtaufactor[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		expbrtaufactor[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		hnldecay[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesonparam[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesonparam[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesondecay[d431,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesondecay[d411,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesondecaybar[d431bar,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		mesondecaybar[d411bar,n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		tauparam[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]],
		expbrtau[n,Sqrt[fmixe[n]],Sqrt[fmixmu[n]],Sqrt[fmixtau[n]]]
		]]];
,{n,paso,2-paso,paso}];
AppendTo[names,""];
Export["./mathbrs/confignames.dat",names];
Print["All done"]
