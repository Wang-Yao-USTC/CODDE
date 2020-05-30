#include <iostream>
#include <armadillo>
#include <complex>
#include <cstdio>

#include "json11.hpp"
#include "rem.hpp"

using namespace std;
using namespace arma;
using namespace json11;


static double entropy (const cx_mat& rho) {
	double ferr=0.00000001;
    vec eval = eig_sym(rho);
    int nsys=rho.n_rows;
    if (any(eval<0)) {
        return -9527;
    } else {
	double    entrop=0.0;
	for(int i=0;i<nsys;++i){
		if(abs(eval(i))>ferr){
			entrop -=eval(i)*log(eval(i));
		}
	}
        return entrop;
    }
}

cx_double trace(const cx_mat& rho){

    int n=rho.n_rows;
    cx_double tracevalue=cx_double(0.0,0.0);
    for (int i=0; i<n; ++i){
        tracevalue +=rho(i,i);
    }
    return tracevalue;
}

int main(){

    
    // *** const *** //
    const cx_double deom_ci = cx_double(0.0,1.0);
    const cx_double deom_c1 = cx_double(1.0,0.0);
    const double deom_pi = 3.14159265358979; 

    // *** initial *** //

    cx_mat hams;
    if (hams.load ("inp_hams.mat", arma_ascii)) {
        hams.print("inp_hams.mat");
    } else {
        printf ("Fail to load hams !\n");
    }
     
    cx_cube qmds;  
    if (qmds.load ("inp_qmds.mat", arma_ascii)) {
        qmds.print("inp_qmds.mat");
    } else {
        printf ("Fail to load qmds ! \n");
    }

    mat sdip;
    if (sdip.load ("inp_sdip.mat", arma_ascii)) {
        sdip.print("inp_sdip.mat");
    } else {
        printf ("Fail to load sdip ! \n");
    }
    
    cube pdip;
    if (pdip.load ("inp_pdip.mat", arma_ascii)) {
        pdip.print("inp_pdip.mat");
    } else {
        printf ("Fail to load pdip ! \n");
    }

    vec bdip;
    if (bdip.load ("inp_bdip.mat", arma_ascii)) {
        bdip.print("inp_bdip.mat");
    } else {
        printf ("Fail to load bdip ! \n");
    }

    cx_vec etal;
    if (etal.load ("inp_etal.mat", arma_ascii)) {
        etal.print("inp_etal.mat");
    } else {
        printf ("Fail to load etal ! \n");
    }

    cx_vec etar;
    if (etar.load ("inp_etar.mat", arma_ascii)) {
        etar.print("inp_etar.mat");
    } else {
        printf ("Fail to load etar ! \n");
    }

    vec etaa;
    if (etaa.load ("inp_etaa.mat", arma_ascii)) {
        etaa.print("inp_etaa.mat");
    } else {
        printf ("Fail to load etaa ! \n");
    }

    vec delr; // useless temporarily
    if (delr.load ("inp_delr.mat", arma_ascii)) {
        delr.print("inp_delr.mat");
    } else {
        printf ("Fail to load delr ! \n");
    }

    cx_vec expn;
    if (expn.load ("inp_expn.mat", arma_ascii)) {
        expn.print("inp_expn.mat");
    } else {
        printf ("Fail to load expn ! \n");
    }

    ivec mode;
    if (mode.load ("inp_mode.mat", arma_ascii)) {
        mode.print("inp_mode.mat");
    } else {
        printf ("Fail to load mode ! \n");
    }

    int nsys,nmod,nper,nmax;
    nsys = qmds.n_rows;
    nmod = qmds.n_slices;
    nper = etal.n_rows;
    nmax = 2*nper+1;

    cout << "nsys:"<< nsys << endl;
    cout << "nmod:"<< nmod << endl;
    cout << "nper:"<< nper << endl;
    cout << "nmax:"<< nmax << endl;

    // *** preprop *** //

    ifstream jsonFile("input.json");
    stringstream strStream;
    strStream << jsonFile.rdbuf();
    string jsonStr = strStream.str();
    string err;

    const Json json = Json::parse(jsonStr,err);
    if (!err.empty()) {
        printf ("Error in parsing input file: %s\n", err.c_str());
        return 0;
    }

    const int    inistate = json["td-rhot"]["inistate"].int_value();
	const double dt = json["td-rhot"]["dt"].number_value();
	const double ti = json["td-rhot"]["ti"].number_value();
	const double tf = json["td-rhot"]["tf"].number_value();
    const int    swit = json["td-rhot"]["pulse"]["swit"].int_value();
    const double ampl = json["td-rhot"]["pulse"]["ampl"].number_value();
    const double freq = json["td-rhot"]["pulse"]["freq"].number_value();
    const double sigm = json["td-rhot"]["pulse"]["sigm"].number_value();
 
    cout << "inistate:"<< inistate << endl;
	cout << "dt:"<< dt << endl;
	cout << "ti:"<< ti << endl;
	cout << "tf:"<< tf << endl;
    cout << "swit:"<< swit << endl;     
    cout << "ampl:"<< ampl << endl;
    cout << "freq:"<< freq << endl;
    cout << "sigm:"<< sigm << endl;

    // *** propagation *** //

    int nt_i = ceil(abs(ti)/dt);
    int nt_f = ceil(abs(tf)/dt);
    int nt = nt_i+nt_f;

    double t = -(nt_i-1)*dt;

    cx_cube ddos = zeros<cx_cube>(nsys,nsys,nmax);
    cx_cube K1 = zeros<cx_cube>(nsys,nsys,nmax);
    cx_cube K2 = zeros<cx_cube>(nsys,nsys,nmax);
    cx_cube K3 = zeros<cx_cube>(nsys,nsys,nmax);
    cx_cube K4 = zeros<cx_cube>(nsys,nsys,nmax);

    ddos(inistate,inistate,0) = 1.0;

    cx_cube Vh=zeros<cx_cube>(nsys,nsys,nper);

    for (int n=0; n<nper; ++n){
        int nm=mode(n);
        cx_mat qmdstemp=qmds.slice(nm);
        Vh.slice(n)=InverseLiou(qmdstemp,expn(n),hams);
    }
    
    FILE *frho = fopen("td-rhot.dat","w");
    FILE *fent = fopen("entropy.dat","w");
    FILE *fpol = fopen("polar.dat","w");
    cx_vec pt=zeros<cx_vec>(nt);

    for (int it=0; it<nt; ++it) {
        
        printf ("%s: %5.1f%% \n", "td-rhot",100*it/static_cast<double>(nt));

        mcodde(ddos,K1,Vh,t,dt,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr,swit,ampl,freq,sigm);
        mcodde(ddos+K1*dt/2.0,K2,Vh,t+dt/2.0,dt,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr,swit,ampl,freq,sigm);
        mcodde(ddos+K2*dt/2.0,K3,Vh,t+dt/2.0,dt,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr,swit,ampl,freq,sigm);
        mcodde(ddos+K3*dt,K4,Vh,t+dt,dt,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr,swit,ampl,freq,sigm);
        ddos=ddos+(K1+2.0*K2+2.0*K3+K4)*dt/6.0;
        t += dt;

        fprintf(frho,"%16.6e",t);
        for (int i=0; i<nsys; ++i) {
            fprintf(frho,"%16.6e",real(ddos(i,i,0)));
        }
        fprintf(frho,"\n");

        cx_mat rho=ddos.slice(0);
        double tmp = entropy(rho);
        if (tmp < 0) {
            fprintf (fent, "%16.6e nan\n", t);
        } else {
            fprintf (fent, "%16.6e%16.6e", t, tmp);
        }
        fprintf (fent, "\n");

        fprintf(fpol,"%16.6e",t);
        cx_double tmp1=cx_double(0.0,0.0);
        cx_mat sdipe;
        sdipe = sdip*deom_c1;
        for (int n=0; n<nper;++n){
            sdipe-= deom_ci*(etal(n)-etar(n))*bdip(n)*Vh.slice(n);
        }
        tmp1+=trace(sdipe*ddos.slice(0));
        for (int iado=1; iado<=nper;++iado){
            tmp1-=bdip(iado-1)*trace(ddos.slice(iado));
            tmp1+=bdip(iado-1)*trace(ddos.slice(iado+nper));
        }
        fprintf(fpol,"%16.6e%16.6e",real(tmp1),imag(tmp1));
        fprintf(fpol,"\n");
        pt(it)=tmp1;
    }

    fclose(frho);
    fclose(fent);
    fclose(fpol);

    // 1D FFT 

    pt-=pt(nt-1);

    FILE *fpoliw = fopen("polar-iw.dat","w");
    cx_vec ftm=zeros<cx_vec>(nt_i);
    cx_vec ftp=zeros<cx_vec>(nt_f);
    for (int i=0;i<nt_i;++i){
       ftm(i) = pt(nt_i-1-i,0);  
    }
    for (int i=0;i<nt_f;++i){
       ftp(i) = pt(nt_i+i,0);  
    }
    ftm.row(0) *= 0.5;
    ftp.row(0) *= 0.5;
    const cx_vec& fw = ifft(deom_c1*real(ftp))*nt_f*dt+fft(deom_c1*real(ftm))*dt;
    
    // write freq-domain signal
    const double dw1 = 2*deom_pi/nt_f/dt;
    const vec& fw_w1 = linspace(0.0,dw1*(nt_f-1),nt_f);
    const mat& fw_re = real(fw.rows(0,nt_f-1));
    const mat& fw_im = imag(fw.rows(0,nt_f-1));
    for (int iw=0; iw<nt_f/2; ++iw) {
        fprintf (fpoliw,"%16.6e%16.6e\n",fw_w1(iw),fw_im(iw));
    }
    fclose(fpoliw);
    return 0;
}
