#include <iostream>
#include <armadillo>
#include <complex>
#include <cstdio>
#include <cstdlib>

#include "json11.hpp"
#include "rem.hpp"

using namespace std;
using namespace arma;
using namespace json11;

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

    cout << "------------------------------   para   ------------------------------"<< endl;
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

    const int    inistate = json["corr"]["inistate"].int_value();
	const double dt = json["corr"]["dt"].number_value();
	const int nt = json["corr"]["nt"].int_value();
	const int n4eq = json["corr"]["n4eq"].int_value();

 
    cout << "inistate:"<< inistate << endl;
	cout << "dt:"<< dt << endl;
	cout << "nt:"<< nt << endl;
	cout << "n4eq:"<< n4eq << endl;

    
    // *** propagation *** //

    // int per=nper/nmod;


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
    
    FILE *fres = fopen("resp.dat","w");
    cx_vec pt=zeros<cx_vec>(nt);

    for (int it=0; it<n4eq; ++it) {
        
        printf ("%s: %5.1f%% \n", "td-rhot",100*it/static_cast<double>(nt));


        mcodde(ddos,K1,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K1*dt/2.0,K2,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K2*dt/2.0,K3,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K3*dt,K4,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        ddos=ddos+(K1+2.0*K2+2.0*K3+K4)*dt/6.0;
    }

    
    cx_cube ddos_tmp=ddos;

    ddos.slice(0)=sdip*ddos_tmp.slice(0);

    for (int iado=1; iado<=nper;++iado){
        ddos.slice(iado)=deom_ci*etal(iado-1)*(sdip*Vh.slice(iado-1)-Vh.slice(iado-1)*sdip)*ddos_tmp.slice(0);
        ddos.slice(iado)+=sdip*ddos_tmp.slice(iado);
        ddos.slice(iado+nper)=bdip(iado-1)*etal(iado-1)*ddos_tmp.slice(0);
        ddos.slice(iado+nper)+=sdip*ddos_tmp.slice(iado+nper);
    }
    
  
    double t=0.0;
    for (int it=0; it<nt; ++it) {
        
        printf ("%s: %5.1f%% \n", "td-rhot",100*it/static_cast<double>(nt));

        fprintf(fres,"%16.6e",t);
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

        //tmp1=tmp1-conj(tmp1);

        fprintf(fres, "%16.6e%16.6e",real(tmp1),imag(tmp1));
        fprintf(fres,"\n");
        pt(it)=tmp1;

        mcodde(ddos,K1,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K1*dt/2.0,K2,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K2*dt/2.0,K3,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);
        mcodde(ddos+K3*dt,K4,Vh,hams,qmds,sdip,pdip,bdip,mode,etal,etar,etaa,expn,delr);

        ddos=ddos+(K1+2.0*K2+2.0*K3+K4)*dt/6.0;
        t += dt;
   
    }
    fclose(fres);

    // 1D FFT
    //cx_vec ft = shift(pt,-nt_i);
    cx_vec ft = pt;
    ft(0) *= 0.5;
    const vec& w = linspace(0,2*deom_pi/dt,nt);
    const cx_vec& fw = ifft(ft)*nt*dt;
    vec fwr=real(fw);
    vec fwi=imag(fw);
    w.save("w.dat",raw_ascii);
	fwr.save("freq-corr-real.dat",raw_ascii); 
    fwi.save("freq-corr-imag.dat",raw_ascii);   
     
    return 0;
}
