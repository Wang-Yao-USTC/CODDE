#include <armadillo>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>
#include "rem.hpp"

using namespace arma;
using namespace std;

double delta (int m, int n) {
    if (m == n)
        return 1.0;
    else 
        return 0.0;
}

cx_mat InverseLiou(cx_mat& Q, cx_double gamma, const cx_mat& hams){

    const cx_double deom_ci = cx_double(0.0,1.0);
    
    int    n=Q.n_rows;
    cx_mat V=zeros<cx_mat>(n,n);
    cx_mat M=zeros<cx_mat>(n*n,n*n);
    cx_double T[n][n][n][n],temp;

        

    for(int i=0;i<n;++i){ 
    for(int j=0;j<n;++j){ 
    for(int k=0;k<n;++k){ 
    for(int l=0;l<n;++l){ 
        T[i][j][k][l]=deom_ci*hams(i,k)*delta(j,l)-deom_ci*hams(l,j)*delta(i,k)+gamma*delta(i,k)*delta(j,l);
    }
    }
    }
    }

    for(int i=0;i<n*n;++i){
    for(int j=0;j<n*n;++j){
        int itemp=i/n;int jtemp=i%n;
        int ktemp=j/n;int ltemp=j%n;
        M(i,j)=T[itemp][jtemp][ktemp][ltemp];
    }
    }

    M=M.i();

    for(int i=0;i<n;++i){ 
    for(int j=0;j<n;++j){ 
    for(int k=0;k<n;++k){ 
    for(int l=0;l<n;++l){ 
        T[i][j][k][l]=M(i*n+j,k*n+l);
    }
    }
    }
    }

    for (int i=0; i<n; ++i){
    for (int j=0; j<n; ++j){
        temp = cx_double(0.0,0.0);
        
        for (int p=0; p<n; ++p){
        for (int q=0; q<n; ++q){
            temp+=T[i][j][p][q]*Q(p,q);                           
        }
        } 
        V(i,j)=temp;            
    }
    } 

    return V;
}

double e(double t, const double ampl, const double freq, const double sigm){
    double et;
    const double deom_pi = 3.14159265358979;
    et = (ampl*deom_pi/(sqrt(2.0*deom_pi)*sigm))*exp(-0.5*t*t/(sigm*sigm))*cos(freq*t);
    return et;
}


void mcodde(const cx_cube& ddos, cx_cube& ddos1,const cx_cube& Vh,const double t, const double dt, const cx_mat& hams, const cx_cube& qmds, const mat& sdip, const cube& pdip,const vec& bdip, const ivec& mode,const cx_vec& etal, const cx_vec& etar, const vec& etaa, const cx_vec& expn, const vec& delr, const int swit, const double ampl, const double freq, const double sigm){

    const cx_double deom_ci = cx_double(0.0,1.0);
    const cx_double deom_c1 = cx_double(1.0,0.0);   

    int nsys = qmds.n_rows;
    int nmod = qmds.n_slices;
    int nper = etal.n_rows;
    int per= bdip.n_rows/nmod;

    cx_mat  hamt(hams); 
    
    // V
    cx_cube V= Vh;

    if (swit == 1) {
        // Ht
        hamt -= sdip*e(t,ampl,freq,sigm)*deom_c1;     
    }

    // =-i[H_t,rho]
    ddos1.slice(0)=-deom_ci*(hamt*ddos.slice(0)-ddos.slice(0)*hamt);
    
    // -R rho
    for (int m=0; m<nmod; ++m){
        cx_mat Q =zeros<cx_mat>(nsys,nsys);
        // \tilde Q
        for (int n=0; n<per; ++n){
            Q +=etal(m*per+n)*V.slice(m*per+n);    
        }
        ddos1.slice(0)-=qmds.slice(m)*(Q*ddos.slice(0)-ddos.slice(0)*Q.t())-(Q*ddos.slice(0)-ddos.slice(0)*Q.t())*qmds.slice(m);
    }
    

    //-i sum[Q,varrho]
    for (int iado=1; iado<=nper; ++iado){
        int nm=mode(iado-1);
        ddos1.slice(0)+=deom_ci*(qmds.slice(nm)*ddos.slice(iado)-ddos.slice(iado)*qmds.slice(nm));
        ddos1.slice(0)-=deom_ci*(qmds.slice(nm)*ddos.slice(iado+nper)-ddos.slice(iado+nper)*qmds.slice(nm));
    }

    for (int iado=1; iado<=nper; ++iado){
        // -i[Ht, varrho]-gamma*varrho
        ddos1.slice(iado)=-deom_ci*(hamt*ddos.slice(iado)-ddos.slice(iado)*hamt)-expn(iado-1)*ddos.slice(iado);
        ddos1.slice(iado+nper)=-deom_ci*(hamt*ddos.slice(iado+nper)-ddos.slice(iado+nper)*hamt)-expn(iado-1)*ddos.slice(iado+nper);
        
        if (swit == 1){
            ddos1.slice(iado)-=e(t,ampl,freq,sigm)*(etal(iado-1)*(sdip*V.slice(iado-1)-V.slice(iado-1)*sdip)*ddos.slice(0)-etar(iado-1)*ddos.slice(0)*(sdip*V.slice(iado-1)-V.slice(iado-1)*sdip));
            ddos1.slice(iado+nper)+=deom_ci*e(t,ampl,freq,sigm)*bdip(iado-1)*(etal(iado-1)-etar(iado-1))*ddos.slice(0);     
        }

    }
}

void mcodde(const cx_cube& ddos, cx_cube& ddos1,const cx_cube& Vh, const cx_mat& hams, const cx_cube& qmds, const mat& sdip, const cube& pdip,const vec& bdip, const ivec& mode,const cx_vec& etal, const cx_vec& etar, const vec& etaa, const cx_vec& expn, const vec& delr){

    const cx_double deom_ci = cx_double(0.0,1.0);   

    int nsys = qmds.n_rows;
    int nmod = qmds.n_slices;
    int nper = etal.n_rows;
    int per= bdip.n_rows/nmod;

    cx_mat  hamt(hams); 
    
    // V
    cx_cube V= Vh;

    // =-i[H_t,rho]
    ddos1.slice(0)=-deom_ci*(hamt*ddos.slice(0)-ddos.slice(0)*hamt);
    
    // -R rho
    for (int m=0; m<nmod; ++m){
        cx_mat Q =zeros<cx_mat>(nsys,nsys);
        // \tilde Q
        for (int n=0; n<per; ++n){
            Q +=etal(m*per+n)*V.slice(m*per+n);    
        }
        ddos1.slice(0)-=qmds.slice(m)*(Q*ddos.slice(0)-ddos.slice(0)*Q.t())-(Q*ddos.slice(0)-ddos.slice(0)*Q.t())*qmds.slice(m);
    }
    

    //-i sum[Q,varrho]
    for (int iado=1; iado<=nper; ++iado){
        int nm=mode(iado-1);
        ddos1.slice(0)+=deom_ci*(qmds.slice(nm)*ddos.slice(iado)-ddos.slice(iado)*qmds.slice(nm));
        ddos1.slice(0)-=deom_ci*(qmds.slice(nm)*ddos.slice(iado+nper)-ddos.slice(iado+nper)*qmds.slice(nm));
    }

    for (int iado=1; iado<=nper; ++iado){
        // -i[Ht, varrho]-gamma*varrho
        ddos1.slice(iado)=-deom_ci*(hamt*ddos.slice(iado)-ddos.slice(iado)*hamt)-expn(iado-1)*ddos.slice(iado);
        ddos1.slice(iado+nper)=-deom_ci*(hamt*ddos.slice(iado+nper)-ddos.slice(iado+nper)*hamt)-expn(iado-1)*ddos.slice(iado+nper);
    }
}
