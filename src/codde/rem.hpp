#ifndef REM_H_
#define REM_H_  

#include<armadillo>

using namespace arma;

double delta (int m, int n);

cx_mat InverseLiou(cx_mat& Q, cx_double gamma, const cx_mat& hams);

double e(double t, const double ampl, const double freq, const double sigm);

void mcodde(const cx_cube& ddos, cx_cube& ddos1,const cx_cube& Vh,const double t, const double dt, const cx_mat& hams, const cx_cube& qmds, const mat& sdip, const cube& pdip,const vec&  bdip, const ivec& mode,const cx_vec& etal, const cx_vec& etar, const vec& etaa, const cx_vec& expn, const vec& delr, const int swit, const double ampl, const double freq, const double sigm);

void mcodde(const cx_cube& ddos, cx_cube& ddos1,const cx_cube& Vh, const cx_mat& hams, const cx_cube& qmds, const mat& sdip, const cube& pdip,const vec& bdip, const ivec& mode,const cx_vec& etal, const cx_vec& etar, const vec& etaa, const cx_vec& expn, const vec& delr);



#endif
