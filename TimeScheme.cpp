#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_df(data_file), _fin_vol(adv), _t(_df->Get_t0()), _sol(adv->Initial_condition())
{
}

EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
  std::cout << "Build time scheme class." << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Renvoie _sol (pratique pour vérifier la résolution)
const VectorXd & TimeScheme::Get_sol() const
{
  return _sol;
}

// Euler Explicite
void EulerScheme::Advance()
{//todo
  _fin_vol->Build_flux_mat_and_rhs(_t);//attention au _t
  SparseMatrix<double> A=_fin_vol->Get_flux_matrix();
  VectorXd S=_fin_vol->Source_term(_t);
  double dt=_df->Get_dt();
  VectorXd b=_fin_vol->Get_BC_RHS();
  //VectorXd phi=_fin_vol->Initial_condition();

  _sol=_sol-dt*(A*_sol+b)+dt*S;
  //_t=_t+dt;
}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{//todo
  _fin_vol->Build_flux_mat_and_rhs(_t);//attention au _t
  SparseMatrix<double> A=_fin_vol->Get_flux_matrix();
  double dt=_df->Get_dt();
  VectorXd S=_fin_vol->Source_term(_t+dt);
  SparseMatrix<double> I(_sol.size(),_sol.size());
  I.setIdentity();
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  VectorXd b=_fin_vol->Get_BC_RHS();
  //double dt=_df->Get_dt();
  solver.analyzePattern(dt*A+I);

  solver.factorize(dt*A+I);

  _sol= solver.solve(_sol-dt*b+dt*S);
}

#define _TIME_SCHEME_CPP
#endif
