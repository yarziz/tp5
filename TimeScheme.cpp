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
{
  _fin_vol->Build_flux_mat_and_rhs(_t);
  
  _sol+=-(_df->Get_dt())*( _fin_vol->Get_flux_matrix()*_sol+ _fin_vol->Get_BC_RHS())+(_df->Get_dt())*( _fin_vol->Source_term(_t));
  _t+=(_df->Get_dt());
}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
  _fin_vol->Build_flux_mat_and_rhs(_t);
  Eigen::SparseMatrix<double> ID;
  VectorXd b;
  ID.resize(_fin_vol->Get_flux_matrix().size(),_fin_vol->Get_flux_matrix().size());
  ID.setIdentity();
  _solver_direct.compute(ID+(_df->Get_dt())*_fin_vol->Get_flux_matrix());
  b=_sol-(_df->Get_dt())*(_fin_vol->Get_BC_RHS()-_fin_vol->Source_term(_t+(_df->Get_dt())));
  _sol= _solver_direct.solve(b);
}

#define _TIME_SCHEME_CPP
#endif
