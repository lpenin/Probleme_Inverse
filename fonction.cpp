#include "fonction.h"

double dot(std::vector<double> uloc, std::vector<double> vloc)
{
  // //Calcul le produit scalaire entre 2 vecteurs et envoie le résultat à tous les processeurs.


  int Np = uloc.size();
  double y,cloc;

  //Calcul d'une partie de y
  cloc = 0;
  for (int i=0; i<Np; i++)
  {
    cloc += uloc[i]*vloc[i];
  }

  return y;
}

std::vector<double> prodMatVec (std::vector<std::vector<double>> Aloc, std::vector<double> bloc)
{
  int nx= bloc.size();
  int ny= Aloc.size()/nx;
  std::vector<double> cloc(ny);
  for (int i=0; i<ny; i++)
  {
    cloc[i]=dot(Aloc[i],bloc);
  }
  return cloc;
}


std::vector<double> prodM5V (std::vector<std::vector<double>> A, std::vector<double> bloc, int nx, int ny)
{
  int np= nx*ny;
  std::vector<double> cloc(np);

  cloc[0] = A[2][0]*bloc[0] + A[3][0]*bloc[1] + A[4][0]*bloc[nx];
  cloc[np-1] = A[2][np-1]*bloc[np-1] + A[3][np-1]*bloc[np-1] + A[4][np-1]*bloc[nx];

  for (int i= 1; i<nx; i++)
  {
    cloc[i] = A[1][i]*bloc[i-1] + A[2][i]*bloc[i] + A[3][i]*bloc[i+1] + A[4][i]*bloc[i+nx];
  }
  for (int i= nx; i<np-nx; i++)
  {
    cloc[i] = A[0][i]*bloc[i-nx] + A[1][i]*bloc[i-1] + A[2][i]*bloc[i] + A[3][i]*bloc[i+1] + A[4][i]*bloc[i+nx];
  }
  for (int i= np-nx; i<np-1; i++)
  {
    cloc[i] = A[0][i]*bloc[i-nx] + A[1][i]*bloc[i-1] + A[2][i]*bloc[i] + A[3][i]*bloc[i+1];
  }

  return cloc;
}



std::vector<double> CG (std::vector<std::vector<double> > Aloc, std::vector<double> bloc, std::vector<double> x0loc , double err, int kmax, int nx, int ny)
{
  // // Algorithme du gradient conjugué parallèle qui prend en argument uniquement des vecteurs locaux et renvoie un vecteur local.

  int k;

  double norm_r, nr_carre, nr2_carre;
  std::vector<double> wloc, rloc, r_1loc, ploc, dloc, xloc;

  int Np = nx*ny;

  wloc.resize(Np); rloc.resize(Np); ploc.resize(Np); dloc.resize(Np);

  xloc = x0loc;


  // wloc = prodMVC(Aloc,xloc,nx,ny);
  wloc = prodM5V(Aloc, xloc, nx, ny);

  for (int i = 0; i < Np; i++)
  {
    rloc[i] = bloc[i] - wloc[i];
  }
  ploc = rloc;

  k = 0;

  nr_carre = dot(rloc,rloc);
  norm_r = sqrt(nr_carre);


  while ((norm_r > err) and (k<kmax))
    {
      // dloc = prodMVC(Aloc,ploc,nx,ny);
      dloc = prodM5V(Aloc, ploc, nx, ny);

      double alpha = nr_carre/dot(ploc,dloc);

      for (int i = 0; i < Np; i++)
      {
        xloc[i] += alpha*ploc[i];
        rloc[i] -= alpha*dloc[i]; //rk devient r(k+1)
      }

      nr2_carre = dot(rloc,rloc); //norme de r(k+1)

      double beta = nr2_carre/nr_carre;

      for (int i = 0; i < Np; i++)
      {
        ploc[i] = rloc[i] + beta*ploc[i];
      }

      nr_carre = nr2_carre;
      norm_r = sqrt(nr_carre);

      k++;
    }

  return xloc;

}
