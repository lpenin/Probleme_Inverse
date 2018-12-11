#include "fonction.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include "spline.h"

using namespace std;


int main(int argc, char * argv[])
{
  // PARAMETRE
  int Nx = 50;
  int Ny = 50;
  int Np = Nx*Ny;
  int M=5;
  double xmin = 0.;
  double xmax = 1.;
  double ymin = 0.;
  double ymax = 1.;
  double hx = (xmax-xmin)/(Nx+1);
  double hy = (ymax-ymin)/(Ny+1);

  double D = 1.;
  double deltaT =0.9*hx*hy/(4.*D);
  double tfinal = 1.;
  int Nt = int (tfinal/deltaT) +1;

  double Coef = D;
  double eps = 0.0001;
  double norm;

  //POINT D'INTERET U*
  std::vector<int> Xpos(M);
  std::vector<double> Xvalues(M);
  std::vector<double> Yvalues(M);

  for (int i=0; i<M; i++)
  {
    Xvalues[i] =(i+1)*(xmax-xmin)/(M+1);
    int k=0;
    while (k*hy < Xvalues[i])
    {
      k++;
    }
    Xpos[i]=k;
  }
  int Ypos = 10;
  cout << "calcul value ok"<< endl;

  // INITIALISATION DE LA MATRICE

  double alpha = 1 + 2*D*deltaT/(hx*hx) + 2*D*deltaT/(hy*hy);
  double beta = -D*deltaT/(hx*hx);
  double gamma = -D*deltaT/(hy*hy);

  vector<vector <double> > LapMat(5);

  // LapMat.resize(Nx*Ny);
  // for (int i=0; i<Nx*Ny; i++)
  // {
  //     LapMat[i].resize(Nx*Ny);
  // }
  //
  // for (int i=0; i<Nx*Ny; i++)
  // {
  //   for (int j=0; j<Ny*Nx; j++)
  //   {
  //     if (i==j)
  //     {
  //       LapMat[i][j]=alpha;
  //     }
  //     if (((i==j+1)&&(i%(Nx-1)!=0)) || ((j=i+1)&&(i%(Nx+1)))) LapMat[i][j]=gamma;
  //
  //     if ((i==j+2*Nx) || (j==i+2*Nx)) LapMat[i][j]=beta;
  //
  //     cout << " j atteint = "<<j<<endl;
  //   }
  //   cout << " i atteint = "<<i<<endl;
  // }
  // cout << "fin initialisation" << endl;

  for (int i=0; i<5; i++)
  {
    LapMat[i].resize(Np);
  }

  for (int i=0; i<Np; i++)
  {
    LapMat[0][i]=gamma;
    LapMat[1][i]=beta;
    LapMat[2][i]=alpha;
    LapMat[3][i]=beta;
    LapMat[4][i]=gamma;

    if (i%Nx==0)
    LapMat[1][i]=0.;

    if (i%Nx==Nx-1)
    LapMat[3][i]=0.;

    if(i<Nx)
    LapMat[0][i]=0.;

    if(i>(Ny-1)*Nx-1)
    LapMat[4][i]=0.;
  }
  cout << "fin initialisation matrice" << endl;



  // CONDITIONS INITIALS
  vector<double> sol(Np);
  for (int i=0; i<Np;i++)
  {
    sol[i] = 0.;
  }

  // CONDITIONS limites
  for (int j = 0; j < Nx ; j++) //HAUT
  {
    sol[j] = sol[j+Nx];
  }

  for (int j = Nx*(Ny-1); j < Nx*Ny ; j++) //BAS
  {
    sol[j] = sol[j-Nx];
  }

  for (int i = 0; i < Nx*(Ny-1); i=i+Nx) //GAUCHE
  {
    sol[i] = 0;
  }

  for (int i = Nx-1; i < Ny*Nx; i=i+Nx) //DROITE
  {
    sol[i] = (1-i*hy)*i*hy; //TODO fonction g_th
  }
  cout << "fin initialisation CL" << endl;

  // Calcul Solution
  double t=0.;
  for( int it=0 ; it<Nt; it++)
  {
    int kmax = Nx*Ny +100; //Pour une matrice de taille n le GC met max n étapes en théorie,
    sol = CG(LapMat,sol,sol,0.001,kmax,Nx,Ny);
    // CONDITIONS limites
    for (int j = 0; j < Nx ; j++) //HAUT
    {
      sol[j] = sol[j+Nx];
    }

    for (int j = Nx*(Ny-1); j < Nx*Ny ; j++) //BAS
    {
      sol[j] = sol[j-Nx];
    }

    for (int i = 0; i < Nx*(Ny-1); i=i+Nx) //GAUCHE
    {
      sol[i] = 0;
    }

    for (int i = Nx-1; i < Ny*Nx; i=i+Nx) //DROITE
    {
      sol[i] = (1-i*hy)*i*hy; //TODO fonction g_th
    }
    t+= deltaT;
    // cout << "temps atteint " << t << endl;
  }
  cout << "fin calcul solution" << endl;

    // RECUPERATION DES DONNEES UTILES
    std::vector<double> Uth(M);
    for (int i=0; i<M; i++)
    {
      int position =Xpos[i]*Nx+Ypos;
      Uth[i] = sol[position];
    }
    cout << "fin recuperation données" << endl;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUT PROBLEME INVERSE !!!!!!!!!!!!!!!!!!!!


    // SOURCE CORRESPONDAND AU Uexp
    std::vector<double> Uexp(M);
    norm = 0.;
    for (int i=0; i<M; i++)
    {
      Uexp[i]=0.;
      Yvalues[i]=0.;
      norm = norm + abs (Uth[i]-Uexp[i]);
    }
    tk::spline gexp;
    gexp.set_points(Xvalues,Yvalues);
    cout << "fin interpolation" << endl;
    for (int i = Nx-1; i < Ny*Nx; i=i+Nx) //DROITE
    {
      sol[i] = gexp(i*hy/(Nx+1));
      cout << i*hy/(Nx+1) << "     " << gexp(i*hy) << endl;
    }


    // CONDITIONS INITIALS
    for (int i=0; i<Np;i++)
    {
      sol[i] = 0.;
    }

    int ik=0;
    int kmax = 1;
    while (norm > eps || ik<kmax)
    {
      // Calcul Solution
      for( int it=0 ; it<Nt; it++)
      {
        int kmax = Nx*Ny +100; //Pour une matrice de taille n le GC met max n étapes en théorie,
        sol = CG(LapMat,sol,sol,0.001,kmax,Nx,Ny);
        // CONDITIONS limites
        for (int j = 0; j < Nx ; j++) //HAUT
        {
          sol[j] = sol[j+Nx];
        }

        for (int j = Nx*(Ny-1); j < Nx*Ny ; j++) //BAS
        {
          sol[j] = sol[j-Nx];
        }

        for (int i = 0; i < Nx*(Ny-1); i=i+Nx) //GAUCHE
        {
          sol[i] = 0;
        }

        for (int i = Nx-1; i < Ny*Nx; i=i+Nx) //DROITE
        {
          sol[i] = gexp(i*hy/(Nx+1));
        }
      }
      cout << "fin calcul exp "<< ik  << endl;
      ik++;

      // RECUPERATION DES DONNEES UTILES
      norm = 0.;
      for (int i=0; i<M; i++)
      {
        Uexp[i]=sol[Xpos[i]*Nx+Ypos];
        Yvalues[i]=Yvalues[i] + Coef * (Uth[i] - Uexp[i]);
        norm = norm + abs (Uth[i]-Uexp[i]);
      }
      tk::spline gexp;
      gexp.set_points(Xvalues,Yvalues);
      cout << "fin calcul norm "<< norm << " et interpolation"  << endl;

    }


    return 0;
  }
