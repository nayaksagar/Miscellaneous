#include<iostream>
#include<vector>
#include<array>
#include<cmath>
#include<fstream>
#include<cstring>
#include<chrono>
#include<numeric>

#define nx 50
#define ny 50
#define Q 9  //D2Q9

// Global variables
// double f[9][nx][ny], f_eq[9][nx][ny], f_stream[9][nx][ny], f_coll[9][nx][ny];
double rho0=1.0;
double cssq = 1./3;
double tau = 0.933;
double nu=1.0/3 * (tau-0.5);
double sum_rho = 0.0, sum_ux = 0.0, sum_uy = 0.0;
double Fx = 1.0e-5, Fy = 0.0;
double L2_ux = 1.0, L2_uy = 1.0, tol = 1.0e-6;
double uw = 0.0;
double dx = 1.0, dy = 1.0;
// Discrete velocities and corresponding weights of D2Q9 set
int cx[9]={0, 1, 0, -1, 0, 1, -1, -1, 1};
int cy[9]={0, 0, 1, 0, -1, 1, 1, -1, -1};
double w[9] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};

class Lattice{
  public:
    std::array<double,2> pos;
    double rho, T;
    std::array<double,2> vel, vel_old;
    std::array<double,Q> f, f_eq, f_stream, f_coll;
    int tag;

};

// Functions
void equilibrium();
void collision();
void streaming();
void Boundary_conditions();

inline int indx(int i, int j){return j*nx + i;}

std::vector<Lattice> lattice;

int main(){
  for(int i=0;i<nx;++i){
      double x = i*dx;
      for(int j=0;j<ny;++j){
        double y = j*dy;
        
        Lattice node;
        node.pos[0] = x;
        node.pos[1] = y;
        node.rho = rho0;
        node.vel[0] = 0.0;
        node.vel[1] = 0.0;
        for(int a = 0; a < Q; a++)
          node.f[a] = w[a] * node.rho;

        if(i==0 && j==0)
          node.tag = 0;
        else if(i==nx-1 && j==0)
          node.tag = 1;
        else if(i==nx-1 && j==ny-1)
          node.tag = 2;
        else if(i==0 && j==ny-1)
          node.tag = 3;
        else if(i==0)
          node.tag = 4;
        else if(j==0)
          node.tag = 5;
        else if(i==nx-1)
          node.tag = 6;
        else if(j==ny-1)
          node.tag = 7;
        else
          node.tag = 8;

        lattice.push_back(node);     
      }
  }

  auto start = std::chrono::high_resolution_clock::now();


  int t=1;
  while (t<100){
    // Macroscopic variables
    double Err_ux = 0.0, Err_uy = 0.0;
    double sq_ux = 0.0, sq_uy = 0.0;
    
    for(int i = 0; i < lattice.size(); i++){
      sum_rho = 0.0;
      sum_ux = 0.0;
      sum_uy = 0.0;

      for(int a = 0; a < Q; a++){
        sum_rho += lattice[i].f[a];
        sum_ux  += lattice[i].f[a]*cx[a];
        sum_uy  += lattice[i].f[a]*cy[a];
      }

      //save old velocity
      lattice[i].vel_old[0] = lattice[i].vel[0];
      lattice[i].vel_old[1] = lattice[i].vel[1];

      lattice[i].rho = sum_rho;
      lattice[i].vel[0] = (sum_ux + 0.5*Fx)/lattice[i].rho;
      lattice[i].vel[1] = (sum_uy + 0.5*Fy)/lattice[i].rho;
      
      // Error calculation
      Err_ux += pow((lattice[i].vel_old[0] - lattice[i].vel[0]),2);
      Err_uy += pow((lattice[i].vel_old[1] - lattice[i].vel[1]),2);
      sq_ux += pow(lattice[i].vel_old[0],2);
      sq_uy += pow(lattice[i].vel_old[1],2);
    }
    
    L2_ux = sqrt(Err_ux/sq_ux);
    // if(t%1000==0){

      // std::cout<<t<<"\t"<<L2_ux<<"\t"<<std::endl;
    // }
    t++;

    equilibrium();
    collision();
    streaming();
    Boundary_conditions();
  }

  double analytical;
  std::ofstream lineplot;
  lineplot.open("lineplot.dat");
  double delta = 0.5;
  double y;
  int H = ny-1;
  for (int j=1;j<ny;j++){
    y = j-delta;
    analytical = (0.5*Fx*y*(H-y)/nu) + uw*((y/H)-0.5);
    lineplot<<j<<"\t"<<lattice[indx(nx/2,j)].vel[0]<<"\t"<<analytical<<std::endl;
  }
  lineplot.close();

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;
  return 0;

}

void equilibrium(){
  double term1, term2, term3;
  for(int i = 0; i < lattice.size(); i++){
    for (int a = 0;a < Q; a++){
      term1 = (lattice[i].vel[0]*cx[a] + lattice[i].vel[1]*cy[a])/cssq;
      term2 = 0.5*pow(term1,2);
      term3 = 0.5*(lattice[i].vel[0]*lattice[i].vel[0] + lattice[i].vel[1]*lattice[i].vel[1])/cssq;
      lattice[i].f_eq[a] = w[a]*lattice[i].rho*(1.0 + term1 + term2 - term3);
    }
  }
}

void collision(){
  double source;
  for (int i = 0; i < lattice.size(); i++){
    for (int a = 0; a < Q; a++){
      source = (1.0-0.5/tau) * w[a] * ( (cx[a]-lattice[i].vel[0]) /cssq + (cx[a]*(cx[a]*lattice[i].vel[0]+cy[a]*lattice[i].vel[1])) /(cssq*cssq)) * Fx;
      lattice[i].f_coll[a] = lattice[i].f[a] - (lattice[i].f[a] - lattice[i].f_eq[a])/tau + source;
      std::cout << a << "\t" << lattice[i].f_coll[a] << std::endl; 
    }
  }
}

void streaming(){
  for (int i = 0; i < nx; i++){
    for (int j = 0; j < ny; j++){
      int ia, ja;
      for (int a = 0; a < Q; a++){
        ia = i + cx[a];
        ja = j + cy[a];
        if (ia<0){           // periodic
          ia = nx-1;
        }
        if (ia>nx-1){
          ia = 0;
        }
        lattice[indx(i,j)].f[a] = lattice[indx(i,j)].f_coll[a];
        // f[a][ia][ja] = f_coll[a][i][j];
      }

    }
  }
}

void Boundary_conditions(){
  for (int i=0;i<nx;i++){
    // bottom wall
    int j = 1;
    lattice[indx(i,j)].f[2]= lattice[indx(i,j)].f_coll[4];
    lattice[indx(i,j)].f[5]= lattice[indx(i,j)].f_coll[7] + (uw*cx[5]/12.);
    lattice[indx(i,j)].f[6]= lattice[indx(i,j)].f_coll[8] + (uw*cx[6]/12.);

    // top wall
    j = ny-1;
    lattice[indx(i,j)].f[4] = lattice[indx(i,j)].f[2];
    lattice[indx(i,j)].f[7] = lattice[indx(i,j)].f[5] + (uw*cx[7]/12.);
    lattice[indx(i,j)].f[8] = lattice[indx(i,j)].f[6] + (uw*cx[8]/12.);
  }
}
