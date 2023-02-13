// ============================================================
// This is a event-driven molecular dynamics simulations of
// hard spheres particles in a box with periodic boundary
// the mases are the same. 
//============================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;

#define PI 3.1415926535897932



int main(){

    //------------------------------------------------------------
    //Declaration of variables
    double x0,y0,z0,v0x,v0y,v0z,t0,tf,dt,t;
    double Lx, Ly, Lz, R, alpha, omega, Nx, Ny, Nz;
    int  i, N;
    string buf;

    //------------------------------------------------------------
    // NOTE: if we desire the same random number at each run, 
    // use the following two lines

    // default_random_engine generator;
    // uniform_real_distribution<double> dist(-1.0,1.0);
    //------------------------------------------------------------

   //------------------------------------------------------------
    // if we desire different random numbers at each run, 
    // use the following 3 lines

    random_device generator;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(generator()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> dist(-1.0,1.0);
 
    //------------------------------------------------------------
    //Ask user for input:

    cout << "# Enter R, Nx, Ny, Nz:\n";
    cin  >> R >> Nx >> Ny >> Nz;           getline(cin,buf);

    cout << "# Enter t0,tf:\n";
    cin  >> t0 >> tf;  getline(cin,buf);

    cout << "# Enter alpha, omega:\n";
    cin  >> alpha >> omega;  getline(cin,buf);


    // cout << "# Enter Lz:\n";
    // cin  >> Lz ;  getline(cin,buf);

    //------------------------------------------------------------
    
    N = Nx*Ny*Nz;

    
    //Lz = 29*(2*R);
    
    // note: for Lz = 8*(2*R) eight particles are lined in the vertical direction touching each other 
    // and touching both wall which cause an issued in the collision time calculation. To overcome 
    // this we have to modify the initial number of particles in x, y, and z. This means that 
    // Nx, Ny and Nz must change to fix the number of particles. This is valid for Lz <= 8*(2*R)
    Lz = 2*(2*R);  
    
    double n = 0.02/pow(2*R,3);                 // particle density
    double phi = n*(4.0/3.0)*PI*pow(R, 3);
    Lx = sqrt((N)/((Lz-2*R)*n));
    Ly = sqrt((N)/((Lz-2*R)*n));
                               
    double r[N][3];
    double v[N][3];

    double d_x =  (Lx - 2*R)/(1.0*Nx);        
    double d_y =  (Ly - 2*R)/(1.0*Ny);      
    double d_z =  (Lz - 2*R)/(1.0*Nz); 

    int np = 0;

    for (int k=1; k<=Nz; k++)
    {    
        

        z0 = k*d_z;
        
        for (int i=1; i<=Ny; i++)
        {  
            y0 = i*d_y;

            for (int j=1; j<=Nx; j++)
            {
                x0 = j*d_x;

                r[np][0] = x0;
                r[np][1] = y0;
                r[np][2] = z0;

                v[np][0] = 0.5*dist(generator);
                v[np][1] = 0.5*dist(generator);
                v[np][2] = 0.5*dist(generator);

                np = np + 1;

            }
        }
    }
   
    // if( R <= 0.0){cerr << "R <=0\n"; exit(1);}
    // if( L <= 0.0){cerr << "L <=0\n"; exit(1);}
    // if( z0 <  0.0){cerr << "z0<=0\n"; exit(1);}
    // if( z0 >  L ){cerr << "z0> L\n"; exit(1);}
    // if( rxy > R){cerr << "rxy > R\n"; exit(1);}
    // if( v0x*v0x + v0y*v0y + v0z*v0z <  0.0){cerr << "v0=0\n"; exit(1);}


    t = t0;
    int count = 0;
    int Iteration = 0;
    //int increment = 100;        // for H=29sigma
    // int increment = 1500;        // for H=8sigma
    int increment = 10000;        // for H=8sigma
    
    int file=0.0;
    
    int nc;
    double tij1;
    double tij2;
    int indexi;
    int indexj;

    double rxij, ryij, rzij;
    double vxij, vyij, vzij;
    double dvdr, dr2, dv2, bij;



    while( t <= tf )
    {   

        if ( count == Iteration )
        {   
            cout << "t = " << t << endl;

            ofstream myfile("3Dgas_"+ to_string(file) + ".dat" );
            myfile.precision(16);
           
            myfile << "t" << " " << "N" << " \t " << "Lx" << " \t\t " << "Ly" << " \t " << "Lz" << " \t " << "R" << endl; 
            myfile << t << " " <<  N << " " << Lx << " " << Ly << " " << Lz << " " << R << endl; 
            
            myfile << "x" << " \t\t " << "y" << " \t\t\t " << "z" << " \t\t "
            << "vx" << " \t\t " << "vy" << " \t\t\t " << "vz" << endl;

            for (int m = 0; m < N; m++)
            {   
                myfile << r[m][0] << " " << setw(15) << r[m][1] << setw(9) << " " << r[m][2] << setw(9) << " "
                    << v[m][0] << " " << setw(18) << v[m][1] << setw(9) << " " << v[m][2] << endl; 
            }   

            myfile.close();
                        
            Iteration = Iteration + increment;
            file = file + 1;

            cout << file;

        }

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // COLLISION TIME CALCULATION 
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        tij1 = pow(10,5);
        
        for(int a = 0; a < N; a++)
        {   
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // COLLISION TIME WITH WALLS
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            // z-direction bottom wall
            if( v[a][2] < 0) 
            {
                tij2 = (R - r[a][2])/v[a][2];

                if ( tij2 < tij1 )
                {                            
                    tij1 = tij2;
                    indexi = a;
                    indexj = N;
                }
            }

            //z-direction top wall
            if( v[a][2] > 0) 
            {
                tij2 = (Lz - R - r[a][2])/v[a][2];

                if ( tij2 < tij1 )
                {
                    tij1 = tij2;
                    indexi = a;
                    indexj = N + 1;
                }
            }       
            
            
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // COLLISION TIME BETWEEN PARTICLES
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            for(int b = a + 1; b < N; b++)
            {
                rxij =  r[a][0] - r[b][0];
                ryij =  r[a][1] - r[b][1];
                rzij =  r[a][2] - r[b][2];

                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // MINIMUM IMAGE CONVENTION CALCULATION
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                rxij = rxij - Lx*round(rxij/Lx);
                ryij = ryij - Ly*round(ryij/Ly);

                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                vxij = v[a][0] - v[b][0];
                vyij = v[a][1] - v[b][1];
                vzij = v[a][2] - v[b][2];

                dr2 =  rxij*rxij + ryij*ryij + rzij*rzij;
                dv2 =  vxij*vxij + vyij*vyij + vzij*vzij; 

                dvdr = rxij*vxij + ryij*vyij + rzij*vzij;
                
                if ( dvdr < 0.0 ) // this ensures that particles collide
                {   
                    bij = dvdr*dvdr - dv2*(dr2 - 4*R*R);

                    if ( bij > 0.0 )
                    {     
                        tij2 = (- dvdr - sqrt(bij))/dv2;  //positive sign before sqrt gives overlap of particles, no possible for hard spheres.
                        
                        if ( tij2 < tij1 )
                        {   
                            tij1 = tij2;
                            indexi = a;
                            indexj = b;   
                        }
                    }
                }
            }
        }

        count++;

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //MOVE ALL PARTICLES AT THE MINIMUM TIME FOUND
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        for(int a = 0; a < N; a++)
        {
            for(int j = 0; j < 3; j++) 
            {    r[a][j] = r[a][j] + v[a][j]*tij1; }

            //+++++++++++++++++++++++++++++++++++++++++++++++
            // APPLY PBC
            //+++++++++++++++++++++++++++++++++++++++++++++++

            r[a][0] = r[a][0] - Lx*round((r[a][0] - 0.5*Lx)/Lx);
            r[a][1] = r[a][1] - Ly*round((r[a][1] - 0.5*Ly)/Ly);

        }
            
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //APPLY COLLISION LAW OF PARTICLE INDEX WITH WALLS IN Z 
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if( indexj == N )     v[indexi][2] = - v[indexi][2] + 2*omega;
        if( indexj == N + 1 ) v[indexi][2] = - v[indexi][2] - 2*omega;
                                    
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //  APPLY COLLISION LAW TO PARTICLES IN CONTACT
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if ( indexj < N )
        {  

            rxij = r[indexi][0] - r[indexj][0];
            ryij = r[indexi][1] - r[indexj][1];
            rzij = r[indexi][2] - r[indexj][2];

            rxij = rxij - Lx*round(rxij/Lx);
            ryij = ryij - Ly*round(ryij/Ly);
            
            vxij = v[indexi][0] - v[indexj][0];
            vyij = v[indexi][1] - v[indexj][1];
            vzij = v[indexi][2] - v[indexj][2];

            dvdr = rxij*vxij + ryij*vyij + rzij*vzij; 

            v[indexi][0] = v[indexi][0] - 0.5*(1 + alpha)*(dvdr/(2*R))*(rxij/(2*R));
            v[indexj][0] = v[indexj][0] + 0.5*(1 + alpha)*(dvdr/(2*R))*(rxij/(2*R));
            v[indexi][1] = v[indexi][1] - 0.5*(1 + alpha)*(dvdr/(2*R))*(ryij/(2*R));
            v[indexj][1] = v[indexj][1] + 0.5*(1 + alpha)*(dvdr/(2*R))*(ryij/(2*R));
            v[indexi][2] = v[indexi][2] - 0.5*(1 + alpha)*(dvdr/(2*R))*(rzij/(2*R));
            v[indexj][2] = v[indexj][2] + 0.5*(1 + alpha)*(dvdr/(2*R))*(rzij/(2*R));
        }        
                
        t = t + tij1;

        // cout << "t = " << t << endl;

        // cout << "nc = " << count << endl;

    }

}
