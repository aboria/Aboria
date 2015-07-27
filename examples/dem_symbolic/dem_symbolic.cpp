/*
 * dem_symbolic.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */

#include <random>
typedef std::mt19937 generator_type;
generator_type generator;

#include "Aboria.h"
using namespace Aboria;

#include "Visualisation.h"


int main(int argc, char **argv) {
    const double tol = GEOMETRY_TOLERANCE;

    ABORIA_VARIABLE(radius,double,"radius")
    ABORIA_VARIABLE(mass,double,"mass")
    ABORIA_VARIABLE(velocity,Vect3d,"velocity")
    ABORIA_VARIABLE(acceleration,Vect3d,"acceleration")

    typedef Particles<radius,velocity,acceleration,mass> dem_type;
    dem_type dem;

    const int timesteps = 30000;
    const int nout = 200;
    const int timesteps_per_out = timesteps/nout;
    const double L = 31.0/1000.0;
    const int N = 1000;


     /* dem parameters
     */
    const double dem_diameter = 0.0011;
    const double dem_gamma = 0.0004;
    const double dem_k = 1.0e01;
    const double dem_dens = 1160.0;
    const double dem_mass_min = (1.0/6.0)*PI*pow(0.5*dem_diameter,3)*dem_dens;
    const double dem_mass_max = (1.0/6.0)*PI*pow(dem_diameter,3)*dem_dens;
    const double dem_min_reduced_mass = dem_mass_min*dem_mass_max/(dem_mass_min+dem_mass_max);
    const double dt = (1.0/50.0)*PI/sqrt(dem_k/dem_min_reduced_mass-pow(0.5*dem_gamma/dem_min_reduced_mass,2));

    dem.init_neighbour_search(Vect3d(0,0,-dem_diameter),Vect3d(L/3,L/3,L+dem_diameter),dem_diameter,Vect3b(true,true,false));

    std::uniform_real_distribution<double> uni(0,1);
    for (int i = 0; i < N; ++i) {
        bool free_position = false;
        dem_type::value_type p(Vect3d(uni(generator)*L/3,uni(generator)*L/3,uni(generator)*(L-dem_diameter)+dem_diameter/2));
        set<radius>(p,(uni(generator)+1)*dem_diameter/4);
        while (free_position == false) {
            set<position>(p,Vect3d(uni(generator)*L/3,uni(generator)*L/3,uni(generator)*(L-dem_diameter)+dem_diameter/2));
            free_position = true;
            for (auto tpl: dem.get_neighbours(get<position>(p))) {
                const Vect3d& dx = std::get<1>(tpl);
                const dem_type::value_type& j = std::get<0>(tpl);
                if (norm(dx) < get<radius>(j) + get<radius>(p)) {
                    free_position = false;
                    break;
                }
            }
        }
        dem.push_back(p);
    }


    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    dem.copy_to_vtk_grid(grid);
    Visualisation::vtkWriteGrid("vis/demInit",0,grid);

    auto p = get_vector<position>(dem);
    auto r = get_vector<radius>(dem);
    auto m = get_vector<mass>(dem);
    auto v = get_vector<velocity>(dem);
    auto dvdt = get_vector<acceleration>(dem);
    auto id_ = get_vector<id>(dem);


    Label<0> a;
    Label<1> b;
    Dx dx;
    Accumulate<std::plus<Vect3d> > sum;
    
    /*
     * timestepping
     */
    v = 0;
    m = (1.0/6.0)*PI*8*r*r*r*dem_dens;
    for (int io = 0; io < nout; ++io) {
        std::cout << "." << std::flush;
        for (int i = 0; i < timesteps_per_out; i++) {
            p += v*dt;

            dvdt = (// spring force between dem particles
                    sum(b=dem, id_[a]!=id_[b] && norm_(dx)<r[a]+r[b], 
                          dem_k*((r[a]+r[b])/norm_(dx)-1)*dx  + dem_gamma*(v[b]-v[a]))
                
                    // spring force between particles and bottom wall
                    + if_else(r-p[2] > 0, dem_k*(r-p[2]), 0.0)*Vect3d(0,0,1) 

                    // gravity force
                    + Vect3d(0,0,-9.81)*m

                    )/m;

            v += dvdt*dt;
        }
        dem.copy_to_vtk_grid(grid);
        Visualisation::vtkWriteGrid("vis/dem",io,grid);
    }
    std::cout << std::endl;
}
