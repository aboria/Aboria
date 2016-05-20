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

#ifdef HAVE_VTK
#include "Visualisation.h"
#endif

#include <boost/math/constants/constants.hpp>

const double PI = boost::math::constants::pi<double>();


int main(int argc, char **argv) {
    ABORIA_VARIABLE(radius,double,"radius")
    ABORIA_VARIABLE(mass,double,"mass")
    ABORIA_VARIABLE(velocity,double3,"velocity")
    ABORIA_VARIABLE(acceleration,double3,"acceleration")

    typedef Particles<std::tuple<radius,velocity,acceleration,mass>> dem_type;
    typedef position_d<3> position;
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
    const double dem_mass_min = (1.0/6.0)*PI*std::pow(0.5*dem_diameter,3)*dem_dens;
    const double dem_mass_max = (1.0/6.0)*PI*std::pow(dem_diameter,3)*dem_dens;
    const double dem_min_reduced_mass = dem_mass_min*dem_mass_max/(dem_mass_min+dem_mass_max);
    const double dt = (1.0/50.0)*PI/sqrt(dem_k/dem_min_reduced_mass-std::pow(0.5*dem_gamma/dem_min_reduced_mass,2));

    dem.init_neighbour_search(double3(0,0,-dem_diameter),double3(L/3,L/3,L+dem_diameter),dem_diameter,bool3(true,true,false));

    std::uniform_real_distribution<double> uni(0,1);
    for (int i = 0; i < N; ++i) {
        bool free_position = false;
        dem_type::value_type p;
        set<position>(p,double3(uni(generator)*L/3,uni(generator)*L/3,uni(generator)*(L-dem_diameter)+dem_diameter/2));
        set<radius>(p,(uni(generator)+1)*dem_diameter/4);
        while (free_position == false) {
            set<position>(p,double3(uni(generator)*L/3,uni(generator)*L/3,uni(generator)*(L-dem_diameter)+dem_diameter/2));
            free_position = true;
            for (auto tpl: dem.get_neighbours(get<position>(p))) {
                const double3& dx = std::get<1>(tpl);
                const dem_type::value_type& j = std::get<0>(tpl);
                if (dx.norm() < get<radius>(j) + get<radius>(p)) {
                    free_position = false;
                    break;
                }
            }
        }
        dem.push_back(p);
    }

#ifdef HAVE_VTK
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    dem.copy_to_vtk_grid(grid);
    Visualisation::vtkWriteGrid("vis/demInit",0,grid);
#endif

    Symbol<position> p;
    Symbol<radius> r;
    Symbol<mass> m;
    Symbol<velocity> v;
    Symbol<acceleration> dvdt;
    Symbol<id> id_;
    Label<0,dem_type> a(dem);
    Label<1,dem_type> b(dem);
    Dx dx;
    Accumulate<std::plus<double3> > sum;
    
    /*
     * timestepping
     */
    v[a] = 0;
    m[a] = (1.0/6.0)*PI*8*r*r*r*dem_dens;
    for (int io = 0; io < nout; ++io) {
        std::cout << "." << std::flush;
        for (int i = 0; i < timesteps_per_out; i++) {
            p[a] += v*dt;

            dvdt[a] = (// spring force between dem particles
                    sum(b, id_[a]!=id_[b] && norm(dx)<r[a]+r[b], 
                          dem_k*((r[a]+r[b])/norm(dx)-1)*dx  + dem_gamma*(v[b]-v[a]))
                
                    // spring force between particles and bottom wall
                    + if_else(r-p[2] > 0, dem_k*(r-p[2]), 0.0)*double3(0,0,1) 

                    // gravity force
                    + double3(0,0,-9.81)*m

                    )/m;

            v[a] += dvdt*dt;
        }
#ifdef HAVE_VTK
        dem.copy_to_vtk_grid(grid);
        Visualisation::vtkWriteGrid("vis/dem",io,grid);
#endif
    }
    std::cout << std::endl;
}
