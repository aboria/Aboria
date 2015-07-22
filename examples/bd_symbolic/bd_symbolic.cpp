/*
 * bd_symbolic.cpp
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

    Particles<radius> spheres;

    const double L = 10.0;
    const double D = 1.0;
    const double dt = 0.1;
    const double timesteps = 1000;

    spheres.push_back(Vect3d(0,0,0));
    spheres[0].set<radius>(1.0);
    spheres.push_back(Vect3d(5,0,0));
    spheres[1].set<radius>(2.0);
    spheres.push_back(Vect3d(0,-5,0));
    spheres[2].set<radius>(1.5);
    spheres.push_back(Vect3d(0,0,5));
    spheres[3].set<radius>(1.0);

    Particles<> points;
    std::uniform_real_distribution<double> uni(-L/5,L/5);
    for (int i = 0; i < 10000; ++i) {
        points.push_back(Vect3d(uni(generator),uni(generator),uni(generator)));
    }


    points.init_neighbour_search(Vect3d(-L/5,-L/5,-L/5),Vect3d(L/5,L/5,L/5),4,Vect3b(true,true,true));
    spheres.init_neighbour_search(Vect3d(-L,-L,-L),Vect3d(L,L,L),4,Vect3b(false,false,false));

    auto spheres_position = get_vector<position>(spheres);
    auto spheres_radius = get_vector<radius>(spheres);

    auto points_position = get_vector<position>(points);
    auto points_alive = get_vector<alive>(points);


    Label<0> a;
    Label<1> b;
    Dx dx;
    Normal N;
    GeometriesSymbolic<Sphere> spheres_;        
    VectorSymbolic<double> vector;      

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    int count_before=0;
    for(auto i: points) {
        if (norm(get<position>(i)-get<position>(spheres[0])) < get<radius>(spheres[0])) {
            count_before++;
        }
    }

    /*
     * Kill any points within spheres
     */
    points_alive = !sum_(b=spheres, norm_(dx) < spheres_radius[b],true,false);


    int count_after=0;
    for(auto i: points) {
        if (norm(get<position>(i)-get<position>(spheres[0])) < get<radius>(spheres[0])) {
            count_after++;
        } 
    }

    std::cout <<" found "<<count_before<<" before and "<<count_after<<" after"<<std::endl;

    points.copy_to_vtk_grid(grid);
    Visualisation::vtkWriteGrid("vis/pointsInit",0,grid);

    /*
     * Diffusion step for points and reflect off spheres
     */
    for (int i = 1; i < timesteps; ++i) {
        if (i%10==0) {
            std::cout << "." << std::flush;
            points.copy_to_vtk_grid(grid);
            Visualisation::vtkWriteGrid("vis/points",i/10,grid);
        }
        points_position += std::sqrt(2*D*dt)*vector(N,N,N) | spheres_(b=spheres,spheres_radius[b]);
    }
    std::cout << std::endl;
}
