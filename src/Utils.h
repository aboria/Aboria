/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef UTILS_H_
#define UTILS_H_

#ifdef HAVE_VTK
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

#include "Aboria.h"
#include <boost/math/constants/constants.hpp>

#ifdef HAVE_CAIRO
#include <cairo-svg.h>
#endif

#include <random>

namespace Aboria {

/*
#ifdef HAVE_THRUST
template <typename T>
struct is_device_reference {
    typedef std::false_type type;
    typedef T value_type;
};

template <typename T>
struct is_device_reference<thrust::device_reference<T>>{
    typedef std::true_type type;
    typedef T value_type;
};

template <typename T>
typename is_device_reference<typename
std::remove_reference<T>::type>::value_type copy_to_host_impl(T&&
arg,std::true_type) { return typename is_device_reference<typename
std::remove_reference<T>::type>::value_type(arg);
}

template <typename T>
auto copy_to_host_impl(T&& arg,std::false_type) -> decltype(arg) {
    return arg;
}

template <typename T>
auto copy_to_host(T&& arg)
-> decltype(copy_to_host_impl(arg,typename is_device_reference<typename
std::remove_reference<T>::type>::type())) { return
copy_to_host_impl(arg,typename is_device_reference<typename
std::remove_reference<T>::type>::type());
}
#else
template <typename T>
auto copy_to_host(T&& arg) -> decltype(arg) {
    return arg;
}
#endif
*/

#ifdef HAVE_VTK
/// Write a vtk unstructured grid to a filename of the form:
/// printf("name%05d",timestep). User can optionally provide \p fields, which is
/// a set of key/value pairs (e.g. "time"->0.1) written as additional fields to
/// the file
///
/// \param name the base filename
/// \param timestep the index in the filename
/// \param grid the grid to write to the file
/// \param fields a `std::map` of additional fields to write
///
template <typename T = void>
void vtkWriteGrid(const char *name, int timestep,
                  vtkSmartPointer<vtkUnstructuredGrid> grid,
                  const std::map<std::string, double> &fields = {}) {
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

#if VTK_MAJOR_VERSION > 5
  writer->SetInputData(grid);
#else
  writer->SetInput(grid);
#endif

  for (auto &i : fields) {
    vtkFloatArray *t = vtkFloatArray::New();
    t->SetName(i.first.c_str());
    t->SetNumberOfTuples(1);
    t->SetTuple1(0, i.second);
    grid->GetFieldData()->AddArray(t);
  }

  writer->SetDataModeToBinary();

  char buffer[100];
  sprintf(buffer, "%s%05d.vtu", name, timestep);

  writer->SetFileName(buffer);
  writer->Write();
}
#endif

template <typename PTYPE, typename GTYPE, typename RNDGTYPE>
void create_n_particles_with_rejection_sampling(const unsigned int n,
                                                PTYPE &particles,
                                                const GTYPE &generator_fun,
                                                const Vector<double, 3> &low,
                                                const Vector<double, 3> &high,
                                                RNDGTYPE &generator) {

  // assume 3d for now...
  typedef Vector<double, 3> Vect3d;

  std::uniform_real_distribution<double> uni(0, 1);
  const Vect3d h = (high - low) / 1000;
  double max_generator_fun = 0;
  assert(h[0] > 0);
  for (double x = low[0] + h[0] / 2; x < high[0]; x += h[0]) {
    double yinc = h[1];
    if (yinc == 0)
      yinc = h[0];
    for (double y = low[1] + h[1] / 2; y <= high[1]; y += yinc) {
      double zinc = h[2];
      if (zinc == 0)
        zinc = h[0];
      for (double z = low[2] + h[2] / 2; z <= high[2]; z += zinc) {
        if (generator_fun(Vect3d(x, y, z)) > max_generator_fun)
          max_generator_fun = generator_fun(Vect3d(x, y, z));
      }
    }
  }
  for (size_t i = 0; i < n; ++i) {
    Vect3d r =
        Vect3d(uni(generator), uni(generator), uni(generator)) * (high - low) +
        low;
    while (uni(generator) > generator_fun(r) / max_generator_fun) {
      r = Vect3d(uni(generator), uni(generator), uni(generator)) *
              (high - low) +
          low;
    }
    particles.push_back(r);
  }
}

template <typename T>
void radial_distribution_function(const T &particles, const double min,
                                  const double max, const int n,
                                  std::vector<double> &out) {

  // assume 3d for now...
  typedef Vector<double, 3> Vect3d;

  out.resize(n, 0);
  const double bsep = (max - min) / n;
  const double old_size = particles->get_lengthscale();
  particles->reset_neighbour_search(max);

  for (typename T::value_type &i : particles) {
    for (auto tpl : i.get_neighbours(particles)) {
      const Vect3d &dx = std::get<1>(tpl);
      const typename T::value_type &j = std::get<0>(tpl);
      if (i.get_id() == j.get_id())
        continue;
      const double r = dx.norm();
      const int index = std::floor((r - min) / bsep);
      if ((index >= 0) && (index < n))
        out[index] += 1.0;
    }
  }

  const double volume = (particles->get_high() - particles->get_low()).prod();
  const double rho = particles->size() / volume;

  const double PI = boost::math::constants::pi<double>();
  for (size_t i = 0; i < n; ++i) {
    out[i] /=
        particles->size() * (4.0 / 3.0) * PI *
        (std::pow((i + 1) * bsep + min, 3) - std::pow(i * bsep + min, 3)) * rho;
  }

  particles->reset_neighbour_search(old_size);
}

/// Assumes 2D particle set since we are drawing to a 2D canvas
template <typename Particles_t, typename Transform = IdentityTransform>
void draw_particles_with_search(
    std::string filename, const Particles_t &particles,
    const std::vector<Vector<double, 2>> &search_points, double search_radius,
    const Transform &transform = Transform()) {
#ifdef HAVE_CAIRO
  using position = typename Particles_t::position;
  const int image_size = 512;
  cairo_surface_t *surface =
      cairo_svg_surface_create(filename.c_str(), image_size, image_size);
  cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
  cairo_t *cr = cairo_create(surface);

  auto query = particles.get_query();
  vdouble2 min = particles.get_min();
  vdouble2 max = particles.get_max();
  auto dx = max - min;
  cairo_scale(cr, image_size / dx[0], image_size / dx[1]);
  cairo_translate(cr, -min[0], -min[1]);

  cairo_set_source_rgba(cr, 0.5, 0, 0, 0.5);
  vdouble2 lw(0.005, 0.005);
  // cairo_device_to_user_distance(cr, &lw[0], &lw[1]);

  // draw leaf buckets
  cairo_set_line_width(cr, lw[0]);
  for (auto i = query.get_subtree(); i != false; ++i) {
    auto ci = i.get_child_iterator();
    if (query.is_leaf_node(*ci)) {
      // draw its outline
      auto bounds = query.get_bounds(ci);
      vdouble2 rect[4];
      rect[0] = transform(vdouble2(bounds.bmin[0], bounds.bmin[1]));
      rect[1] = transform(vdouble2(bounds.bmax[0], bounds.bmin[1]));
      rect[2] = transform(vdouble2(bounds.bmax[0], bounds.bmax[1]));
      rect[3] = transform(vdouble2(bounds.bmin[0], bounds.bmax[1]));
      cairo_move_to(cr, rect[0][0], rect[0][1]);
      for (int i = 1; i < 4; ++i) {
        cairo_line_to(cr, rect[i][0], rect[i][1]);
      }
      cairo_close_path(cr);
      cairo_stroke(cr);
    }
  }

  // draw particles
  const double PI = boost::math::constants::pi<double>();
  for (size_t i = 0; i < particles.size(); ++i) {
    bool within_search = false;
    vdouble2 pos = transform(get<position>(particles)[i]);
    for (auto &search_point : search_points) {
      vdouble2 transformed_search_point = transform(search_point);

      if ((pos - transformed_search_point).squaredNorm() <=
          std::pow(search_radius, 2)) {
        within_search = true;
      }
    }
    if (within_search) {
      cairo_set_source_rgba(cr, 1.0, 0.0, 0.0, 0.5);
    } else {
      cairo_set_source_rgba(cr, 0, 0, 0, 0.5);
    }
    cairo_arc(cr, pos[0], pos[1], lw[0], 0, 2 * PI);
    cairo_fill(cr);
  }

  // draw search region and point
  for (auto &search_point : search_points) {
    cairo_set_source_rgba(cr, 0, 0, 0.5, 0.5);
    vdouble2 transformed_search_point = transform(search_point);
    cairo_arc(cr, transformed_search_point[0], transformed_search_point[1],
              search_radius, 0, 2 * PI);
    cairo_stroke(cr);
    cairo_arc(cr, transformed_search_point[0], transformed_search_point[1],
              lw[0], 0, 2 * PI);
    cairo_fill(cr);

    cairo_set_source_rgba(cr, 0, 0, 0.5, 0.2);
    for (auto i = query.template get_buckets_near_point<2>(
             search_point, search_radius, transform);
         i != false; ++i) {
      auto ci = i.get_child_iterator();
      auto bounds = query.get_bounds(ci);
      vdouble2 rect[4];
      rect[0] = transform(vdouble2(bounds.bmin[0], bounds.bmin[1]));
      rect[1] = transform(vdouble2(bounds.bmax[0], bounds.bmin[1]));
      rect[2] = transform(vdouble2(bounds.bmax[0], bounds.bmax[1]));
      rect[3] = transform(vdouble2(bounds.bmin[0], bounds.bmax[1]));

      // colour in search buckets
      cairo_move_to(cr, rect[0][0], rect[0][1]);
      for (int i = 1; i < 4; ++i) {
        cairo_line_to(cr, rect[i][0], rect[i][1]);
      }
      cairo_close_path(cr);
      cairo_fill(cr);
      for (auto j = query.get_bucket_particles(*ci); j != false; ++j) {
        vdouble2 pos = transform(get<position>(*j));
        cairo_arc(cr, pos[0], pos[1], lw[0], 0, 2 * PI);
        cairo_fill(cr);
      }
    }
    for (auto i =
             euclidean_search(query, search_point, search_radius, transform);
         i != false; ++i) {
      vdouble2 pos = transform(get<position>(*i));
      cairo_arc(cr, pos[0], pos[1], 2 * lw[0], 0, 2 * PI);
      cairo_fill(cr);
    }
  }

  cairo_destroy(cr);
  cairo_surface_destroy(surface);
#endif // HAVE_CAIRO
}

} // namespace Aboria

#endif /* UTILS_H_ */
