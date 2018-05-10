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

#ifndef DETAIL_H2_LIB_H_
#define DETAIL_H2_LIB_H_

#include "detail/FastMultipoleMethod.h"

#ifdef HAVE_H2LIB

//#include <GL/glut.h>
extern "C" {
#include <h2arith.h>
#include <h2matrix.h>
#include <h2update.h>
#include <harith.h>
#include <hcoarsen.h>
#undef I
}

namespace Aboria {

namespace detail {

bool admissible_max_cluster(pcluster rc, pcluster cc, void *data) {

  double eta = *static_cast<double *>(data);

  if (eta < 0) {
    double vol = 1.0;
    for (size_t i = 0; i < rc->dim; i++) {
        const double side = std::max(rc->bmax[i] - cc->bmin[i],
                                     cc->bmax[i] - rc->bmin[i]);
        vol *= side > 0 ? side : 0;
    }
    if (vol > std::numeric_limits<double>::epsilon()) {
        return false;
    } else {
        return true;
    }
  }

  const double diamt = getdiam_max_cluster(rc);
  const double diams = getdiam_max_cluster(cc);
  const double dist = getdist_max_cluster(rc, cc);

  const double a = REAL_MAX(diamt, diams);

  bool i;
  if (a <= eta * dist) {
    i = true;
  } else {
    i = false;
  }

  return i;
}

template <typename Particles, typename Expansions>
static void assemble_h2matrix_row_clusterbasis(pcclusterbasis rbc, uint rname,
                                               void *data) {
  auto &data_cast = *static_cast<std::tuple<Expansions *, Particles *> *>(data);

  const Expansions &expansions = *std::get<0>(data_cast);
  const Particles &particles = *std::get<1>(data_cast);
  pclusterbasis rb = (pclusterbasis)rbc;

  if (rb->sons > 0) {
    // find orders of sons
    std::vector<size_t> orders(rb->sons);
    size_t max_order = 0;
    for (size_t i = 0; i < rb->sons; ++i) {
      orders[i] = static_cast<size_t>(std::round(std::pow(
          rb->son[i]->k / expansions.block_rows, 1.0 / Expansions::dimension)));
      if (orders[i] > max_order)
        max_order = orders[i];
    }
    // increment max order by beta
    const size_t order = max_order + expansions.m_beta;
    const uint k =
        std::pow(order, Expansions::dimension) * expansions.block_rows;
    resize_clusterbasis(rb, k);
    for (size_t i = 0; i < rb->sons; ++i) {
      expansions.L2L_amatrix(&rb->son[i]->E, rb->son[i]->t, rb->t, orders[i],
                             order);
    }
  } else {
    const uint k = std::pow(expansions.m_order, Expansions::dimension) *
                   expansions.block_rows;
    resize_amatrix(&rb->V, rb->t->size, k);
    expansions.L2P_amatrix(&rb->V, rb->t, rb->t->idx, rb->t->size, particles);
    rb->k = k;
    update_clusterbasis(rb);
  }
}

template <typename Particles, typename Expansions>
static void assemble_h2matrix_col_clusterbasis(pcclusterbasis rbc, uint rname,
                                               void *data) {
  auto &data_cast = *static_cast<std::tuple<Expansions *, Particles *> *>(data);

  const Expansions &expansions = *std::get<0>(data_cast);
  const Particles &particles = *std::get<1>(data_cast);
  pclusterbasis rb = (pclusterbasis)rbc;

  if (rb->sons > 0) {
    // find orders of sons
    std::vector<size_t> orders(rb->sons);
    size_t max_order = 0;
    for (size_t i = 0; i < rb->sons; ++i) {
      orders[i] = static_cast<size_t>(std::round(std::pow(
          rb->son[i]->k / expansions.block_cols, 1.0 / Expansions::dimension)));
      if (orders[i] > max_order)
        max_order = orders[i];
    }
    // increment max order by beta
    const size_t order = max_order + expansions.m_beta;
    const uint k =
        std::pow(order, Expansions::dimension) * expansions.block_cols;

    resize_clusterbasis(rb, k);
    for (size_t i = 0; i < rb->sons; ++i) {
      // should this be transposed?
      expansions.M2M_trans_amatrix(&rb->son[i]->E, rb->t, rb->son[i]->t, order,
                                   orders[i]);
    }
  } else {
    const uint k = std::pow(expansions.m_order, Expansions::dimension) *
                   expansions.block_cols;
    resize_amatrix(&rb->V, rb->t->size, k);
    expansions.P2M_trans_amatrix(&rb->V, rb->t, rb->t->idx, rb->t->size,
                                 particles);
    rb->k = k;
    update_clusterbasis(rb);
  }
}

template <typename RowParticles, typename ColParticles, typename Expansions,
          typename Kernel>
static void assemble_block_h2matrix(pcblock b, uint bname, uint rname,
                                    uint cname, uint pardepth, void *data) {
  auto &data_cast =
      *static_cast<std::tuple<Expansions *, Kernel *, RowParticles *,
                              ColParticles *, ph2matrix *> *>(data);

  const Expansions &expansions = *std::get<0>(data_cast);
  const Kernel &kernel = *std::get<1>(data_cast);
  const RowParticles &row_particles = *std::get<2>(data_cast);
  const ColParticles &col_particles = *std::get<3>(data_cast);
  ph2matrix *enum_h2 = std::get<4>(data_cast);
  ph2matrix h2 = enum_h2[bname];

  if (h2->u) {
    const uint kr = h2->u->rb->k;
    const uint kc = h2->u->cb->k;
    const size_t orderr = static_cast<size_t>(std::round(
        std::pow(kr / expansions.block_rows, 1.0 / Expansions::dimension)));
    const size_t orderc = static_cast<size_t>(std::round(
        std::pow(kc / expansions.block_cols, 1.0 / Expansions::dimension)));
    resize_amatrix(&h2->u->S, kr, kc);
    expansions.M2L_amatrix(&h2->u->S, h2->u->rb->t, h2->u->cb->t, orderr,
                           orderc);

  } else if (h2->f) {
    detail::P2P_amatrix(h2->f, h2->rb->t->idx, h2->rb->t->size, h2->cb->t->idx,
                        h2->cb->t->size, row_particles, col_particles, kernel);
  }
}

} // namespace detail

} // namespace Aboria

#endif // HAVE_H2LIB

#endif
