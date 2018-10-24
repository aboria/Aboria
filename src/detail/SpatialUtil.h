
#ifndef DETAIL_SPATIAL_UTIL_H_
#define DETAIL_SPATIAL_UTIL_H_

#include "CudaInclude.h"
#include "Distance.h"
#include "Log.h"
#include "Vector.h"

#include <bitset>  // std::bitset
#include <iomanip> // std::setw
#include <limits>

namespace Aboria {
namespace detail {

CUDA_HOST_DEVICE
constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1) {
  return exp < 1
             ? result
             : ipow(base * base, exp / 2, (exp % 2) ? result * base : result);
}

template <typename T> CUDA_HOST_DEVICE constexpr T get_max() {
  return std::numeric_limits<T>::max();
}
#ifdef __CUDA_ARCH__
template <> CUDA_HOST_DEVICE constexpr unsigned int get_max<unsigned int>() {
  return NPP_MAX_16U;
}
template <> CUDA_HOST_DEVICE constexpr double get_max<double>() {
  return NPP_MAXABS_32F;
}
#endif

template <unsigned int D> struct bucket_index {
  typedef Vector<double, D> double_d;
  typedef Vector<unsigned int, D> unsigned_int_d;
  typedef Vector<int, D> int_d;

  unsigned_int_d m_size;

  CUDA_HOST_DEVICE
  bucket_index(){};

  CUDA_HOST_DEVICE
  bucket_index(const unsigned_int_d &size) : m_size(size) {}

  inline CUDA_HOST_DEVICE int collapse_index_vector(const int_d &vindex) const {
    int index = 0;
    unsigned int multiplier = 1.0;
    for (int i = D - 1; i >= 0; --i) {
      if (i != D - 1) {
        multiplier *= m_size[i + 1];
      }
      index += multiplier * vindex[i];
    }
    return index;
  }

  inline CUDA_HOST_DEVICE unsigned int
  collapse_index_vector(const unsigned_int_d &vindex) const {
    unsigned int index = 0;
    unsigned int multiplier = 1.0;
    for (int i = D - 1; i >= 0; --i) {
      if (i != D - 1) {
        multiplier *= m_size[i + 1];
      }
      index += multiplier * vindex[i];
    }
    return index;
  }

  inline CUDA_HOST_DEVICE int_d reassemble_index_vector(const int index) const {
    int_d vindex;
    int i = index;
    for (int d = D - 1; d >= 0; --d) {
      double div = (double)i / m_size[d];
      vindex[d] = std::round((div - std::floor(div)) * m_size[d]);
      i = std::floor(div);
    }
    return vindex;
  }

  inline CUDA_HOST_DEVICE unsigned_int_d
  reassemble_index_vector(const unsigned int index) const {
    unsigned_int_d vindex;
    unsigned int i = index;
    for (int d = D - 1; d >= 0; --d) {
      double div = (double)i / m_size[d];
      vindex[d] = std::round((div - std::floor(div)) * m_size[d]);
      i = std::floor(div);
    }
    return vindex;
  }
};

template <unsigned int D> struct point_to_bucket_index {
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  typedef Vector<unsigned int, D> unsigned_int_d;

  bucket_index<D> m_bucket_index;
  double_d m_bucket_side_length;
  double_d m_inv_bucket_side_length;
  bbox<D> m_bounds;

  CUDA_HOST_DEVICE
  point_to_bucket_index(){};

  CUDA_HOST_DEVICE
  point_to_bucket_index(const unsigned_int_d &size,
                        const double_d &bucket_side_length,
                        const bbox<D> &bounds)
      : m_bucket_index(size), m_bucket_side_length(bucket_side_length),
        m_inv_bucket_side_length(1.0 / bucket_side_length), m_bounds(bounds) {}

  CUDA_HOST_DEVICE
  int operator()(const double_d &v) const { return find_bucket_index(v); }

  inline CUDA_HOST_DEVICE int collapse_index_vector(const int_d &vindex) const {
    return m_bucket_index.collapse_index_vector(vindex);
  }

  inline CUDA_HOST_DEVICE int_d
  find_bucket_index_vector(const double_d &r) const {
    // find the raster indices of p's bucket
    // std::cout << "r = "<<r<<" indexv =
    // "<<floor((r-m_bounds.bmin)/m_bucket_side_length)<<std::endl;
    return floor((r - m_bounds.bmin) * m_inv_bucket_side_length)
        .template cast<int>();
  }

  // hash a point in the unit square to the index of
  // the grid bucket that contains it
  inline CUDA_HOST_DEVICE int find_bucket_index(const double_d &r) const {
    return m_bucket_index.collapse_index_vector(find_bucket_index_vector(r));
  }

  inline CUDA_HOST_DEVICE int find_max_bucket_point(const int_d &vindex) const {
    return (vindex + 1) * m_bucket_side_length + m_bounds.bmin;
  }

  inline CUDA_HOST_DEVICE double_d
  find_bucket_centre(const int_d &vindex) const {
    return (vindex + 0.5) * m_bucket_side_length + m_bounds.bmin;
  }

  CUDA_HOST_DEVICE
  int get_min_index_by_quadrant(const double r, const int i,
                                const bool up) const {
    // std::cout << "up = "<<up<<"r = "<<r<<" i = "<<i << std::endl;
    return std::floor(
        (r + (up ? 0.5 : -0.5) * m_bucket_side_length[i] - m_bounds.bmin[i]) *
        m_inv_bucket_side_length[i]);
  }

  /*
  CUDA_HOST_DEVICE
  double_d get_dist_to_bucket_unsigned(const double_d &r,
                                       const int_d &target_index) const {

    double_d dx =
        (target_index + 0.5) * m_bucket_side_length + m_bounds.bmin - r;
    for (int i = 0; i < D; ++i) {
      dx[i] = std::max(std::abs(dx[i]) - 0.5 * m_bucket_side_length[i], 0.0);
    }
    return dx;
  }

  CUDA_HOST_DEVICE
  double_d get_dist_to_bucket_signed(const double_d &r,
                                     const int_d &target_index) const {

    double_d dx =
        (target_index + 0.5) * m_bucket_side_length + m_bounds.bmin - r;
    for (int i = 0; i < D; ++i) {
      dx[i] = std::copysign(
          std::max(std::abs(dx[i]) - 0.5 * m_bucket_side_length[i], 0.0),
          dx[i]);
    }
    return dx;
  }

  CUDA_HOST_DEVICE
  double_d get_dist_to_bucket_vertex(const double_d &r,
                                     const int_d &target_index,
                                     const int vertex) const {
    double_d dx;
    for (int i = 0; i < D; ++i) {
      const bool jth_bit = (1 == ((vertex >> i) & 1));
      dx[i] = (target_index[i] + jth_bit) * m_bucket_side_length[i] +
              m_bounds.bmin[i] - r[i];
    }
    return dx;
  }

  CUDA_HOST_DEVICE
  double get_dist_to_bucket(const double r, const int my_index,
                            const int target_index, const int i) const {
    if (my_index < target_index) {
      // compare point to lower edge of bucket, return a positive distance
      return target_index * m_bucket_side_length[i] + m_bounds.bmin[i] - r;
    } else if (my_index > target_index) {
      // compare point to upper edge of bucket, return a negative distance
      return ((target_index + 1) * m_bucket_side_length[i] + m_bounds.bmin[i]) -
             r;
    } else
      // same index, return 0.0
      return 0.0;
  }

  CUDA_HOST_DEVICE
  double get_dist_by_quadrant(const double r, const int index, const int i,
                              const bool up) const {
    if (up) {
      // compare point to lower edge of bucket, return a positive distance
      return index * m_bucket_side_length[i] + m_bounds.bmin[i] - r;
    } else {
      // compare point to upper edge of bucket, return a positive distance
      return r - ((index + 1) * m_bucket_side_length[i] + m_bounds.bmin[i]);
    }
  }
  */
};

// Utility functions to encode leaves and children in single int
inline CUDA_HOST_DEVICE bool is_empty(int id) { return id == (int)0xffffffff; }

inline CUDA_HOST_DEVICE bool is_node(int id) { return id > 0; }

inline CUDA_HOST_DEVICE bool is_leaf(int id) { return id < 0; }

inline CUDA_HOST_DEVICE constexpr int get_empty_id() { return 0xffffffff; }

inline CUDA_HOST_DEVICE int get_leaf_id(int offset) {
  return 0x80000000 | offset;
}

inline CUDA_HOST_DEVICE int get_leaf_offset(int id) { return 0x80000000 ^ id; }

inline CUDA_HOST_DEVICE int child_tag_mask(int tag, int which_child, int level,
                                           int max_level, unsigned int D) {
  int shift = (max_level - level) * D;
  return tag | (which_child << shift);
}

template <int CODE> struct is_a {
  typedef int result_type;
  inline CUDA_HOST_DEVICE int operator()(int code) {
    return code == CODE ? 1 : 0;
  }
};

template <unsigned int D>
CUDA_HOST_DEVICE int point_to_tag(const Vector<double, D> &p, bbox<D> box,
                                  int max_level) {
  typedef Vector<double, D> double_d;
  typedef Vector<int, D> int_d;
  int result = 0;

  for (int level = 1; level <= max_level; ++level) {
    double_d mid;
    int_d hi_half;

    for (size_t i = 0; i < D; ++i) {
      // Classify in i-direction
      mid[i] = 0.5f * (box.bmin[i] + box.bmax[i]);
      hi_half[i] = (p[i] < mid[i]) ? 0 : 1;

      // Push the bit into the result as we build it
      result |= hi_half[i];
      result <<= 1;
    }

    // Shrink the bounding box, still encapsulating the point
    for (size_t i = 0; i < D; ++i) {
      box.bmin[i] = (hi_half[i]) ? mid[i] : box.bmin[i];
      box.bmax[i] = (hi_half[i]) ? box.bmax[i] : mid[i];
    }
  }
  // Unshift the last
  result >>= 1;

  return result;
}

template <unsigned int D> void print_tag(int tag, int max_level) {
  for (int level = 1; level <= max_level; ++level) {
    std::bitset<D> bits = tag >> (max_level - level) * D;
    std::cout << bits << " ";
  }
}

template <unsigned int D, typename Vector>
void print_active_nodes(const Vector &active_nodes, int max_level) {
  std::cout << "Active nodes:\n      ";
  for (int i = 1; i <= max_level; ++i) {
    for (int j = 0; j < D; j++) {
      if (j == 0) {
        std::cout << "x";
      } else if (j == 1) {
        std::cout << "y";
      } else if (j == 2) {
        std::cout << "z";
      }
    }
    std::cout << " ";
  }
  std::cout << std::endl;
  for (size_t i = 0; i < active_nodes.size(); ++i) {
    std::cout << std::setw(4) << i << ": ";
    print_tag<D>(active_nodes[i], max_level);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <unsigned int D, typename Vector>
void print_children(const Vector &children, int max_level) {
  std::cout << "Children:\n      ";
  for (int i = 1; i <= max_level; ++i) {
    for (int j = 0; j < D; j++) {
      if (j == 0) {
        std::cout << "x";
      } else if (j == 1) {
        std::cout << "y";
      } else if (j == 2) {
        std::cout << "z";
      }
    }
    std::cout << " ";
  }
  std::cout << std::endl;
  for (size_t i = 0; i < children.size(); ++i) {
    std::cout << std::setw(4) << i << ": ";
    print_tag<D>(children[i], max_level);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <typename Vector>
void print_child_bounds(const Vector &lower_bounds,
                        const Vector &upper_bounds) {
  std::cout << "Child bounds:\n      [ lower upper count ]\n";
  for (size_t i = 0; i < lower_bounds.size(); ++i) {
    std::cout << std::setw(4) << i << ": [ ";
    std::cout << std::setw(4) << lower_bounds[i] << "  ";
    std::cout << std::setw(4) << upper_bounds[i] << "  ";
    std::cout << std::setw(4) << upper_bounds[i] - lower_bounds[i] << "  ]";
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// Markers
enum { NODE = 1, LEAF = 2, EMPTY = 4 };

template <typename Vector>
void print_child_node_kind(const Vector &child_node_kind) {
  std::cout << "child_node_kind:\n";
  for (size_t i = 0; i < child_node_kind.size(); ++i) {
    std::cout << std::setw(4) << i << ": [ ";
    std::cout << std::setw(5) << std::right;
    switch (child_node_kind[i]) {
    case EMPTY:
      std::cout << "EMPTY ]";
      break;
    case LEAF:
      std::cout << "LEAF ]";
      break;
    case NODE:
      std::cout << "NODE ]";
      break;
    default:
      std::cout << "ERROR ]";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <typename Vector>
void print_child_enumeration(const Vector &child_node_kind,
                             const Vector &nodes_on_this_level,
                             const Vector &leaves_on_this_level) {
  std::cout << "Node/leaf enumeration:\n      [ nodeid leafid ]\n";
  for (size_t i = 0; i < child_node_kind.size(); ++i) {
    std::cout << std::setw(4) << i << ": [ ";
    switch (child_node_kind[i]) {
    case EMPTY:
      std::cout << std::setw(4) << "."
                << "   " << std::setw(4) << "."
                << "   ]";
      break;
    case LEAF:
      std::cout << std::setw(4) << "."
                << "   " << std::setw(4) << leaves_on_this_level[i] << "   ]";
      break;
    case NODE:
      std::cout << std::setw(4) << nodes_on_this_level[i] << "   "
                << std::setw(4) << "."
                << "   ]";
      break;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <typename Vector>
void print_nodes(const Vector &nodes, const unsigned int D) {
  std::cout << "Octtree nodes:\n";
  std::cout << "          [ nodeid  leafid ]\n";

  size_t next_level = 0;
  size_t children_at_next_level = 1 << D;

  for (size_t i = 0; i < nodes.size(); ++i) {
    if (i == next_level) {
      std::cout << "          [================]\n";
      next_level += children_at_next_level;
    } else if (i % children_at_next_level == 0) {
      std::cout << "          [----------------]\n";
    }

    if (is_empty(nodes[i])) {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << "."
                << "    ";
      std::cout << std::setw(4) << "."
                << "   ]\n";
    } else if (is_leaf(nodes[i])) {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << "."
                << "    ";
      std::cout << std::setw(4) << get_leaf_offset(nodes[i]) << "   ]\n";
    } else {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << nodes[i] << "    ";
      std::cout << std::setw(4) << "."
                << "   ]\n";
    }
  }
  std::cout << "          [================]\n";
}

template <typename Vector> void print_leaves(const Vector &leaves) {
  std::cout << "Octtree leaves:\n";
  std::cout << "          [ lower    upper ]\n";

  for (size_t i = 0; i < leaves.size(); ++i) {
    std::cout << std::setw(7) << i << " : [ ";
    std::cout << std::setw(4) << static_cast<vint2>(leaves[i])[0] << "    ";
    std::cout << std::setw(4) << static_cast<vint2>(leaves[i])[1] << "   ]\n";
  }
}

} // namespace detail
} // namespace Aboria
#endif // SPATIAL_UTIL_H_
