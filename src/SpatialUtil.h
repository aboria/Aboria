
#ifndef SPATIAL_UTIL_H_ 
#define SPATIAL_UTIL_H_ 

#include "Vector.h"
#include "CudaInclude.h"

#include <bitset>         // std::bitset
#include <iomanip>      // std::setw

namespace Aboria {

#define FLT_MAX 1.0


template <typename T>
struct plus {
    T operator()(const T& t1, const T& t2) { return t1+t2; }
};

template<unsigned int D>
struct bbox {
    typedef Vector<double,D> double_d;
    double_d bmax;
    double_d bmin;

    inline CUDA_HOST_DEVICE
    bbox() : bmin(double_d(FLT_MAX)), bmax(double_d(-FLT_MAX))
    {}

    inline CUDA_HOST_DEVICE
    bbox(const double_d &p) : bmin(p),bmax(p)
    {}

    inline CUDA_HOST_DEVICE
    bbox operator+(const bbox &arg) {
        bbox bounds;
        for (int i=0; i<D; i++) {
            bounds.bmin[i] = std::min(bmin[i], arg.bmin[i]);
            bounds.bmax[i] = std::max(bmax[i], arg.bmax[i]);
        }
        return bounds;
    }

    inline CUDA_HOST_DEVICE
    bool operator<(const bbox &arg) {
        bbox bounds;
        bool within = true;
        for (int i=0; i<D; i++) {
            within |= bmin[i] >= arg.bmin[i];
            within |= bmax[i] < arg.bmax[i];
        }
        return within;
    }

    inline CUDA_HOST_DEVICE
    bool operator<=(const bbox &arg) {
        bbox bounds;
        bool within = true;
        for (int i=0; i<D; i++) {
            within |= bmin[i] >= arg.bmin[i];
            within |= bmax[i] <= arg.bmax[i];
        }
        return within;
    }

    inline CUDA_HOST_DEVICE
    bool is_empty() {
        for (int i=0; i<D; i++) {
            if (bmax[i] < bmin[i]) return true;
        }
        return false;
    }
 
};


// Utility functions to encode leaves and children in single int
inline CUDA_HOST_DEVICE
bool is_empty(int id) { return id == 0xffffffff; }

inline CUDA_HOST_DEVICE
bool is_node(int id) { return id > 0; }

inline CUDA_HOST_DEVICE
bool is_leaf(int id) { return id < 0; }

inline CUDA_HOST_DEVICE
int get_empty_id() { return 0xffffffff; }

inline CUDA_HOST_DEVICE
int get_leaf_id(int offset) { return 0x80000000 | offset; }

inline CUDA_HOST_DEVICE
int get_leaf_offset(int id) { return 0x80000000 ^ id; }

inline CUDA_HOST_DEVICE
int child_tag_mask(int tag, int which_child, int level, int max_level)
{
  int shift = (max_level - level) * 2;
  return tag | (which_child << shift);
}

template <int CODE>
struct is_a
{
  typedef int result_type;
  inline CUDA_HOST_DEVICE
  int operator()(int code) { return code == CODE ? 1 : 0; }
};





CUDA_HOST_DEVICE
template<unsigned int D>
int point_to_tag(const Vector<double,D> &p, bbox<D> box, int max_level) {
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;
    int result = 0;
  
    for (int level = 1 ; level <= max_level ; ++level) {
        double_d mid;
        int_d hi_half;
    
        for (int i=0; i<D; i++) {
            // Classify in i-direction
            mid[i] = 0.5f * (box.bmin[i] + box.bmax[i]);
            hi_half[i] = (p[i] < mid[i]) ? 0 : 1;
  
            // Push the bit into the result as we build it
            result |= hi_half[i];
            result <<= 1;
        }
  
        // Shrink the bounding box, still encapsulating the point
        for (int i=0; i<D; i++) {
            box.bmin[i] = (hi_half[i]) ? mid[i] : box.bmin[i];
            box.bmax[i] = (hi_half[i]) ? box.bmax[i] : mid[i];
        }

  }
  // Unshift the last
  result >>= 1;

  return result;
}

void print_tag(int tag, int max_level)
{
  for (int level = 1 ; level <= max_level ; ++level)
  {
    std::bitset<2> bits = tag >> (max_level - level) * 2;
    std::cout << bits << " ";
  }
}

template<typename Vector>
void print_active_nodes(const Vector &active_nodes, int max_level)
{
  std::cout << "Active nodes:\n      ";
  for (int i = 1 ; i <= max_level ; ++i)
  {
    std::cout << "xy ";
  }
  std::cout << std::endl;
  for (int i = 0 ; i < active_nodes.size() ; ++i)
  {
    std::cout << std::setw(4) << i << ": ";
    print_tag(active_nodes[i], max_level);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<typename Vector>
void print_children(const Vector &children, int max_level)
{
  std::cout << "Children:\n      ";
  for (int i = 1 ; i <= max_level ; ++i)
  {
    std::cout << "xy ";
  }
  std::cout << std::endl;
  for (int i = 0 ; i < children.size() ; ++i)
  {
    std::cout << std::setw(4) << i << ": ";
    print_tag(children[i], max_level);
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<typename Vector>
void print_child_bounds(const Vector &lower_bounds,
                        const Vector &upper_bounds)
{
  std::cout << "Child bounds:\n      [ lower upper count ]\n";
  for (int i = 0 ; i < lower_bounds.size() ; ++i)
  {
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

template<typename Vector>
void print_child_node_kind(const Vector &child_node_kind)
{
  std::cout << "child_node_kind:\n";
  for (int i = 0 ; i < child_node_kind.size() ; ++i)
  {
    std::cout << std::setw(4) << i << ": [ ";
    std::cout << std::setw(5) << std::right;
    switch (child_node_kind[i])
    {
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

template<typename Vector>
void print_child_enumeration(const Vector &child_node_kind,
                             const Vector &nodes_on_this_level,
                             const Vector &leaves_on_this_level)
{
  std::cout << "Node/leaf enumeration:\n      [ nodeid leafid ]\n";
  for (int i = 0 ; i < child_node_kind.size() ; ++i)
  {
    std::cout << std::setw(4) << i << ": [ ";
    switch (child_node_kind[i])
    {
    case EMPTY:
      std::cout << std::setw(4) << "." << "   " << std::setw(4) << "." << "   ]";
      break;
    case LEAF:
      std::cout << std::setw(4) << "." << "   " << std::setw(4) << leaves_on_this_level[i] << "   ]";
      break;
    case NODE:
      std::cout << std::setw(4) << nodes_on_this_level[i] << "   " << std::setw(4) << "." << "   ]";
      break;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<typename Vector>
void print_nodes(const Vector &nodes)
{
  std::cout << "Quadtree nodes:\n";
  std::cout << "          [ nodeid  leafid ]\n";
  
  int next_level = 0;
  int children_at_next_level = 4;

  for (int i = 0 ; i < nodes.size() ; ++i)
  {
    if (i == next_level)
    {
      std::cout << "          [================]\n";
      next_level += children_at_next_level;
      children_at_next_level = 0;
    }
    else if (i % 4 == 0)
    {
      std::cout << "          [----------------]\n";
    }

    if (is_empty(nodes[i]))
    {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << "." << "    ";
      std::cout << std::setw(4) << "." << "   ]\n";
    }
    else if (is_leaf(nodes[i]))
    {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << "." << "    ";
      std::cout << std::setw(4) << get_leaf_offset(nodes[i]) << "   ]\n";
    }
    else
    {
      std::cout << std::setw(7) << i << " : [ ";
      std::cout << std::setw(4) << nodes[i] << "    ";
      std::cout << std::setw(4) << "." << "   ]\n";
    }
  }
  std::cout << "          [================]\n";
}

template<typename Vector>
void print_leaves(const Vector &leaves)
{
  std::cout << "Quadtree leaves:\n";
  std::cout << "          [ lower    upper ]\n";
  
  for (int i = 0 ; i < leaves.size() ; ++i)
  {
    std::cout << std::setw(7) << i << " : [ ";
    std::cout << std::setw(4) << leaves[i].x << "    ";
    std::cout << std::setw(4) << leaves[i].y << "   ]\n";
  }
}

}
#endif //SPATIAL_UTIL_H_ 
