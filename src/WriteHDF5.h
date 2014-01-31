/*
 * WriteHDF5.h
 *
 *  Created on: 31 Jan 2014
 *      Author: mrobins
 */

#ifndef WRITEHDF5_H_
#define WRITEHDF5_H_

#include "Particles.h"
#include <hdf5>

template<typename ParticlesType>
class WriteParticle {
public:
	WriteParticle (hid_t dataset, hid_t datatype
      , hid_t dataspace, hid_t memspace)
    : m_dataset (dataset), m_datatype (datatype)
    , m_dataspace (dataspace), m_memspace (memspace)
    , m_pos () {}

  void operator ()(ParticlesType::Value const & v) {
    // Select the file position, 1 record at position 'pos'
    hsize_t count[] = { 1 } ;
    hsize_t offset[] = { m_pos++ } ;
    H5Sselect_hyperslab( m_dataspace
      , H5S_SELECT_SET
      , offset
      , NULL
      , count
      , NULL );

    const char * s = v.c_str ();
    H5Dwrite (m_dataset
      , m_datatype
      , m_memspace
      , m_dataspace
      , H5P_DEFAULT
      , &s );
    }

private:
  hid_t m_dataset;
  hid_t m_datatype;
  hid_t m_dataspace;
  hid_t m_memspace;
  int m_pos;
};

void writeVector (hid_t group, std::vector<std::string> const & v)
{
  hsize_t     dims[] = { m_files.size ()  } ;
  hid_t dataspace = H5Screate_simple(sizeof(dims)/sizeof(*dims)
                                    , dims, NULL);

  dims[0] = 1;
  hid_t memspace = H5Screate_simple(sizeof(dims)/sizeof(*dims)
                                    , dims, NULL);

  hid_t datatype = H5Tcopy (H5T_C_S1);
  H5Tset_size (datatype, H5T_VARIABLE);

  hid_t dataset = H5Dcreate1 (group, "files", datatype
                             , dataspace, H5P_DEFAULT);

  //
  // Select the "memory" to be written out - just 1 record.
  hsize_t offset[] = { 0 } ;
  hsize_t count[] = { 1 } ;
  H5Sselect_hyperslab( memspace, H5S_SELECT_SET, offset
                     , NULL, count, NULL );

  std::for_each (v.begin ()
      , v.end ()
      , WriteStrings (dataset, datatype, dataspace, memspace));

  H5Dclose (dataset);
  H5Sclose (dataspace);
  H5Sclose (memspace);
  H5Tclose (datatype);
}




#endif /* WRITEHDF5_H_ */
