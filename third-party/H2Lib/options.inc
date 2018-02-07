
# --- Build options configuration

# Brief output during build process
BRIEF_OUTPUT=1

# Build for debugging
#DEBUG=1

# Build optimized
OPT=1

# Use complex numbers
#USE_COMPLEX=1

# Use BLAS/LAPACK for linear algebra
USE_BLAS=1

# Use Cairo for visualization
#USE_CAIRO=0

# Use OpenMP for parallelization
USE_OPENMP=1

# Use OpenCL for heterogeneous computing
#USE_OPENCL=1

# Use NetCDF for writing and reading numerical data
#USE_NETCDF=1

# Use zlib for compressing data files
#USE_ZLIB=1

# Use GTK+ 3.x for visualization
#USE_GTK3=1

# Use FreeGLUT for visualization
#USE_FREEGLUT=1

# Do nothing if a rank-zero matrix is added to another matrix
# (i.e., do not recompress the target matrix)
HARITH_RKMATRIX_QUICK_EXIT=1

# Do nothing if a zero matrix is added to another matrix
# (i.e., do not recompress the target matrix)
HARITH_AMATRIX_QUICK_EXIT=1
