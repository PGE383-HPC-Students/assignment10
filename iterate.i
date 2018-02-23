%module iterate
%{
#define SWIG_FILE_WITH_INIT
#include "iterate.h"
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* u, int nx, int ny)};
%include "iterate.h"
