//---------------------------------------------------------------------
// IO_.CXX: a series of I/O functions, mostly convenient wrappers for
//          existing external routines.
//
// Author: Vlad Popovic <popovici@recetox.muni.cz>
//---------------------------------------------------------------------
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <boost/python.hpp>
#include <numpy/ndarrayobject.h>
#include <string>
#include <iostream>


#include <openslide.h>

namespace bp = boost::python;

template <typename T>
void reference_contiguous_array(PyObject* in, PyArrayObject* &in_con, T* &ptr, long unsigned &count)
{
    in_con = PyArray_GETCONTIGUOUS((PyArrayObject*)in);
    ptr = (T*)PyArray_DATA(in_con);
    int num_dim = PyArray_NDIM(in_con);
    npy_intp* pdim = PyArray_DIMS(in_con);
    count = 1;
    for (int i = 0; i < num_dim; i++) {
        count *= static_cast<long unsigned>(pdim[i]);
    }
}


void dereference(PyObject* o)
{
    Py_DECREF(o);
}


PyObject* entry_square_matrix(PyObject* input_matrix)
{
    // get the input array
    double* ptr;
    long unsigned count;

    PyArrayObject* input_contigous_array;
    reference_contiguous_array(input_matrix, input_contigous_array, ptr, count);

    // create the output array
    npy_intp dst_dim[1];
    dst_dim[0] = count;
    PyObject* out_matrix = PyArray_SimpleNew(1, dst_dim, NPY_FLOAT64);
    double* ptr_out;
    PyArrayObject* output_contigous_array;
    reference_contiguous_array(out_matrix, output_contigous_array, ptr_out, count);
    for (int i = 0; i < count; i++) {
        ptr_out[i] = ptr[i] * ptr[i];
    }
    dereference((PyObject*)input_contigous_array);
    dereference((PyObject*)output_contigous_array);
    return out_matrix;
}

// OSL_READ_REGION
// Use OpenSlide's read_region function to read a rectangular region of a layer in a multiresolution image.
// The function is meant to be called from Python with a numpu.ndarray as destination of the read region.
// Every call causes the image file to be opened and its header parsed, therefore the overhead might be
// important in the case of small regions. However, for large regions it is a better choice than the
// wrapper provided by OpenSlide library, since it does not use Python Imaging Library under the hood
// and, thus, is not limited to image sizes that would fit a 32bit integer. The function does not allocate
// the memory for the image, but expects a pointer to a memory destination. The required size for the
// buffer is 4 x width x height bytes.
//
// Args:
//  filename (string)
//  dst (PyObject): a pointer to a numpy.ndarray PRE-ALLOCATED
//  x, y (long unsigned): the coordinates of the top-left corner of the region, in layer-0.
//  width, height (long unsigned)
//  level (unsigned)
//
// Returns:
//  0: success
// -1: cannot access buffer
// -2: cannot open file
// -3: region coordinates or size out of boundaries
// -4: buffer size mismatch

int osl_read_region(const std::string& filename,
                    PyObject* dst,   // destination of the read region
                    long unsigned x, long unsigned y,  // level 0 coordinates of the top-left corner
                    long unsigned width, long unsigned height,  // width and height of the region to read
                    unsigned level)
{
    unsigned int* buf = 0;
    PyArrayObject* dst_buf = 0;
    long unsigned buf_size = 0;

    reference_contiguous_array(dst, dst_buf, buf, buf_size);
    if (!buf || !dst_buf)
        // cannot access buffer
        return -1;

    openslide_t* osl_reader = openslide_open(filename.c_str());
    if (!osl_reader)
        // cannot open file
        return -2;
    long int img_w, img_h;
    long int w, h;
    int nlevels = openslide_get_level_count(osl_reader);

    openslide_get_level0_dimensions(osl_reader, &w, &h);
    openslide_get_level_dimensions(osl_reader, level, &img_w, &img_h);
    if (x > w || y > h || width > img_w || height > img_h)
        // region specification error
        return -3;

    if (buf_size != 4*width*height)
        // buffer size mismatch
        return -4;

    openslide_read_region(osl_reader, buf, x, y, level, width, height);

    openslide_close(osl_reader);
    dereference((PyObject*)dst_buf);

    return 0;
}


BOOST_PYTHON_MODULE(io_)
{
    import_array();
//    bp::def("square_matrix", entry_square_matrix); // for experiments
    bp::def("osl_read_region_", osl_read_region);
}
