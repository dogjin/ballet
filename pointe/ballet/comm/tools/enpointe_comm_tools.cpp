#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <ballet/comm/bi2de.h>
#include <ballet/comm/de2bi.h>
#include <iostream>

using namespace std;

static PyObject *
ballet_bi2de(PyObject *self, PyObject *args, PyObject *kwds)
{

    PyObject * arg1;
    PyObject * arg2;

    static char * kwlist[] = {"b","flg",NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds,"OO",kwlist,&arg1,
        &arg2))
    {
        return NULL;
    }

    if (!PyArray_Check(arg1))
    {
        PyErr_Format(PyExc_TypeError,
            "Input must be a numpy.ndarray type, not %.500s",
            Py_TYPE(arg1)->tp_name);
        return 0;
    }

    // cast to PyArrayObject type
    PyArrayObject * b = (PyArrayObject*)arg1;

    // get pointer to the data buffer
    char * b_data = PyArray_BYTES(b);

    // get the number of array dimensions
    npy_intp size = PyArray_Size((PyObject*)b);
    npy_intp * shape = PyArray_SHAPE(b);
    npy_intp dim1 = shape[0];
    npy_intp dim2 = shape[1];
 
    // create armadillo matrix
    arma::imat b1(dim1,dim2);

    for (npy_intp r=0; r<dim1; ++r)
    {

        for (npy_intp c=0; c<dim2; ++c)
        {
            char * tmp = b_data + r*PyArray_STRIDE(b,0) + c*PyArray_STRIDE(b,1);
            b1(r,c) = *((npy_int64*)tmp);
        }

    }

    // run c++ bi2de
    arma::imat d1 = ballet::bi2de(b1);

    // create return array
    int d2_ndim = 2;
    npy_intp d2_dims[d2_ndim];
    d2_dims[0] = d1.n_rows;
    d2_dims[1] = d1.n_cols;
    PyObject * d2 = (PyObject*)PyArray_NewFromDescr(&PyArray_Type,
        PyArray_DescrFromType(NPY_INT64),
        d2_ndim,
        d2_dims,
        NULL,
        NULL,
        NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE,
        NULL);

    // get pointer to array data buffer
    char * d2_data = PyArray_BYTES((PyArrayObject*)d2);
    npy_intp d2s0 = PyArray_STRIDE((PyArrayObject*)d2,0);
    npy_intp d2s1 = PyArray_STRIDE((PyArrayObject*)d2,1);

    // copy result from armadillo c++ to numpy array
    for (npy_intp r=0; r<d1.n_rows; ++r)
    {
        for (npy_intp c=0; c<d1.n_cols; ++c)
        {
            char * tmp = d2_data + r*d2s0 + c*d2s1;
            (*((npy_int64*)tmp)) = d1(r,c);
        }
    }

    Py_INCREF(d2);
    return d2;

}

static PyObject *
ballet_de2bi(PyObject *self, PyObject *args, PyObject *kwds)
{

    PyObject * arg001;
    PyObject * arg002;
    PyObject * arg003;

    static char * kwlist[] = {"d","n","flg",NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds,"OOO",kwlist,&arg001,&arg002,
        &arg003))
    {
        return NULL;
    }

    if (!PyArray_Check(arg001))
    {
        PyErr_Format(PyExc_TypeError,
            "Input must be a numpy.ndarray type, not %.500s",
            Py_TYPE(arg001)->tp_name);
        return 0;
    }

    long n;
    if (PyLong_Check(arg002) || PyInt_Check(arg002))
    {

        n = PyLong_AsLongLong(arg002);

    }
    else
    {
        PyErr_SetString(PyExc_Exception,
            "Specified number of columns must be integer type.");
        return 0;
    }

    // cast to PyArrayObject type
    PyArrayObject * d = (PyArrayObject*)arg001;

    // get pointer to the data buffer
    char * d_data = PyArray_BYTES(d);

    // get the number of array dimensions
    npy_intp * shape = PyArray_SHAPE(d);
    npy_intp dim1 = shape[0];
    npy_intp dim2 = shape[1];
 
    // create armadillo matrix
    arma::imat d1(dim1,dim2);

    for (npy_intp r=0; r<dim1; ++r)
    {

        for (npy_intp c=0; c<dim2; ++c)
        {
            char * tmp = d_data + r*PyArray_STRIDE(d,0) + c*PyArray_STRIDE(d,1);
            d1(r,c) = *((npy_int64*)tmp);
        }

    }

    // run c++ bi2de
    arma::imat b = ballet::de2bi(d1,n);

    // create return array
    int b1_ndim = 2;
    npy_intp b1_dims[b1_ndim];
    b1_dims[0] = b.n_rows;
    b1_dims[1] = b.n_cols;
    PyObject * b1 = (PyObject*)PyArray_NewFromDescr(&PyArray_Type,
        PyArray_DescrFromType(NPY_INT64),
        b1_ndim,
        b1_dims,
        NULL,
        NULL,
        NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE,
        NULL);

    // get pointer to array data buffer
    char * b1_data = PyArray_BYTES((PyArrayObject*)b1);
    npy_intp b1s0 = PyArray_STRIDE((PyArrayObject*)b1,0);
    npy_intp b1s1 = PyArray_STRIDE((PyArrayObject*)b1,1);

    // copy result from armadillo c++ to numpy array
    for (npy_intp r=0; r<b.n_rows; ++r)
    {
        for (npy_intp c=0; c<b.n_cols; ++c)
        {
            npy_int64 * tmp = (npy_int64*)(b1_data + r*b1s0 + c*b1s1);
            *tmp = b(r,c);
        }
    }

    Py_INCREF(b1);
    return b1;

}

static PyMethodDef commtools_methods[] = {
    {"de2bi", (PyCFunction)ballet_de2bi, METH_VARARGS | METH_KEYWORDS,
    "Convert decimal numbers to binary vectors"},
    {"bi2de", (PyCFunction)ballet_bi2de, METH_VARARGS | METH_KEYWORDS,
    "Convert binary vectors to decimal numbers"},
    {NULL}
};

PyDoc_STRVAR(commtools_documentation,"Python wrapper for ballet::CRCBase");

PyMODINIT_FUNC
init_tools(void)
{

    PyObject *m;
    
    m = Py_InitModule3("_tools",
            commtools_methods,
            commtools_documentation);

    if (m == NULL)
        return;

    import_array();

}
