#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <ballet/comm/crcdetector.h>
#include <ballet/comm/de2bi.h>

typedef struct
{
    PyObject_HEAD
    ballet::CRCDetector *d;
} crcdetector;

static PyObject *
crcdetector_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{

    crcdetector *self;
    
    self = (crcdetector *)type->tp_alloc(type,0);
    if (self != NULL)
    {
        self->d = new ballet::CRCDetector;
    }

    return (PyObject*)self;

}

static int
crcdetector_init(crcdetector *self, PyObject *args, PyObject *kwds)
{ return 0; }

static void
crcdetector_dealloc(crcdetector *self)
{
    delete self->d;
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
crcdetector_setproperty(crcdetector *self, PyObject *args, PyObject *kwds)
{

    PyObject * propertyName;
    PyObject * propertyValue;
    char * property_name;

    static char * kwlist[] = {"name","value",NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds,"OO",kwlist,&propertyName,
        &propertyValue))
    {
        return NULL;
    }

    if (!PyString_Check(propertyName))
    {
        PyErr_Format(PyExc_TypeError,
            "PropertyName must be a string, not %.500s",
            Py_TYPE(propertyName)->tp_name);
        return 0;
    }

    // convert :propertyName: to c string
    property_name = PyString_AsString(propertyName);

    if (!strcmp(property_name,"Polynomial"))
    {

        if (!PyArray_Check(propertyValue))
        {
            PyErr_Format(PyExc_TypeError,
                "CRCDetector.Polynomial property must"
                "by a numpy.ndarray type, not %.500s",
                Py_TYPE(propertyValue)->tp_name);
            return 0;
        }

        // cast to PyArrayObject type
        PyArrayObject * p = (PyArrayObject*)propertyValue;

        // get pointer to the data buffer
        char * bytes = PyArray_BYTES(p);

        // get the number of array dimensions
        npy_intp size = PyArray_Size(propertyValue);
        npy_intp * shape = PyArray_SHAPE(p);
        npy_intp dim1 = shape[0];
        npy_intp dim2 = shape[1];

        // get array stride length
        npy_intp s0 = PyArray_STRIDE(p,0);
        npy_intp s1 = PyArray_STRIDE(p,1);
     
        // resize CRCDetector.Polynomial property
        self->d->Polynomial.set_size(dim1,dim2);

        for (npy_intp r=0; r<dim1; ++r)
        {

            for (npy_intp c=0; c<dim2; ++c)
            {
                char * tmp = bytes + r*s0 + c*s1;
                self->d->Polynomial(r,c) = *((npy_int64*)tmp);
            }

        }

    }
    else if (!strcmp(property_name,"ReflectChecksums"))
    {

        self->d->ReflectChecksums = PyObject_IsTrue(propertyValue);

    }
    else
    {
        PyErr_Format(PyExc_AttributeError,
            "CRCDetector has no property %s",
            property_name);
        return 0;
    }

    Py_RETURN_NONE;

}

static PyObject *
crcdetector_getproperty(crcdetector *self, PyObject *args, PyObject *kwds)
{

    PyObject * propertyName;
    char * property_name;

    static char * kwlist[] = {"name",NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds,"O",kwlist,&propertyName))
    {
        return NULL;
    }

    if (!PyString_Check(propertyName))
    {
        PyErr_Format(PyExc_TypeError,
            "PropertyName must be a string, not %.500s",
            Py_TYPE(propertyName)->tp_name);
        return 0;
    }

    // convert :propertyName: to c string
    property_name = PyString_AsString(propertyName);

    if (!strcmp(property_name,"Polynomial"))
    {

        // create return array
        int ndim = 2;
        npy_intp dims[ndim];
        dims[0] = self->d->Polynomial.n_rows;
        dims[1] = self->d->Polynomial.n_cols;
        PyObject * poly = (PyObject*)PyArray_NewFromDescr(&PyArray_Type,
            PyArray_DescrFromType(NPY_INT64),
            ndim,
            dims,
            NULL,
            NULL,
            NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE,
            NULL);

        // get pointer to array data buffer
        char * bytes = PyArray_BYTES((PyArrayObject*)poly);
        npy_intp s0 = PyArray_STRIDE((PyArrayObject*)poly,0);
        npy_intp s1 = PyArray_STRIDE((PyArrayObject*)poly,1);

        // copy result from armadillo c++ to numpy array
        for (npy_intp r=0; r<dims[0]; ++r)
        {
            for (npy_intp c=0; c<dims[1]; ++c)
            {
                char * tmp = bytes + r*s0 + c*s1;
                (*((npy_int64*)tmp)) = self->d->Polynomial(r,c);
            }
        }

        Py_INCREF(poly);
        return poly;

    }
    
    else if (!strcmp(property_name,"ReflectChecksums"))
    {

        return PyBool_FromLong(self->d->ReflectChecksums);

    }

    else
    {

        PyErr_Format(PyExc_AttributeError,
            "CRCDetector has no property %s",
            property_name);
        return 0;

    }

}

static PyObject *
crcdetector_islocked(crcdetector *self)
{

    // get the status of the system object
    bool status = self->d->isLocked();
    return PyBool_FromLong(status);

}

static PyObject *
crcdetector_release(crcdetector *self)
{

    // release system object
    self->d->release();
    Py_RETURN_NONE;

}

static PyObject *
crcdetector_reset(crcdetector *self)
{

    // release system object
    self->d->reset();
    Py_RETURN_NONE;

}

static PyObject *
crcdetector_step(crcdetector *self, PyObject *args)
{

    PyObject * arg001;

    if (!PyArg_ParseTuple(args,"O",&arg001))
    {
        return 0;
    }

    if (!PyArray_Check(arg001))
    {
        PyErr_Format(PyExc_TypeError,
            "Input must be a numpy.ndarray type, not %.500s",
            Py_TYPE(arg001)->tp_name);
        return 0;
    }

    // cast to PyArrayObject type
    PyArrayObject * X = (PyArrayObject*)arg001;

    // get pointer to the data buffer
    char * bytes = PyArray_BYTES(X);

    // get the number of array dimensions
    npy_intp * shape = PyArray_SHAPE(X);
    npy_intp dim1 = shape[0];
    npy_intp dim2 = shape[1];

    // get array strides
    npy_intp s0 = PyArray_STRIDE(X,0);
    npy_intp s1 = PyArray_STRIDE(X,1);
 
    // create armadillo matrix
    arma::imat X001(dim1,dim2);

    // copy data from numpy input array to
    // armadillo matrix for processing
    for (npy_intp r=0; r<dim1; ++r)
    {

        for (npy_intp c=0; c<dim2; ++c)
        {
            char * tmp = bytes + r*s0 + c*s1;
            X001(r,c) = *((npy_int64*)tmp);
        }

    }

    // run CRCDetector::detect()
    int err = self->d->detect(X001);

    // release system object
    return PyInt_FromLong(err);

}

static PyMethodDef crcdetector_methods[] = {
    {"setProperty", (PyCFunction)crcdetector_setproperty, METH_VARARGS | METH_KEYWORDS,
    "Set the property value"},
    {"getProperty", (PyCFunction)crcdetector_getproperty, METH_VARARGS | METH_KEYWORDS,
    "Get the property value"},
    {"isLocked", (PyCFunction)crcdetector_islocked, METH_NOARGS,
    "Locked status for input attributes and nontunable properties"},
    {"release", (PyCFunction)crcdetector_release, METH_NOARGS,
    "Allow property value and input characteristics changes"},
    {"reset", (PyCFunction)crcdetector_reset, METH_NOARGS,
    "Reset states of CRC detector object"},
    {"step", (PyCFunction)crcdetector_step, METH_VARARGS,
    "Detect errors in input data using CRC"},
    {NULL,NULL,0,NULL}          /* sentinel */
};

static PyTypeObject crcdetector_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                                  /*ob_size*/
    "CRCDetector",                                      /*tp_name*/
    sizeof(crcdetector),                                /*tp_basicsize*/
    0,                                                  /*tp_itemsize*/
    (destructor)crcdetector_dealloc,                    /*tp_dealloc*/
    0,                                                  /*tp_print*/
    0,                                                  /*tp_getattr*/
    0,                                                  /*tp_setattr*/
    0,                                                  /*tp_compare*/
    0,                                                  /*tp_repr*/
    0,                                                  /*tp_as_number*/
    0,                                                  /*tp_as_sequence*/
    0,                                                  /*tp_as_mapping*/
    0,                                                  /*tp_hash */
    0,                                                  /*tp_call*/
    0,                                                  /*tp_str*/
    0,                                                  /*tp_getattro*/
    0,                                                  /*tp_setattro*/
    0,                                                  /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,           /*tp_flags*/
    "CRCDetector object",                               /* tp_doc */
    0,                                                  /* tp_traverse */
    0,                                                  /* tp_clear */
    0,                                                  /* tp_richcompare */
    0,                                                  /* tp_weaklistoffset */
    0,                                                  /* tp_iter */
    0,                                                  /* tp_iternext */
    crcdetector_methods,                                /* tp_methods */
    0,                                                  /* tp_members */
    0,                                                  /* tp_getset */
    0,                                                  /* tp_base */
    0,                                                  /* tp_dict */
    0,                                                  /* tp_descr_get */
    0,                                                  /* tp_descr_set */
    0,                                                  /* tp_dictoffset */
    (initproc)crcdetector_init,                         /* tp_init */
    0,                                                  /* tp_alloc */
    crcdetector_new,                                    /* tp_new */
};

static PyObject *
ballet_convertHexToBinary(PyObject *self, PyObject *args, PyObject *kwds)
{

    PyObject * arg001;

    static char * kwlist[] = {"hexval",NULL};

    if (!PyArg_ParseTupleAndKeywords(args,kwds,"O",kwlist,&arg001))
    {
        return NULL;
    }

    if (!PyString_Check(arg001))
    {
        PyErr_Format(PyExc_TypeError,
            "Input must be a str type, not %.500s",
            Py_TYPE(arg001)->tp_name);
        return 0;
    }

    // get c string
    char * hexString = PyString_AsString(arg001);

    // remove '0x'
    hexString = hexString + 2;

    // convert hex string to decimal integer
    unsigned int hexval = static_cast<unsigned int>(strtol(hexString,NULL,16));

    // prepend zeros
    double p = std::floor(std::log(static_cast<double>(hexval)) /
        std::log(16.0)) + 1.0;
    double nextpow = pow(16.0,p);
    hexval = hexval + static_cast<unsigned int>(nextpow);

    // convert to binary vector
    arma::imat binaryVector = ballet::de2bi(hexval);

    // flip binary vector (left-msb)
    binaryVector = arma::fliplr(binaryVector);

    // create return array
    int bvndim = 2;
    npy_intp bvdims[bvndim];
    bvdims[0] = binaryVector.n_rows;
    bvdims[1] = binaryVector.n_cols;
    PyObject * bv = (PyObject*)PyArray_NewFromDescr(&PyArray_Type,
        PyArray_DescrFromType(NPY_INT64),
        bvndim,
        bvdims,
        NULL,
        NULL,
        NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE,
        NULL);

    // get pointer to array data buffer
    char * bvdata = PyArray_BYTES((PyArrayObject*)bv);
    npy_intp bvs0 = PyArray_STRIDE((PyArrayObject*)bv,0);
    npy_intp bvs1 = PyArray_STRIDE((PyArrayObject*)bv,1);

    // copy result from armadillo c++ to numpy array
    for (npy_intp r=0; r<binaryVector.n_rows; ++r)
    {
        for (npy_intp c=0; c<binaryVector.n_cols; ++c)
        {
            char * tmp = bvdata + r*bvs0 + c*bvs1;
            (*((npy_int64*)tmp)) = binaryVector(r,c);
        }
    }

    Py_INCREF(bv);
    return bv;

}


static PyMethodDef crcdetector_module_methods[] = {
    {"convertHexToBinary", (PyCFunction)ballet_convertHexToBinary, METH_VARARGS | METH_KEYWORDS,
    "Convert hex string to binary number"},
    {NULL}
};

PyDoc_STRVAR(crcdetector_module_documentation,"Python wrapper for ballet::CRCDetector");

PyMODINIT_FUNC
init_crcdetector(void)
{

    PyObject *m;
    
    if (PyType_Ready(&crcdetector_type) < 0)
        return;

    m = Py_InitModule3("_crcdetector",
            crcdetector_module_methods,
            crcdetector_module_documentation);

    if (m == NULL)
        return;

    Py_INCREF(&crcdetector_type);
    PyModule_AddObject(m,"CRCDetector",(PyObject*)&crcdetector_type);

    import_array();

}

