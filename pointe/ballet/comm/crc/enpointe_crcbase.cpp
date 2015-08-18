#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <ballet/comm/crcbase.h>

typedef struct
{
    PyObject_HEAD
    ballet::CRCBase *d;
} crcbase;

static PyObject *
crcbase_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{

    crcbase *self;
    
    self = (crcbase *)type->tp_alloc(type,0);
    if (self != NULL)
    {
        self->d = new ballet::CRCBase;
    }

    return (PyObject*)self;

}

static int
crcbase_init(crcbase *self, PyObject *args, PyObject *kwds)
{ return 0; }

static void
crcbase_dealloc(crcbase *self)
{
    delete self->d;
    self->ob_type->tp_free((PyObject*)self);
}

static PyMethodDef crcbase_methods[] = {
    {NULL,NULL,0,NULL}          /* sentinel */
};

static PyTypeObject crcbase_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                                  /*ob_size*/
    "CRCBase",                                          /*tp_name*/
    sizeof(crcbase),                                    /*tp_basicsize*/
    0,                                                  /*tp_itemsize*/
    (destructor)crcbase_dealloc,                        /*tp_dealloc*/
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
    "CRCBase object",                                    /* tp_doc */
    0,                                                  /* tp_traverse */
    0,                                                  /* tp_clear */
    0,                                                  /* tp_richcompare */
    0,                                                  /* tp_weaklistoffset */
    0,                                                  /* tp_iter */
    0,                                                  /* tp_iternext */
    crcbase_methods,                                    /* tp_methods */
    0,                                                  /* tp_members */
    0,                                                  /* tp_getset */
    0,                                                  /* tp_base */
    0,                                                  /* tp_dict */
    0,                                                  /* tp_descr_get */
    0,                                                  /* tp_descr_set */
    0,                                                  /* tp_dictoffset */
    (initproc)crcbase_init,                             /* tp_init */
    0,                                                  /* tp_alloc */
    crcbase_new,                                        /* tp_new */
};

static PyMethodDef crcbase_module_methods[] = {
    {NULL}
};

PyDoc_STRVAR(crcbase_module_documentation,"Python wrapper for ballet::CRCBase");

PyMODINIT_FUNC
init_crcbase(void)
{

    PyObject *m;
    
    if (PyType_Ready(&crcbase_type) < 0)
        return;

    m = Py_InitModule3("_crcbase",
            crcbase_module_methods,
            crcbase_module_documentation);

    if (m == NULL)
        return;

    Py_INCREF(&crcbase_type);
    PyModule_AddObject(m,"CRCBase",(PyObject*)&crcbase_type);

    import_array();

}

