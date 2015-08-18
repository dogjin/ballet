#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include <ballet/comm/crcgenerator.h>

typedef struct
{
    PyObject_HEAD
    ballet::CRCGenerator *d;
} crcgenerator;

static PyObject *
crcgenerator_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{

    crcgenerator *self;
    
    self = (crcgenerator *)type->tp_alloc(type,0);
    if (self != NULL)
    {
        self->d = new ballet::CRCGenerator;
    }

    return (PyObject*)self;

}

static int
crcgenerator_init(crcgenerator *self, PyObject *args, PyObject *kwds)
{ return 0; }

static void
crcgenerator_dealloc(crcgenerator *self)
{
    delete self->d;
    self->ob_type->tp_free((PyObject*)self);
}

static PyMethodDef crcgenerator_methods[] = {
    {NULL,NULL,0,NULL}          /* sentinel */
};

static PyTypeObject crcgenerator_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                                  /*ob_size*/
    "CRCGenerator",                                     /*tp_name*/
    sizeof(crcgenerator),                               /*tp_basicsize*/
    0,                                                  /*tp_itemsize*/
    (destructor)crcgenerator_dealloc,                   /*tp_dealloc*/
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
    "CRCGenerator object",                              /* tp_doc */
    0,                                                  /* tp_traverse */
    0,                                                  /* tp_clear */
    0,                                                  /* tp_richcompare */
    0,                                                  /* tp_weaklistoffset */
    0,                                                  /* tp_iter */
    0,                                                  /* tp_iternext */
    crcgenerator_methods,                               /* tp_methods */
    0,                                                  /* tp_members */
    0,                                                  /* tp_getset */
    0,                                                  /* tp_base */
    0,                                                  /* tp_dict */
    0,                                                  /* tp_descr_get */
    0,                                                  /* tp_descr_set */
    0,                                                  /* tp_dictoffset */
    (initproc)crcgenerator_init,                        /* tp_init */
    0,                                                  /* tp_alloc */
    crcgenerator_new,                                   /* tp_new */
};

static PyMethodDef crcgenerator_module_methods[] = {
    {NULL}
};

PyDoc_STRVAR(crcgenerator_module_documentation,"Python wrapper for ballet::CRCGenerator");

PyMODINIT_FUNC
init_crcgenerator(void)
{

    PyObject *m;
    
    if (PyType_Ready(&crcgenerator_type) < 0)
        return;

    m = Py_InitModule3("_crcgenerator",
            crcgenerator_module_methods,
            crcgenerator_module_documentation);

    if (m == NULL)
        return;

    Py_INCREF(&crcgenerator_type);
    PyModule_AddObject(m,"CRCGenerator",(PyObject*)&crcgenerator_type);

    import_array();

}

