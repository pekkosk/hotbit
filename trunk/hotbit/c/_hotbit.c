/*
 * HOTBIT native extensions
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#include <numpy/arrayobject.h>

#include "geig.h"
#include "slako.h"


static PyMethodDef hotbit_methods[] = {
    { "geig", py_geig, METH_VARARGS,
      "Solve the generalized eigenvalue equation." },
    { "free_geig_workspace", py_free_geig_workspace, METH_NOARGS,
      "Free the workspace of the generalized eigensolver." },
    { "fast_slako_transformations", py_fast_slako_transformations, METH_VARARGS,
      "Slater-Koster transformation." },
    { NULL, NULL, 0, NULL }
};


PyMODINIT_FUNC
init_hotbit(void)
{
  PyObject *m;

  import_array();

  m = Py_InitModule3("_hotbit", hotbit_methods,
                     "HOTBIT native extensions.");
  
  if (!m)
    return;
}
