/*
 * Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
 * Please see the accompanying LICENSE file for further information.
 */

/*
 * HOTBIT native extensions module
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#include <numpy/arrayobject.h>

#include "geig.h"
#include "slako.h"
#include "spherical.h"
#include "multipole.h"


static PyMethodDef hotbit_methods[] = {
    { "geig", py_geig, METH_VARARGS,
      "Solve the generalized eigenvalue equation." },
    { "free_geig_workspace", py_free_geig_workspace, METH_NOARGS,
      "Free the workspace of the generalized eigensolver." },
    { "fast_slako_transformations", py_fast_slako_transformations, METH_VARARGS,
      "Slater-Koster transformation." },
    { "cartesian_to_spherical", py_cartesian_to_spherical, METH_VARARGS,
      "Transform cartesian into spherical coordinates." },
    { "solid_harmonic_R", py_solid_harmonic_R, METH_VARARGS,
      "Compute the regular solid harmonics." },
    { "solid_harmonic_I", py_solid_harmonic_I, METH_VARARGS,
      "Compute the irregular solid harmonics." },
    { "multipole_to_multipole", py_multipole_to_multipole, METH_VARARGS,
      "Multipole-to-multipole transformation, i.e. shift the center of the multipole "
      "expansion to a new position." },
    { "multipole_to_local", py_multipole_to_local, METH_VARARGS,
      "Multipole-to-local transformation, i.e. compute the expansion of the field "
      "due to a multipole at a certain distance." },
    { "local_to_local", py_local_to_local, METH_VARARGS,
      "Local-to-local transformation, i.e. shift the expansion origin to a new positions." },
    { "transform_multipole", py_transform_multipole, METH_VARARGS,
      "Transform a multipole moment given a linear cartesian transformation matrix." },
    { NULL, NULL, 0, NULL }
};

/*
 * Module initialization
 */

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

/*
 * Module declaration
 */

#if PY_MAJOR_VERSION >= 3
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
        ob = PyModule_Create(&moduledef);
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        ob = Py_InitModule3(name, methods, doc);
#endif

MOD_INIT(_hotbit)
{
    PyObject *m;

    import_array();

    MOD_DEF(m, "_hotbit", hotbit_methods,
            "HOTBIT native C extensions.");

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
