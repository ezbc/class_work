#include <Python.h>
#include <pspline.h>

/* http://dan.iel.fm/posts/python-c-extensions/ Great module for c wrappers in
 * python */

/* Docstrings */
static char module_docstring[] =
    "This module provides an interface for calling c function pspline"
static char pspline_docstring[] =
    "an O(n) spline smoother with penalty on D^m This version can save intermediate results for reruns with new values of smoothing parameter LAMBDA and can compute the GCV, CV, and DF criteria"

/* */
static PyObject *pspline_pspline(PyObject *self, PyObject *args);

/* Method definition This second line contains all the info that the
 * interpreter needs to link a Python call to the correct C function and call
 * it in the right way. The first string is the name of the function as it will
 * be called from Python, the second object is the C function to link to and
 * the last argument is the docstring for the function. The third argument
 * METH_VARARGS means that the function only accepts positional arguments. If
 * you wanted to support keyword arguments, you would need to change this to
 * METH_VARARGS | METH_KEYWORDS. */
static PyMethodDef module_methods[] = {
    {"pspline", pspline_pspline, METH_KEYWORDS, pspline_docstring},
    {NULL, NULL, 0, NULL}
};


/* Initialize C module */
PyMODINIT_FUNC init_pspline(void)
{
    PyObject *m = Py_InitModule3("_pspline", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
}

/* Write the pspline_pspline function */
static PyObject *pspline_pspline(PyObject *self, PyObject *args)
{
    pspline







