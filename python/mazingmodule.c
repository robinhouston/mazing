#include <Python.h>
#include <gmp.h>

static PyObject *
maze_by_index(PyObject *self, PyObject *args)
{
    PyObject *index_obj;
    if (!PyArg_ParseTuple(args, "O!", PyLong_Type, &index_obj))
        return NULL;
    return index_obj;
}

static PyMethodDef MazingMethods[] = {
    {"maze_by_index", maze_by_index, METH_VARARGS, "Find a maze by index."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initmazing(void)
{
    (void) Py_InitModule("mazing", MazingMethods);
}
