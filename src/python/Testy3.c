#include <Python.h>
#include <stdio.h>
#include "libnumarray.h"
#include "numarray.h"
#define SIZE=5
//line 6
void copyIt(void);

/* Written Nic Wherry circa July 2003 */


void test(void)
{
  printf("This is only a test");
}
//line11
void copyIt(void)
{
  printf("copied");
}
//line16
static PyObject *
Testy3_copyIt(PyObject *self, PyObject *args)
{
  PyObject *firstPtr, *copy;
  PyArrayObject *secPtr;
  
  if(!PyArg_ParseTuple(args, "O", &firstPtr))
    printf("bad");
  
  secPtr=NA_InputArray(firstPtr, tInt32, C_ARRAY);
  copy = PyArray_Copy(secPtr);
  return copy;
}
//line30
static PyObject*
Testy3_test(PyObject *self, PyObject *args)
{
  test();
  return (PyObject*)Py_BuildValue("");
}
//line37
static PyMethodDef
Testy3Methods[] = 
{
  { "copyIt", Testy3_copyIt, METH_VARARGS },
  { "test", Testy3_test, METH_VARARGS },
  { NULL, NULL },
};
//line45
void initTesty3(void)
{
  Py_InitModule("Testy3", Testy3Methods);
  import_libnumarray();
}
