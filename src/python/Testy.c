#include <Python.h>
#include <stdio.h>
#include "libnumarray.h"
#include "numarray.h"

/* Written Nic Wherry circa July 2003 */

#define SIZE 5
//line 6
void changeElem( Int32 * );

void test()
{
  int arr[SIZE] = {1,2,3,4,5}; 
  int *arrPtr;
  int i;
  arrPtr=&arr;
  changeElem(arrPtr);
  for (i=0; i<SIZE; i++)
    printf("%d ", arrPtr[i]);
}
//line19
void changeElem( Int32 *arrPtr )
{
  int i;
  for (i=0; i<SIZE; i++)
  {
    *arrPtr = 2 * *arrPtr;
    arrPtr++;
  }
}
//line29
static PyObject *
Testy_changeElem(PyObject *self, PyObject *args)
{
  PyObject *firstPtr;
  PyArrayObject *secPtr;
  Int32 *data;
  
  if(!PyArg_ParseTuple(args, "O", &firstPtr))
    return NULL;

  secPtr = NA_InputArray(firstPtr, tInt32, C_ARRAY);  
  data = (Int32 *)NA_OFFSETDATA(secPtr);
  changeElem(data);
  PyArray_XDECREF(secPtr);
  return (PyObject *)Py_BuildValue("O", secPtr);
}
//line46
static PyObject *
Testy_test(PyObject *self, PyObject *args)
{
  test();
  return (PyObject*)Py_BuildValue("");
}

static PyMethodDef
TestyMethods[] =
{
  { "changeElem", Testy_changeElem, METH_VARARGS },
  { "test", Testy_test, METH_VARARGS },
  { NULL, NULL },
};

void initTesty()
{
  Py_InitModule("Testy", TestyMethods);
}
