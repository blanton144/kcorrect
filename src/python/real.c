#include <Python.h>
#include <stdio.h>
#include "libnumarray.h"
#include "numarray.h"
#define SIZE 5
//line6
void printArray( Int32 * );

/* Written Nic Wherry circa July 2003 */


void test(void)
{
  Int32 arr[SIZE] = {1, 2, 3, 4, 5};
  printArray(arr);
}
//line14
void printArray( Int32 *arrPtr )
{
  int i;
  for (i=0; i<SIZE; i++)
  {
    printf( "%d ", arrPtr[i] );
  }
}
//line23
static PyObject *
real_printArray(PyObject *self, PyObject *args)
{
  PyObject *firstPtr;
  PyArrayObject *secPtr;
  Int32 *data;

  if(!PyArg_ParseTuple(args, "O", &firstPtr))
  {  
    printf("bad input");
    fflush(stdout);
  }

  secPtr = NA_InputArray( firstPtr, tInt32, C_ARRAY); //HERE!

  data=(Int32 *) NA_OFFSETDATA(secPtr);
  printArray(data);
  return (PyObject *)Py_BuildValue("");
}
//line39
static PyObject *
real_test(PyObject *self, PyObject *args)
{
  test();
  return (PyObject *)Py_BuildValue("");
}
//line46
static PyMethodDef
realMethods[] = 
{
  { "printArray", real_printArray, METH_VARARGS },
  { "test", real_test, METH_VARARGS },
  { NULL, NULL },
};
//line54
void initreal(void)
{
  Py_InitModule("real", realMethods);
  import_libnumarray();
}
