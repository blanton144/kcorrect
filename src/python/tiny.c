#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Written Nic Wherry circa July 2003 */

int add(int m, int n)
{
  return(m+n);
}

int subtract(int m, int n)
{
  return (m-n);
}

void test()
{
  printf("%d + %d = %d", 5, 6, add(5,6));
  printf("%d - %d = %d", 7, 3, subtract(7,3));
}

static PyObject * 
tiny_add(PyObject *self, PyObject *args)
{
  int res;
  int num1;
  int num2;
  PyObject* retval;

  res = PyArg_ParseTuple(args, "ii", &num1, &num2);
  if (!res)
    return NULL;
  res = add(num1, num2);
  retval = (PyObject*)Py_BuildValue("i", res);
  return retval;
}

static PyObject * 
tiny_subtract(PyObject *self, PyObject *args)
{
  int res;
  int num1;
  int num2;
  PyObject* retval;

  res = PyArg_ParseTuple(args, "ii", &num1, &num2);
  if (!res)
    return NULL;
  res = subtract(num1, num2);
  retval = (PyObject*)Py_BuildValue("i", res);
  return retval;
}

static PyObject *
tiny_test(PyObject *self, PyObject *args)
{
  test();
  return (PyObject*)Py_BuildValue("");
}

static PyMethodDef
tinyMethods[] = 
{
  { "add", tiny_add, METH_VARARGS },
  { "subtract", tiny_subtract, METH_VARARGS },
  { "test", tiny_test, METH_VARARGS },
  { NULL, NULL },
};

void inittiny()
{
  Py_InitModule("tiny", tinyMethods);
}
