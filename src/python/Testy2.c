#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Written Nic Wherry circa July 2003 */


double dubb( double elem )
{
  return (2.0*elem);
}

void test()
{
  printf("Goes from %f to %f", 5.0, dubb(5.0));
}
//line15
static PyObject *
Testy2_dubb(PyObject *self, PyObject *args)
{
  double res;
  double elem;
  PyObject *newElem;
  
  res = PyArg_ParseTuple(args, "d", &elem);
  if(!res)
    return NULL;
  res = dubb(elem);
  newElem = (PyObject*)Py_BuildValue("d", res);
  return newElem;
}
//line30
static PyObject*
Testy2_test(PyObject*self, PyObject *args)
{
  test();
  return (PyObject*)Py_BuildValue("");
}

static PyMethodDef
Testy2Methods[] = 
{
  { "dubb", Testy2_dubb, METH_VARARGS },
  { "test", Testy2_test, METH_VARARGS },
  { NULL, NULL },
};

void initTesty2()
{
  Py_InitModule("Testy2", Testy2Methods );
}
      
