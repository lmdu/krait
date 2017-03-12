#include <Python.h>

//search perfect microsatellites
static PyObject *search_ssr(PyObject *self, PyObject *args)
{
	char *seq;
	int mono = 0;
	int di = 0;
	int tri = 0;
	int tetra = 0;
	int penta = 0;
	int hexa = 0;
	int rep[6];

	size_t len;
	int start;
	int length;
	int repeat;
	int i;
	int j;
	char motif[7] = "\0";

	PyObject *result = PyList_New(0);

	if (!PyArg_ParseTuple(args, "s(iiiiii)", &seq, &mono, &di, &tri, &tetra, &penta, &hexa)){
		return NULL;
	}

	rep[0] = mono;
	rep[1] = di;
	rep[2] = tri;
	rep[3] = tetra;
	rep[4] = penta;
	rep[5] = hexa;

	len = strlen(seq);

	for (i=0; i<len; i++)
	{
		if (seq[i] == 78)
		{
			continue;
		}

		for (j=1; j<=6; j++)
		{
			start = i;
			length = j;
			while(start+length<len && seq[i]==seq[i+j] && seq[i]!=78){
				i++;
				length++;
			}
			repeat = length/j;
			if(repeat>=rep[j-1])
			{
				strncpy(motif, seq+start, j);
				motif[j] = '\0';
				PyList_Append(result, Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length));
				//printf("%d,%d,%d,%d\n", j, repeat, start, length);
				//return Py_BuildValue("s", motif);
				i = start + length;
				j = 0;
			}
			else
			{
				i = start;
			}
		}
	}
	return result;
};

//search perfect satellite variable number tandem repeat
static PyObject *search_vntr(PyObject *self, PyObject *args)
{
	char *seq;
	int max;
	int min;
	int mrep;

	size_t len;
	int start;
	int length;
	int repeat;
	int i;
	int j;
	char *motif;

	PyObject *result = PyList_New(0);

	if (!PyArg_ParseTuple(args, "siii", &seq, &max, &min, &mrep)){
		return NULL;
	}

	len = strlen(seq);

	for (i=0; i<len; i++)
	{
		if (seq[i] == 78)
		{
			continue;
		}

		for (j=min; j<=max; j++)
		{
			start = i;
			length = j;
			while(start+length<len && seq[i]==seq[i+j] && seq[i]!=78){
				i++;
				length++;
			}
			repeat = length/j;
			if(repeat>=mrep)
			{
				motif = (char *)malloc(j+1);
				strncpy(motif, seq+start, j);
				motif[j] = '\0';
				PyList_Append(result, Py_BuildValue("(siiii)", motif, j, repeat, start+1, start+length));
				i = start + length;
				j = 0;
			}
			else
			{
				i = start;
			}
		}
	}
	return result;
};


static PyMethodDef add_methods[] = {
	{"search_ssr", search_ssr, METH_VARARGS},
	{"search_vntr", search_vntr, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};


void inittandem()
{
	PyObject *m;
	m = Py_InitModule("tandem", add_methods);
	if(m == NULL){
		return;
	}
}
