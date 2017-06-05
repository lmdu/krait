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

//search imperfect microsatellites
static int extend_right_mis(char *seq, size_t seqlen, int start, int right, int mlen, char *motif, int r, int m){
	int i;
	int j;
	right += r
	for(i=right; i<seqlen; i++){
		j = (i - start + m) % mlen;
		if(seq[i] != motif[j]){
			return --i;
		}
	}
	return i;
}

static int extend_right_del(char *seq, size_t seqlen, int start, int right, int mlen, char * motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 0, 1);
}

static int extend_right_sub(char *seq, size_t seqlen, int start, int right, int mlen, char * motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, 0);
}

static int extend_right_ins(char *seq, size_t seqlen, int start, int right, int mlen, char * motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, -1);
}

static PyObject *search_issr(PyObject *self, PyObject *args)
{
	int i;
	int j;
	char *seq;
	size_t seqlen;
	int start;
	int length;
	int repeat;
	int score;
	int seed_repeats;
	int seed_minlen;
	int continuous_errors;
	int max_errors;
	int required_score;
	char motif[7] = "\0";
	int substitution;
	int deleteion;
	int insertion;

	PyObject *result = PyList_New(0);

	if (!PyArg_ParseTuple(args, "siiiii", &seq, &seed_repeats, &seed_minlen, &max_errors, &required_score)){
		return NULL;
	}

	seqlen = strlen(seq);

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
			if(repeat >= seed_repeats && length >= seed_minlen)
			{
				strncpy(motif, seq+start, j);
				motif[j] = '\0';
				//PyList_Append(result, Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length));
				//printf("%d,%d,%d,%d\n", j, repeat, start, length);
				//return Py_BuildValue("s", motif);
				right = start + length;
				start2 = start;

				// extend right flank sequence
				continuous_errors = 0;
				substitution = 0;
				deleteion = 0;
				insertion = 0;

				while(1){
					continuous_errors++;
					rsub_pos = extend_right_sub(seq, seqlen, start2, right, j, motif);
					rins_pos = extend_right_ins(seq, seqlen, start2, right, j, motif);
					rdel_pos = extend_right_del(seq, seqlen, start2, right, j, motif);

					rsub_match = rsub_pos - right;
					rins_match = rins_pos - right;
					rdel_match = rdel_pos - right;

					if(rsub_match>0 || rins_match>0 || rdel_match>0){
						continuous_errors = 0;
					}

					if(continuous_errors >= max_errors){
						break;
					}

					if(rsub_pos >= rins_pos){
						if(rsub_pos >= rdel_pos){
							right = rsub_pos + 1;
							substitution++;
						}else{
							right = rdel_pos + 1;
							deleteion++;
							start2--;
						}
					}else if(rins_pos > rdel_pos){
						right = rins_pos + 1;
						insertion++;
						start2++;
					}else{
						right = rdel_pos + 1;
						deleteion++;
						start2--;
					}

					if(right >= seqlen){
						break;
					}

				}

				PyList_Append(result, Py_BuildValue("(siiiii)", motif, j, start+1, start+length, length));





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
