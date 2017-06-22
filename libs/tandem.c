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
				length = repeat * j;
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
				length = j*repeat;
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
static int min(int a, int b, int c){
	int d;
	d = a<b?a:b;
	return d<c?d:c;
}

static int** initial_matrix(int size){
	int i;
	int j;
	int **matrix = (int **)malloc(sizeof(int *)*size);
	for(i=0; i<size; i++)
		matrix[i] = (int *)malloc(sizeof(int)*size);

	for(i=0; i<size; i++)
		matrix[i][0] = i;

	for(j=0; j<size; j++)
		matrix[0][j] = j;

	return matrix;
}

static int build_matrix(char *seq, char *motif, int **matrix, int start, int size, int max_error){
	int error = 0;
	int i = 0;
	int j = 0;
	int m = strlen(motif);
	for(i=1; i<size; i++){
		for(j=1; j<i; j++){
			//fill columns in k row
			if(seq[start+j-1] == motif[(i-1)%m]){
				matrix[i][j] = matrix[i-1][j-1];
			}else{
				matrix[i][j] = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) + 1;
			}

			//fill rows in k column
			if(seq[start+i-1] == motif[(j-1)%m]){
				matrix[j][i] = matrix[j-1][i-1];
			}else{
				matrix[j][i] = min(matrix[j-1][i-1], matrix[j-1][i], matrix[j][i-1]) + 1;
			}

		}
		// fill the diagonal edit distance
		if(seq[start+i-1] == motif[(j-1)%m]){
			matrix[i][i] = matrix[i-1][i-1];
			error = 0;
		}else{
			matrix[i][i] = min(matrix[i-1][i-1], matrix[i-1][i], matrix[i][i-1]) + 1;
			error++;
		}

		if(error>max_error){
			break;
		}
	}
	return i-error;
}

static void backtrace_matrix(int **matrix, int diagonal, int *mat, int *sub, int *ins, int *del){
	//static int match[] = {0,0,0,0};
	int i = diagonal;
	int j = diagonal;
	while(i>0 && j>0){
		if(matrix[i][j] == matrix[i][j-1]+1){
			j--;
			*ins += 1;
		}else if(matrix[i][j] == matrix[i-1][j]+1){
			i--;
			*del += 1;
		}else{
			if(matrix[i][j] == matrix[i-1][j-1]){
				*mat += 1;
			}else{
				*sub += 1;
			}
			i--;
			j--;
		}
	}
}


//search imperfect ssr method
static PyObject *search_issr(PyObject *self, PyObject *args)
{
	int i;
	int j;
	char *seq;
	size_t seqlen;
	int seed_start;
	int seed_length;
	int seed_repeat;
	int seed_repeats;
	int seed_minlen;
	int max_errors;
	int required_identity;
	int size;
	char motif[7] = "\0";
	//int start;
	int end;
	int extend_start;
	int extend_len;
	int extend_end;
	int length;
	int matches;
	int substitution;
	int insertion;
	int deletion;
	float identity;

	PyObject *result = PyList_New(0);

	if (!PyArg_ParseTuple(args, "siiiii", &seq, &seed_repeats, &seed_minlen, &max_errors, &required_identity, &size)){
		return NULL;
	}

	//create edit distance matrix
	int **matrix = initial_matrix(size);
	
	seqlen = strlen(seq);

	for (i=0; i<seqlen; i++)
	{
		if (seq[i] == 78)
		{
			continue;
		}

		for (j=1; j<=6; j++)
		{
			seed_start = i;
			seed_length = j;
			while(seed_start+seed_length<seqlen && seq[i]==seq[i+j] && seq[i]!=78)
			{
				i++;
				seed_length++;
			}
			seed_repeat = seed_length/j;
			if(seed_repeat >= seed_repeats && seed_length >= seed_minlen)
			{
				strncpy(motif, seq+seed_start, j);
				motif[j] = '\0';
				//PyList_Append(result, Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length));
				//printf("%d,%d,%d,%d\n", j, repeat, start, length);
				//return Py_BuildValue("s", motif);
				matches = seed_length;
				insertion = 0;
				deletion = 0;
				substitution = 0;

				//extend
				extend_start = seed_start+seed_length;
				extend_len = seqlen - extend_start;
				if(extend_len > size){
					extend_len = size;
				}
				extend_end = build_matrix(seq, motif, matrix, extend_start, extend_len, max_errors);
				backtrace_matrix(matrix, extend_end, &matches, &substitution, &insertion, &deletion);
				end = extend_start + extend_end;
				length = end - seed_start + 1;
				identity = matches/length;
				
				if(identity>=required_identity){
					PyList_Append(result, Py_BuildValue("(siiiiiiiii)", motif, j, seed_start+1, end+1, length, matches, substitution, insertion, deletion, identity));
					i = end;
				}else{
					i = seed_start;
				}

				j = 0;
			}else{
				i = seed_start;
			}
		}
	}
	return result;
};


static PyMethodDef add_methods[] = {
	{"search_ssr", search_ssr, METH_VARARGS},
	{"search_vntr", search_vntr, METH_VARARGS},
	{"search_issr", search_issr, METH_VARARGS},
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
