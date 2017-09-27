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
	PyObject *tmp;

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
				length = repeat*j;
				tmp = Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length);
				PyList_Append(result, tmp);
				Py_DECREF(tmp);
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
	PyObject *tmp;

	if (!PyArg_ParseTuple(args, "siii", &seq, &min, &max, &mrep)){
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
				tmp = Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length);
				PyList_Append(result, tmp);
				Py_DECREF(tmp);
				i = start + length;
				j = min;
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
//static void print_matrix(int **matrix, int size){
//	int i;
//	int j;
//	for(i=0;i<size;i++){
//		for(j=0;j<size;j++){
//			printf("%2d ", *(matrix[i]+j));
//		}
//		printf("\n");
//	}
//}

static int min(int a, int b, int c){
	int d;
	d = a<b?a:b;
	return d<c?d:c;
}

static int** initial_matrix(int size){
	int i;
	int **matrix = (int **)malloc(sizeof(int *)*size);
	for(i=0; i<=size; i++){
		matrix[i] = (int *)malloc(sizeof(int)*size);
	}

	for(i=0; i<=size; i++){
		matrix[i][0] = i;
		matrix[0][i] = i;
	}

	return matrix;
}

static void release_matrix(int **matrix, int size){
	int i;
	for(i=0; i<=size; i++){
		free(matrix[i]);
	}
	free(matrix);
}

static int* build_left_matrix(char *seq, char *motif, int **matrix, int start, int size, int max_error){
	char ref1;
	char ref2;
	int i = 0;
	int j = 0;
	int x = 0;
	int n = 0;
	int prev_min = 0; //previous minimum edit distance (ed)
	int next_min = 0; //current minimum ed
	int top_min = 0; //column minimum ed
	int left_min = 0; //row minimum ed
	int mlen = strlen(motif); //motif length
	int error = 0; //consective errors
	static int res[2]; //result arrary

	for(n=1; n<=size; n++){
		ref1 = seq[start-n];
		ref2 = motif[(mlen-n%mlen)%mlen];
		//min edit
		top_min = matrix[0][n];
		left_min = matrix[n][0];
		for(x=1; x<n; x++){
			if(ref1 == motif[(mlen-x%mlen)%mlen]){
				matrix[x][n] = matrix[x-1][n-1];
			}else{
				matrix[x][n] = min(matrix[x-1][n-1], matrix[x-1][n], matrix[x][n-1]) + 1;
			}

			if(matrix[x][n] < top_min){
				top_min = matrix[x][n];
			}

			if(ref2 == seq[start-x]){
				matrix[n][x] = matrix[n-1][x-1];
			}else{
				matrix[n][x] = min(matrix[n-1][x-1], matrix[n-1][x], matrix[n][x-1]) + 1;
			}

			if(matrix[n][x] < left_min){
				left_min = matrix[n][x];
			}
		}

		if(ref1 == ref2){
			matrix[n][n] = matrix[n-1][n-1];
		}else{
			matrix[n][n] = min(matrix[n-1][n-1], matrix[n-1][n], matrix[n][n-1]) + 1;
		}

		next_min = min(matrix[n][n], top_min, left_min);

		if(next_min > prev_min){
			error++;
		}else{
			error = 0;
		}

		prev_min = next_min;

		if(error > max_error){
			break;
		}
	}

	if(n>size){
		n--;
	}

	n -= error;

	top_min = matrix[0][n];
	left_min = matrix[n][0];
	
	for(x=1; x<n; x++){
		if(matrix[x][n] < top_min){
			top_min = matrix[x][n];
		}

		if(matrix[n][x] < left_min){
			left_min = matrix[n][x];
		}
	}

	next_min = min(top_min, left_min, matrix[n][n]);

	if(next_min == top_min){
		j = n;
		for(x=1; x<n; x++){
			if(matrix[x][n] == next_min){
				i = x;
				break;
			}	
		}
	}else if(next_min == matrix[n][n]){
		i = n;
		j = n;
	}else{
		i = n;
		for(x=n-1; x>0; x--){
			if(matrix[n][x] == next_min){
				j = x;
				break;
			}
		}
	}

	res[0] = i;
	res[1] = j;
	return res;
}

static int* build_right_matrix(char *seq, char *motif, int **matrix, int start, int size, int max_error){
	char ref1;
	char ref2;
	int i = 0;
	int j = 0;
	int x = 0;
	int n = 0;
	int prev_min = 0; //previous minimum edit distance (ed)
	int next_min = 0; //current minimum ed
	int top_min = 0; //column minimum ed
	int left_min = 0; //row minimum ed
	int mlen = strlen(motif); //motif length
	int error = 0; //consective errors
	static int res[2]; //result arrary

	for(n=1; n<=size; n++){
		ref1 = seq[start+n];
		ref2 = motif[(n-1)%mlen];
		//min edit
		top_min = matrix[0][n];
		left_min = matrix[n][0];
		for(x=1; x<n; x++){
			if(ref1 == motif[(x-1)%mlen]){
				matrix[x][n] = matrix[x-1][n-1];
			}else{
				matrix[x][n] = min(matrix[x-1][n-1], matrix[x-1][n], matrix[x][n-1]) + 1;
			}

			if(matrix[x][n] < top_min){
				top_min = matrix[x][n];
			}

			if(ref2 == seq[start+x]){
				matrix[n][x] = matrix[n-1][x-1];
			}else{
				matrix[n][x] = min(matrix[n-1][x-1], matrix[n-1][x], matrix[n][x-1]) + 1;
			}

			if(matrix[n][x] < left_min){
				left_min = matrix[n][x];
			}
		}

		if(ref1 == ref2){
			matrix[n][n] = matrix[n-1][n-1];
		}else{
			matrix[n][n] = min(matrix[n-1][n-1], matrix[n-1][n], matrix[n][n-1]) + 1;
		}

		next_min = min(matrix[n][n], top_min, left_min);

		if(next_min > prev_min){
			error++;
		}else{
			error = 0;
		}

		prev_min = next_min;

		if(error > max_error){
			break;
		}
	}

	if(n>size){
		n--;
	}

	n -= error;

	top_min = matrix[0][n];
	left_min = matrix[n][0];
	
	for(x=1; x<n; x++){
		if(matrix[x][n] < top_min){
			top_min = matrix[x][n];
		}

		if(matrix[n][x] < left_min){
			left_min = matrix[n][x];
		}
	}
	
	next_min = min(top_min, left_min, matrix[n][n]);

	if(next_min == top_min){
		j = n;
		for(x=1; x<n; x++){
			if(matrix[x][n] == next_min){
				i = x;
				break;
			}	
		}
	}else if(next_min == matrix[n][n]){
		i = n;
		j = n;
	}else{
		i = n;
		for(x=n-1; x>0; x--){
			if(matrix[n][x] == next_min){
				j = x;
				break;
			}
		}
	}

	res[0] = i;
	res[1] = j;
	return res;
}

static int backtrace_matrix(int **matrix, int *diagonal, int *mat, int *sub, int *ins, int *del){
	int i = *diagonal;
	int j = *(diagonal+1);
	int cost;
	int r = j;

	while(i>0 && j>0){
		cost = min(matrix[i][j], matrix[i-1][j], matrix[i][j-1]);
		if(cost == matrix[i][j]){
			if(cost == matrix[i-1][j-1]){
				*mat += 1;
			}else{
				*sub += 1;
			}
			i--;
			j--;
		}else if(cost == matrix[i][j-1]){
			*ins += 1;
			j--;
		}else{
			*del += 1;
			i--;
		}
	}

	return r;
}

//search imperfect ssr method
static PyObject *search_issr(PyObject *self, PyObject *args)
{
	int i;
	int j;
	char *seq;
	size_t seqlen;
	int seed_start;
	int seed_end;
	int seed_length;
	int seed_repeat;
	int seed_repeats;
	int seed_minlen;
	int max_errors;
	int size;
	char motif[7] = "\0";
	int start;
	int end;
	int extend_start;
	int extend_len;
	int extend_max_len;
	int *extend_end;
	int length;
	int matches;
	int substitution;
	int insertion;
	int deletion;
	int required_score;
	int mis_penalty;
	int gap_penalty;
	int score;

	PyObject *result = PyList_New(0);
	PyObject *tmp;

	if (!PyArg_ParseTuple(args, "siiiiiii", &seq, &seed_repeats, &seed_minlen, &max_errors, &mis_penalty, &gap_penalty, &required_score, &size)){
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

				//seed end is the same to seed start 0-based coodinates
				seed_end = seed_start + seed_length-seed_length%j - 1;
				matches = seed_length-seed_length%j;
				insertion = 0;
				deletion = 0;
				substitution = 0;

				//extend left
				extend_start = seed_start;
				extend_max_len = extend_start;
				if(extend_max_len > size){
					extend_max_len = size;
				}

				extend_end = build_left_matrix(seq, motif, matrix, extend_start, extend_max_len, max_errors);
				extend_len = backtrace_matrix(matrix, extend_end, &matches, &substitution, &insertion, &deletion);
				start = extend_start - extend_len + 1;

				//extend right
				extend_start = seed_end;
				extend_max_len = seqlen - extend_start - 1;
				if(extend_max_len > size){
					extend_max_len = size;
				}
				//printf("%s\n", "right");
				extend_end = build_right_matrix(seq, motif, matrix, extend_start, extend_max_len, max_errors);
				//printf("%s\n", "back2");
				extend_len = backtrace_matrix(matrix, extend_end, &matches, &substitution, &insertion, &deletion);
				end = extend_start + extend_len + 1;
				
				//printf("%s\n", "extend right yes");

				length = end - start + 1;
				
				score = matches - substitution*mis_penalty - (insertion+deletion)*gap_penalty;
				
				if(score>=required_score){
					tmp = Py_BuildValue("(siiiiiiiii)", motif, j, start, end, length, matches, substitution, insertion, deletion, score);
					PyList_Append(result, tmp);
					Py_DECREF(tmp);
					i = end;
					j = 0;
				}else{
					i = seed_start;
				}
			}else{
				i = seed_start;
			}
		}
	}

	release_matrix(matrix, size);
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
