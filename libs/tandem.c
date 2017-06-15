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
static int max(int a, int b, int c){
	int d;
	d = a>b?a:b;
	return d>c?d:c;
}

static int min(int a, int b, int c){
	int d;
	d = a<b?a:b;
	return d<c?d:c;
}
//extend left deletion, insertion and mismatch
static int extend_left_mis(char *seq, int start, int left, int mlen, char *motif, int l, int m){
	int i = 0;
	int j = 0;
	left += l;
	for(i = left; i>=0; i--){
		j = ((i - start + m) % mlen + mlen) % mlen;
		if(seq[i] != motif[j]){
			return ++i;
		}
	}
	return i;
}

static int extend_left_sub(char *seq, int start, int left, int mlen, char *motif){
	return extend_left_mis(seq, start, left, mlen, motif, -1, 0);
}

static int extend_left_ins(char *seq, int start, int left, int mlen, char *motif){
	return extend_left_mis(seq, start, left, mlen, motif, -1, 1);
}

static int extend_left_del(char *seq, int start, int left, int mlen, char *motif){
	return extend_left_mis(seq, start, left, mlen, motif, 0, -1);
}

//extend right deletion, insertion and mismatch
static int extend_right_mis(char *seq, size_t seqlen, int start, int right, int mlen, char *motif, int r, int m){
	int i=0;
	int j=0;
	right += r;
	for(i=right; i<seqlen; i++){
		j = (i - start + m) % mlen;
		if(seq[i] != motif[j]){
			return --i;
		}
	}
	return i;
}

static int extend_right_del(char *seq, size_t seqlen, int start, int right, int mlen, char *motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 0, 1);
}

static int extend_right_sub(char *seq, size_t seqlen, int start, int right, int mlen, char *motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, 0);
}

static int extend_right_ins(char *seq, size_t seqlen, int start, int right, int mlen, char *motif){
	return extend_right_mis(seq, seqlen, start, right, mlen, motif, 1, -1);
}


//search imperfect ssr method
static PyObject *search_issr(PyObject *self, PyObject *args)
{
	int i;
	int j;
	char *seq;
	size_t seqlen;
	int start;
	int start1;
	int start2;
	int end;
	int seed_start;
	int seed_length;
	int repeat;
	int score;
	int seed_repeats;
	int seed_minlen;
	int continuous_errors;
	int max_errors;
	int max_match;
	int min_pos;
	int max_pos;
	int required_score;
	char motif[7] = "\0";
	int left;
	int right;
	
	int matches;
	int substitution;
	int deleteion;
	int insertion;

	int lsub_pos;
	int lins_pos;
	int ldel_pos;
	int lsub_match;
	int lins_match;
	int ldel_match;

	int rsub_pos;
	int rins_pos;
	int rdel_pos;
	int rsub_match;
	int rins_match;
	int rdel_match;

	PyObject *result = PyList_New(0);

	if (!PyArg_ParseTuple(args, "siiii", &seq, &seed_repeats, &seed_minlen, &max_errors, &required_score)){
		return NULL;
	}

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
			repeat = seed_length/j;
			if(repeat >= seed_repeats && seed_length >= seed_minlen)
			{
				strncpy(motif, seq+seed_start, j);
				motif[j] = '\0';
				//PyList_Append(result, Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length));
				//printf("%d,%d,%d,%d\n", j, repeat, start, length);
				//return Py_BuildValue("s", motif);
				left = seed_start - 1;
				right = seed_start + seed_length;
				start1 = seed_start;
				start2 = seed_start;
				start = seed_start;
				end = right - 1;
				matches = seed_length;

				
				continuous_errors = 0;
				substitution = 0;
				deleteion = 0;
				insertion = 0;

				// extend left flank sequence
				while(1){
					continuous_errors++;
					if(left < 0){
						break;
					}

					lsub_pos = extend_left_sub(seq, start1, left, j, motif);
					lins_pos = extend_left_ins(seq, start1, left, j, motif);
					ldel_pos = extend_left_del(seq, start1, left, j, motif);

					lsub_match = left - lsub_pos;
					lins_match = left - lins_pos;
					ldel_match = left - ldel_pos;

					min_pos = min(lsub_pos, lins_pos, ldel_pos);
					max_match = max(lsub_match, lins_match, ldel_match);

					if(max_match>2){
						continuous_errors = 0;
						matches += max_match;
						start = min_pos;
						if(max_match == lsub_match){
							substitution++;
						}
						else if(max_match == lins_match)
						{
							insertion++;
						}else{
							deleteion++;
						}
					}

					if(continuous_errors > max_errors)
					{
						break;
					}

					if(min_pos == lsub_pos){
						left = lsub_pos - 1;
					}
					else if(min_pos == lins_pos)
					{
						left = lins_pos - 1;
						start1--;
					}
					else
					{
						left = ldel_pos - 1;
						start1++;
					}
				}

				// extend right flank sequence
				continuous_errors = 0;
				while(1){
					if(right >= seqlen){
						break;
					}
					continuous_errors++;
					rsub_pos = extend_right_sub(seq, seqlen, start2, right, j, motif);
					rins_pos = extend_right_ins(seq, seqlen, start2, right, j, motif);
					rdel_pos = extend_right_del(seq, seqlen, start2, right, j, motif);

					rsub_match = rsub_pos - right;
					rins_match = rins_pos - right;
					rdel_match = rdel_pos - right;

					max_pos = max(rsub_pos, rins_pos, rdel_pos);
					max_match = max(rsub_match, rins_match, rdel_match);

					if(max_match > 2)
					{
						continuous_errors = 0;
						matches += max_match;
						end = max_pos;
						if(max_match == rsub_match){
							substitution++;
						}
						else if(max_match == rins_match)
						{
							insertion++;
						}else{
							deleteion++;
						}
					}

					if(continuous_errors > max_errors)
					{
						break;
					}

					if(max_pos == rsub_pos){
						right = rsub_pos + 1;
					}
					else if(max_pos == rins_pos)
					{
						right = rins_pos + 1;
						start2++;
					}
					else
					{
						right = rdel_pos + 1;
						start2--;
					}
				}
				int gap = (insertion + deleteion - 1)*2 + 5;
				score = matches - substitution - gap;
				
				if(score>=required_score)
				{
					PyList_Append(result, Py_BuildValue("(siiiiiiiii)", motif, j, start+1, end+1, end-start+1, matches, substitution, insertion, deleteion, score));
					i = end + 1;
				}
				else
				{
					i = seed_start + 1;
				}

				j = 0;
			}
			else
			{
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
