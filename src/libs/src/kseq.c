#include <zlib.h>
#include <Python.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

kseq_t *SEQ;


static char upper(char c){
	if(c>='a' && c<='z'){
		return 'A'+c-'a';
	}else{
		return c;
	}
}

static int upper_string(char *str){
	int i;
	for(i=0; str[i]; i++){
		str[i] = upper(str[i]);
	}
	return i;
}

static PyObject *open_fasta(PyObject *self, PyObject *args){
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}
	gzFile fp;
	fp = gzopen(fasta_path, "rb");
	SEQ = kseq_init(fp);
	return Py_BuildValue("i", 1);
}

static PyObject *close_fasta(PyObject *self, PyObject *args){
	if(SEQ){
		kseq_destroy(SEQ);
	}
	return Py_BuildValue("i", 1);
}

static PyObject *clean_seq(PyObject *self, PyObject *args){
	char *seq;
	if (!PyArg_ParseTuple(args, "s", &seq)){
		return NULL;
	}
	int i;
	int j = 0;
	for(i=0; seq[i]; i++){
		if(!isspace(seq[i])){
			seq[j++] = upper(seq[i]);
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

static PyObject *sub_seq(PyObject *self, PyObject *args){
	char *seq;
	int start;
	int end;

	if (!PyArg_ParseTuple(args, "sii", &seq, &start, &end)){
		return NULL;
	}
	int i;
	int j = 0;
	int real_pos = 0;
	int flag;
	for(i=0; seq[i]; i++){
		flag = isspace(seq[i]);

		if(!flag){
			real_pos++;
		}

		if(real_pos > end){
			break;
		}

		if(real_pos >= start){
			if(!flag){
				seq[j++] = upper(seq[i]);
			}
		}
	}
	seq[j] = '\0';
	return Py_BuildValue("s", seq);
}

static PyObject *iter_seq(PyObject *self, PyObject *args){
	int l;
	while((l=kseq_read(SEQ))>=0){
		upper_string(SEQ->seq.s);
		return Py_BuildValue("(ss)", SEQ->name.s, SEQ->seq.s);
	}
	return Py_BuildValue("");
}

static PyObject *build_index(PyObject *self, PyObject *args){
	char *fasta_path;
	if (!PyArg_ParseTuple(args, "s", &fasta_path)){
		return NULL;
	}
	int position = 0;
	int start = 0;
	int seqlen = 0;
	int gc = 0;
	int ns = 0;
	int c;
	kseq_t *seq;
	kstream_t *ks;
	gzFile fp;
	PyObject *result = PyList_New(0);
	PyObject *tmp;
	
	fp = gzopen(fasta_path, "rb");
	seq = kseq_init(fp);
	ks = seq->f;
	while((c=ks_getc(ks))!=-1){
		position++;
		if(c == 62){
			if(start){
				tmp = Py_BuildValue("(siiiii)", seq->name.s, start, position-start-1, seqlen, gc, ns);
				PyList_Append(result, tmp);
				Py_DECREF(tmp);
			}
			position += ks_getuntil(ks, 0, &seq->name, &c);
			position++;
			while(c != 10){
				c = ks_getc(ks);
				position++;
			}
			start = position;
			seqlen = 0;
			gc = 0;
			ns = 0;
		}else{
			if(c != 10 && c != 13){
				seqlen++;
				c = toupper(c);
				if(c == 71 || c == 67){
					gc++;
				}else if(c == 78){
					ns++;
				}
			}
		}
	}
	tmp = Py_BuildValue("(siiiii)", seq->name.s, start, position-start, seqlen, gc, ns);
	PyList_Append(result, tmp);
	Py_DECREF(tmp);
	return result;
}

static PyMethodDef kseq_methods[] = {
	{"build_index", build_index, METH_VARARGS},
	{"open_fasta", open_fasta, METH_VARARGS},
	{"close_fasta", close_fasta, METH_VARARGS},
	{"iter_seq", iter_seq, METH_VARARGS},
	{"clean_seq", clean_seq, METH_VARARGS},
	{"sub_seq", sub_seq, METH_VARARGS},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef kseq_definition = {
	PyModuleDef_HEAD_INIT,
	"kseq",
	"Read sequence from fasta file",
	-1,
	kseq_methods
};

PyMODINIT_FUNC PyInit_kseq(void){
	Py_Initialize();
    return PyModule_Create(&kseq_definition);
}
