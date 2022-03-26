/*  Wrapping KDTW similarity function with the Python-C-API. */

// YOU ARE POSSIBLY NEEDED TO CORRECT THE PATHS TO NUMPY HEADERS TO INCLUDE Python.h and arrayobject.h AS I DID BELOW
#include </usr/include/python3.9/Python.h>
#include </usr/local/lib/python3.9/dist-packages/numpy/core/include/numpy/arrayobject.h>
#include <math.h>

/* ==== powered Minkowski Distance   ======================
*/
static double powMinkowski(const double *a, const double *b, int dim, int degree){
double x, out=0.0;
for (int k=0; k<dim; k++){
  x=fabs(a[k]-b[k]);
  out+=pow(x, degree);
  }
return out;
}

static double powMinkowskiZ(const double *a, int dim, int degree){
double x, out=0.0;
for (int k=0; k<dim; k++){
  x=fabs(a[k]);
  out+=pow(x, degree);
  }
return out;
}
/* ==== Evaluate TWED  ======================
*/
double CTWED(double **ta, int la, double *tsa, double **tb, int lb, double *tsb, int dim, double nu, double lambda, int degree) {
	// 
	//    TWED PAMI
	//    
	if(la<0||lb<0){
		fprintf(stderr, "twed: the lengths of the input timeseries should be greater or equal to 0\n");
		exit(-1);
		}
	int r = la;
	int c = lb;
	int deg = degree;
	double disti1, distj1;
	int i,j;
	double dist;

	// allocations
	double **D = (double **)calloc(r+1, sizeof(double*));
	double *Di1 = (double *)calloc(r+1, sizeof(double));
	double *Dj1 = (double *)calloc(c+1, sizeof(double));
	
	double dmin, htrans, dist0;
	
	for(i=0; i<=r; i++) {
		D[i]=(double *)calloc(c+1, sizeof(double));
	}
	// local costs initializations
	for(j=1; j<=c; j++) {
  		distj1=0;
		if(j>1){
			distj1+=powMinkowski(tb[j-2],tb[j-1],dim,deg);	
		}
		else distj1+=powMinkowskiZ(tb[j-1],dim,deg);
   		
		Dj1[j]=distj1;
	}
	
	for(i=1; i<=r; i++) { 
   		disti1=0;
		if(i>1)
			disti1+=powMinkowski(ta[i-2],ta[i-1],dim,deg);
		else disti1+=powMinkowskiZ(ta[i-1],dim,deg);

  		Di1[i]=disti1;
  
  		for(j=1; j<=c; j++) {
  			(dist)=0;
			(dist)+=powMinkowski(ta[i-1],tb[j-1],dim,deg);	
			if(i>1&&j>1)
				(dist)+=powMinkowski(ta[i-2],tb[j-2],dim,deg);

    			D[i][j]=dist;
  		}
	}// for i

	// border of the cost matrix initialization
	D[0][0]=0;
	for(i=1; i<=r; i++)
  		D[i][0]=INFINITY;
	for(j=1; j<=c; j++)
  		D[0][j]=INFINITY;


	for (i=1; i<=r; i++){ 
  		for (j=1; j<=c; j++){
			htrans=fabs((double)(tsa[i-1]-tsb[j-1]));
			if(j>1&&i>1)
				htrans+=fabs((double)(tsa[i-2]-tsb[j-2]));
			dist0=D[i-1][j-1]+D[i][j]+(nu)*htrans;
			dmin=dist0;
   			if(i>1)
   				htrans=((double)(tsa[i-1]-tsa[i-2]));
			else htrans=(double)tsa[i-1];
      			(dist)=Di1[i]+D[i-1][j]+(lambda)+(nu)*htrans;
   			if(dmin>(dist)){
   				dmin=(dist);
			}
   			if(j>1)
   				htrans=((double)(tsb[j-1]-tsb[j-2]));
			else htrans=(double)tsb[j-1]; 
      			(dist)=Dj1[j]+D[i][j-1]+(lambda)+(nu)*htrans; 
   			if(dmin>(dist)){
   				dmin=(dist);
			} 
  			D[i][j] = dmin;
  		}

	}

	dist = D[r][c];
	
	// freeing
	for(i=0; i<=r; i++) {
		free(D[i]);
   	}

	free(D);
	free(Di1);
	free(Dj1);
	return dist;
}


/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
double **ptrvector(long n)  {
	double **v;
	v=(double **)malloc((size_t) (n*sizeof(double)));
	if (!v)   {
		printf("In **ptrvector. Allocation of memory for double array failed.");
		exit(0);  }
	return v;
}

/* ==== Create 1D Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)  {
	return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

/* ==== Create Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory is allocated!                                    */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin, int *L, int *D)  {
	double **c, *a;
	int i,n,m;
	
	n=arrayin->dimensions[0];
	m=arrayin->dimensions[1];
	*L=n;
	*D=m;
	c=ptrvector(n);
	a=(double *) arrayin->data;  /* pointer to arrayin data as double */
	for ( i=0; i<n; i++)  {
		c[i]=a+i*m;  }
	/*for ( i=0; i<n; i++)  { 
		for (int j=0; j<m; j++) 
			printf("%f ", c[i][j]);
		printf("\n");
	}*/
	return c;
}

/* ==== Error tracing ======================
*/
PyObject*  failure(int errid, char* mess){
  printf("%s\n",mess);
  return Py_BuildValue("d", -1.0);
}

/* ==== Python/C bydings for TWED distance computation ======================
*/
static PyObject* distance(PyObject* self, PyObject* args)
{
    int degree = 2;
    double answer=0.0;
    double **A, *tA;
    double **B, *tB;
    double nu, lambda;
    PyArrayObject *seq1, *tseq1, *seq2, *tseq2;
    
    if (!PyArg_ParseTuple(args, "OOOOdd",  &seq1,  & tseq1, &seq2, &tseq2, &nu, &lambda))
        return failure(-1, "TPyArg_ParseTuple error.");

    if (PyArray_DESCR(seq1)->type_num != NPY_DOUBLE)
        return failure(-1, "Type np.float64 expected for p array.");

    if (PyArray_DESCR(seq1)->type_num != NPY_DOUBLE)
        return failure(-1, "Type np.float64 expected for M array.");

    if (PyArray_NDIM(seq1)!=2)
        return failure(-1, "p must be a 2 dimensionnal array.");
    if (PyArray_NDIM(seq2)!=2)
        return failure(-1, "p must be a 2 dimensionnal array.");

    int la, lb, dima, dimb;
    A = pymatrix_to_Carrayptrs(seq1, &la, &dima);
    tA = pyvector_to_Carrayptrs(tseq1);
    B = pymatrix_to_Carrayptrs(seq2, &lb, &dimb);
    tB = pyvector_to_Carrayptrs(tseq2);

    if (dima != dimb){
        printf("dimensions of time series are not equal! \n");
        return Py_BuildValue("f", -1.0);
        }

   answer = CTWED(A, la, tA, B, lb, tB, dima, nu, lambda, degree);

   free(A);
   free(B);

   /*  construct the output, from c double to python float */
   return Py_BuildValue("d", answer);
}

/*  define functions in module */
static PyMethodDef distanceMethods[] =
{
     {"distance", distance, METH_VARARGS, "evaluate kdtw with double 2D-arrays args"},
     {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "TWED", "Some documentation",
    -1,
    distanceMethods
};

PyMODINIT_FUNC
PyInit_TWED(void)
{
    return PyModule_Create(&cModPyDem);
}

#else

/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initTWED(void)
{
    (void) Py_InitModule("TWED", distanceMethods);
}

#endif
