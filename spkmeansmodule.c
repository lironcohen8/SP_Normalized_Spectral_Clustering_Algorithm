#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include "spkmeans.c"

/*
int k, max_iter, dimension, numOfVectors = 0, changes = 1;
float rawK, rawMaxIter;
char *goal;
double **vectors, **centroids;
int **clusters, *clustersSizes;
*/

static void freeMemory() {
    free(vectors);
    free(centroids);
    free(clusters);
    free(clustersSizes);
}

static PyObject* fit(PyObject *self, PyObject *args){
    int i, j;
    int counter = 1;
    PyObject *pyCentroids;
    PyObject *pyVectors;
    PyObject *tempVec = NULL;
    PyObject *tempCentroid = NULL;
    PyObject *resCentroids;

    if (!PyArg_ParseTuple(args,"OiiOsii",&pyCentroids, &k, &max_iter, &pyVectors, &goal, &numOfVectors, &dimension)){
        return NULL;
    }
    
    vectors = (double **)calloc(numOfVectors, dimension*sizeof(double));
    errorAssert(vectors != NULL,0);
    centroids = (double **)calloc(k, dimension*sizeof(double));
    errorAssert(centroids != NULL,0);
    
    for (i = 0; i < k; i++) {
        centroids[i] = (double *)calloc(dimension, sizeof(double));
        errorAssert(centroids[i] != NULL,0);
        tempVec = PyList_GetItem(pyCentroids,i);
        for (j = 0; j < dimension; j++) {
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(tempVec,j));  
        }
    } 

    for (i = 0; i < numOfVectors; i++) {
        vectors[i] = (double *)calloc(dimension, sizeof(double));
        errorAssert(vectors[i] != NULL,0);
        tempVec = PyList_GetItem(pyVectors,i);
        for (j = 0; j < dimension; j++) {
            vectors[i][j] = PyFloat_AsDouble(PyList_GetItem(tempVec,j));  
        }
    } 

    if (strcmp(goal,"spk")==0){
        int calcK = eigengapHeuristic();
        if (k==0) {
            k = calcK;
        }
        createUMatrix();
        assignUToVectors();
        initCentroids();
        clusters = (int **)calloc(k, numOfVectors*sizeof(int));
        errorAssert(clusters != NULL,0);
        while ((counter <= max_iter) && (changes > 0)) {
            assignVectorToCluster();
            updateCentroidValue();
            counter += 1;
        }
    }
    else if (strcmp(goal,"wam")==0){
        printMatrix(weightedAdjacencyMatrix(),numOfVectors,numOfVectors);
        freeMemory();
        Py_RETURN_NONE;
    } 
    else if (strcmp(goal,"ddg")==0){
        printMatrix(diagonalDegreeMatrix(1,1),numOfVectors,numOfVectors);
        freeMemory();
        Py_RETURN_NONE;
    } 
    else if (strcmp(goal,"lnorm")==0){
        printMatrix(laplacianNorm(),numOfVectors,numOfVectors);
        freeMemory();
        Py_RETURN_NONE;
    } 
    else if (strcmp(goal,"jacobi")==0){
        jacobi(vectors, 1);
        freeMemory();
        Py_RETURN_NONE;
    } 
    else{
        errorAssert(0==1,1); /*If the goal is unknown*/
        freeMemory();
        Py_RETURN_NONE;
    }

    resCentroids = PyList_New(0);
    for (i=0; i<k; i++){
        tempCentroid = PyList_New(0);
        for (j=0; j<dimension; j++){
            PyList_Append(tempCentroid,PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_Append(resCentroids, tempCentroid);
    }
      
    freeMemory();
    return resCentroids;
}

static PyMethodDef kmeansMethods[] = {
    {"fit",
    (PyCFunction) fit,
    METH_VARARGS,
    PyDoc_STR("Kmeans")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeasns",     // name of module exposed to Python
    "spkmeasns Python wrapper for custom C extension library.", // module documentation
    -1,
    kmeansMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
