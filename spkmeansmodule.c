#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <assert.h>
#include "spkmeans.c"

int k, max_iter, dimension, numOfVectors = 0, changes = 1;
float rawK, rawMaxIter;
double **vectors, **centroids;
int **clusters, *clustersSizes;

static PyObject* fit(PyObject *self, PyObject *args){
    int i, j;
    int counter = 1;
    PyObject *pyCentroids;
    PyObject *pyVectors;
    PyObject *tempVec = NULL;
    PyObject *tempCentroid = NULL;
    PyObject *resCentroids;

    if (!PyArg_ParseTuple(args,"OiiOii",&pyCentroids, &k, &max_iter, &pyVectors, &numOfVectors, &dimension)){
        return NULL;
    }
    
    vectors = (double **)calloc(numOfVectors, dimension*sizeof(double));
    assert(vectors != NULL);
    centroids = (double **)calloc(k, dimension*sizeof(double));
    assert(centroids != NULL);
    
    for (i = 0; i < k; i++) {
        centroids[i] = (double *)calloc(dimension, sizeof(double));
        assert(centroids[i] != NULL);
        tempVec = PyList_GetItem(pyCentroids,i);
        for (j = 0; j < dimension; j++) {
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(tempVec,j));  
        }
    } 

    for (i = 0; i < numOfVectors; i++) {
        vectors[i] = (double *)calloc(dimension, sizeof(double));
        assert(vectors[i] != NULL);
        tempVec = PyList_GetItem(pyVectors,i);
        for (j = 0; j < dimension; j++) {
            vectors[i][j] = PyFloat_AsDouble(PyList_GetItem(tempVec,j));  
        }
    } 

    clusters = (int **)calloc(k, numOfVectors*sizeof(int));
    assert(clusters != NULL);
    while ((counter<=max_iter) && (changes > 0)) {
        assignVectorToCluster();
        updateCentroidValue();
        counter += 1;
    }

    resCentroids = PyList_New(0);
    for (i=0; i<k; i++){
        tempCentroid = PyList_New(0);
        for (j=0; j<dimension; j++){
            PyList_Append(tempCentroid,PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_Append(resCentroids, tempCentroid);
    }
    
    free(vectors);
    free(centroids);
    free(clusters);
    free(clustersSizes);
    
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
