#define PY_SSIZE_T_CLEAN
/*#include <Python.h>*/
#include <stdio.h>
#include <assert.h>
#include <math.h>

int k, max_iter, dimension, numOfVectors = 0, changes = 1;
float rawK, rawMaxIter;
double **vectors, **centroids, **wam, **ddg, **lnorm;
int **clusters, *clustersSizes;

void *calloc(size_t nitems, size_t size);
void *malloc(size_t size);
void *realloc(void *ptr, size_t size);
void free(void *ptr);
char *strtok(char * str, const char *delim);
double atof(const char * str);
int strcmp (const char* str1, const char* str2);


int calcDimension(char buffer[]) {
    /*
    gets a buffer contains the first vector and calculates 
    */
    int i, dimension; 
    char c;
    dimension = 0;
    i = 0;
    c = buffer[i];
    while (c != '\n') {
        if (c == ',') {
            dimension++;
        }
        c = buffer[++i];
    }
    buffer[i] = '\0';
    return dimension+1;
}

void readFile(FILE *file) {
    /*Reading the input file and put the data into the 'vectors' list*/
    int j, sizeFull = 1;
    char *vectorStr, buffer[1000];
    double **tmp;
    vectorStr = (char *)calloc(dimension, 100*sizeof(char));
    assert(vectorStr != NULL);
    vectors = (double **)malloc(1 * sizeof(*vectors));
    assert(vectors != NULL);
    fgets(buffer,1000,file);
    dimension = calcDimension(buffer);
    do {
        if (numOfVectors == sizeFull) {
            sizeFull *= 2;
            tmp = realloc(vectors, sizeFull * sizeof(*vectors));
            vectors = tmp;
        }
        vectorStr = strtok(buffer, ",");
        j = 0;
        vectors[numOfVectors] = (double *)calloc(dimension, sizeof(double)); 
        assert(vectors[numOfVectors] != NULL);
        while (vectorStr != NULL) {
            vectors[numOfVectors][j] = atof(vectorStr);
            vectorStr = strtok(NULL, ",");
            j++;
        }
        numOfVectors++;
    }
    while (fgets(buffer,1000,file) != NULL);
    free(vectorStr);
}

void initCentroids() {
    /*Initialize the clusters and their centroids from the first K vectors*/
    int i,j;
    centroids = (double **)calloc(k, dimension*sizeof(double));
    assert(centroids != NULL);
    assert(k < numOfVectors);
    for (i = 0; i < k; i++) {
        centroids[i] = (double *)calloc(dimension, sizeof(double)); 
        assert(centroids[i] != NULL);
        for (j = 0; j < dimension; j++) {
            centroids[i][j] = vectors[i][j];
        }
    }
}

double distance(double *vector1, double *vector2) {
    /*'''Claculates the distance between two vectors*/
    double dis = 0;
    int i;
    for (i = 0; i < dimension; i++) {
        dis += (vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
    }
    return dis;
}

int closestCentroid(double *vector) {
    /*Finds the closest centroid to a vector by the distance function*/
    double minDis, dis;
    int minCenInd,i;
    
    minDis = distance(vector, centroids[0]); /*Initiate the minimum distance to be the distance from the first centroid*/
    minCenInd = 0; /*Initiate the closest centroid to be the first one*/
    
    for (i = 0; i < k; i++) { /*For each centroid (there are K)*/
        dis = distance(vector, centroids[i]);
        if (dis < minDis) {
            minDis = dis;
            minCenInd = i;
        }
    }   
    return minCenInd;
}
 
void assignVectorToCluster() {
    /*Finds the closest centroid for each vector and adds
    the vector to the centroids cluster*/
    int i, newCentroidInd, clusterSize;
    int * cluster;
    clustersSizes = (int *)calloc(k, sizeof(int));
    assert(clustersSizes != NULL);

    /*Clearing all clusters (we do not want to remember what was here)*/
    for (i = 0; i < k; i++) { 
        clusters[i] = (int *)calloc(numOfVectors, sizeof(int));
        assert(clusters[i] != NULL);
    }
        
    for (i = 0; i < numOfVectors; i++) {
        newCentroidInd = closestCentroid(vectors[i]); /*Finds the closest centroid*/
        cluster = clusters[newCentroidInd];
        clusterSize = clustersSizes[newCentroidInd];
        cluster[clusterSize] = i; /*Adds the vector to the appropriate cluster*/
        clustersSizes[newCentroidInd]++;
    }
}

double* calcCentroidForCluster(int clusterInd) {
    /*Calculates the centroid for a given cluster*/
    int numOfVectorsInCluster, i, j;
    int * cluster;
    double * sumVector = (double *)calloc(dimension, sizeof(double));
    assert(sumVector != NULL);
    numOfVectorsInCluster = clustersSizes[clusterInd];
    cluster = clusters[clusterInd];
    
    for (i = 0; i < dimension; i++) {
        for (j = 0; j < numOfVectorsInCluster; j++) {
            sumVector[i] += (vectors[cluster[j]])[i];
        }
    }

    for (i = 0; i < dimension; i++) {
        sumVector[i] /= numOfVectorsInCluster; /*Replace the sum with the average*/
    }
    return sumVector;
}

void updateCentroidValue() {
    /*Updates the centroid value for each cluster and checks if 
    it is different then what we had before*/
    int i, j;
    double * newValue;
    changes = 0;
    for (i = 0; i < k; i++) {
        newValue = calcCentroidForCluster(i);
        for (j = 0; j < dimension; j++) {
            if (newValue[j] != centroids[i][j]) { /*If the centroid changed*/
                changes += 1;
            }    
            centroids[i][j] = newValue[j];
        }
    }
}

void printResult() {
    /*Prints the centroids*/
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < dimension; j++) {
            printf("%.4f", centroids[i][j]); /*Format the floats precision to 4 digits*/
            if (j < dimension - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}



double** matrixMultiplication(double** a, double** b){
    int i,j,k;
    double** mul = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(mul != NULL);
    for (i = 0; i < numOfVectors; i++) {
        mul[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(mul[i] != NULL);
    }

    for(i = 0; i < numOfVectors; i++){    
        for(j = 0; j < numOfVectors; j++){    
            mul[i][j]=0;    
            for(k = 0; k < numOfVectors; k++){    
                mul[i][j]+=a[i][k]*b[k][j];    
            }    
        }    
    } 
    return mul;   
}

double calcWeightsForAdjacencyMatrix(double *vector1, double *vector2){
    double dis = 0;
    int i;
    for (i = 0; i < dimension; i++) { /*'''Claculates the euclidean distance between two vectors*/
        dis += pow(vector1[i]-vector2[i],2);
    }
    dis = sqrt(dis);
    dis = -0.5*dis;
    
    return exp(dis);
} 

double** weightedAdjacencyMatrix(){
    int i, j;

    wam = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(wam != NULL);
    for (i = 0; i < numOfVectors; i++) {
        wam[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(wam[i] != NULL);
    }

    for (i = 0; i < numOfVectors; i++){
        double* vector1 = vectors[i];
        for (j = i+1; j < numOfVectors; j++){
            double* vector2 = vectors[j];
            wam[i][j] = calcWeightsForAdjacencyMatrix(vector1, vector2);
            wam[j][i] = wam[i][j];
        }
    }

    return wam;
}

double** diagonalDegreeMatrix(int calcWam, int toPrint){
    int i,j;

    if (calcWam==1){
        wam = weightedAdjacencyMatrix();
    }

    ddg = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(ddg != NULL);
    for (i = 0; i < numOfVectors; i++) {
        ddg[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(ddg[i] != NULL);
    }

    for (i = 0; i < numOfVectors; i++) {
        double sum = 0;
        for (j = 0; j < numOfVectors; j++){
            sum += wam[i][j];
        }

        if (toPrint==1){
            ddg[i][i] = sum;
        }
        else{
            ddg[i][i] = 1/sqrt(sum);
        }  
    }
        
    return ddg;
} 

double** laplacianNorm(){
    int i,j;

    wam = weightedAdjacencyMatrix();
    ddg = diagonalDegreeMatrix(0,0);

    lnorm = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(lnorm != NULL);
    for (i = 0; i < numOfVectors; i++) {
        lnorm[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(lnorm[i] != NULL);
    }

    lnorm =  matrixMultiplication(matrixMultiplication(ddg, wam),ddg);

    for (i = 0; i < numOfVectors; i++){
        for (j = 0; j < numOfVectors; j++){
            if (i==j){
                lnorm[i][j] = 1-lnorm[i][j];
            }
            else{
                lnorm[i][j] = (-1)*lnorm[i][j];
            }

        }
    }
    return lnorm;
}

int* maxOffDiagonalValue(double** mat){
    int i,j;
    int maxRow = 0;
    int maxCol = 1;
    int* res;

    for (i = 0; i < numOfVectors; i++){
        for (j = i+1; j < numOfVectors; j++){
            if (abs(mat[i][j])>abs(mat[maxRow][maxCol])){
                maxRow = i;
                maxCol = j;
            }
        }
    }
    
    res = malloc(2*sizeof(int));
    res[0] = maxRow;
    res[1] = maxCol;
    return res;
}

double calcTheta(int i, int j){
    return (lnorm[j][j]-lnorm[i][i])/(2*lnorm[i][j]);
}

double calcT(double theta){
    int sign = theta < 0 ? -1 : 1;
    double denom = abs(theta)+sqrt(pow(theta,2)+1);
    return sign/denom;
}

double calcC(double t){
    return 1 / (sqrt(pow(t,2)+1));
}

double calcS(double t, double c){
    return t*c;
}

double** createRotationMatrixP(double** A){
    int i;
    int* maxValInd = maxOffDiagonalValue(A);
    int maxRow = maxValInd[0];
    int maxCol = maxValInd[1];
    double theta = calcTheta(maxRow, maxCol);
    double t = calcT(theta);
    double c = calcC(t);
    double s = calcS(t, c);

    double** P = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(P != NULL);
    for (i = 0; i < numOfVectors; i++) {
        P[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(P[i] != NULL);
    }

    for (i = 0; i < numOfVectors; i++) {
        P[i][i] = 1;
    }
    P[maxRow][maxRow] = c;
    P[maxCol][maxCol] = c;
    P[maxRow][maxCol] = s;
    P[maxCol][maxRow] = (-1)*s;

    return P;
}

void updateAPrime(double** A, double** APrime, int i, int j, double c, double s){
    int r;
    for (r = 0; r < numOfVectors; r++){
        if ((r!=i) && (r!=j)){
            APrime[r][i] = c*A[r][i]-s*A[r][j];
            APrime[r][j] = c*A[r][j]+s*A[r][i];
        }
    }
    APrime[i][i] = pow(c,2)*A[i][i]+pow(s,2)*A[j][j]-2*s*c*A[i][j];
    APrime[j][j] = pow(s,2)*A[i][i]+pow(c,2)*A[j][j]+2*s*c*A[i][j];
    APrime[i][j] = 0;
    APrime[j][i] = 0;
}

double calcOffSquared(double** mat){
    int i,j;
    double sum;
    for (i = 0; i < numOfVectors; i++){
        for (j = 0; j < numOfVectors; j++){
            if (i!=j){
                sum += pow(mat[i][j],2);
            }
        }
    }
    return sum;
}

int checkConvergence(double** A, double** APrime){
    double epsilon = 0.001;
    double a = calcOffSquared(A);
    double ap = calcOffSquared(APrime);

    if ((a-ap)<=epsilon){
        return 1;
    }
    return 0;
}

double** jacobi(){
    int count = 0;
    int isConverged = 0;
    int* maxValInd;
    int i, maxRow, maxCol;
    double theta, t, c, s;
    double** A;

    double** APrime = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    assert(APrime != NULL);
    for (i = 0; i < numOfVectors; i++) {
        APrime[i] = (double *)calloc(numOfVectors, sizeof(double));
        assert(APrime[i] != NULL);
    }

    A = laplacianNorm();
    do {
        maxValInd = maxOffDiagonalValue(A);
        maxRow = maxValInd[0];
        maxCol = maxValInd[1];
        theta = calcTheta(maxRow, maxCol);
        t = calcT(theta);
        c = calcC(t);
        s = calcS(t, c);

        updateAPrime(A, APrime, maxRow, maxCol, c, s);
        isConverged = checkConvergence(A, APrime);
        A = APrime;
        count++;
    }
    while ((isConverged==0)&&(count<100));

    return A;
}  

void printMatrix(double** mat) {
    /*Prints a matrix*/
    int i, j;
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < numOfVectors; j++) {
            printf("%.4f", mat[i][j]); /*Format the floats precision to 4 digits*/
            if (j < numOfVectors - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}



void printVectors() {
    /*Prints the centroids*/
    int i, j;
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < dimension; j++) {
            printf("%.4f", vectors[i][j]); /*Format the floats precision to 4 digits*/
            if (j < dimension - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}










int main(int argc, char *argv[]) {
    FILE *file;
    char *goal;

    assert(argc == 4); /*Checks if we have the right amount of args*/ 
    
    assert(sscanf(argv[1], "%f", &rawK) == 1);
    k = (int)rawK;
    assert(rawK - k == 0 && k >= 0); /*checks if k is a non-negative int*/
    
    file = fopen(argv[3],"r");
    readFile(file);

    goal = argv[2];

    if (strcmp(goal,"spk")==0){
        printf("d");
    } 
    else if (strcmp(goal,"wam")==0){
        printMatrix(weightedAdjacencyMatrix());
    } 
    else if (strcmp(goal,"ddg")==0){
        printMatrix(diagonalDegreeMatrix(1,1));
    } 
    else if (strcmp(goal,"lnorm")==0){
        printMatrix(laplacianNorm());
    } 
    else if (strcmp(goal,"jacobi")==0){
        printf("d");
    } 
    else{
        assert(0==1); /*If the goal is unknown*/
    }
        


    /*
    int counter = 1;
    assert(argc == 3 || argc == 2); 

    assert(sscanf(argv[1], "%f", &rawK) == 1);
    k = (int)rawK;
    assert(rawK - k == 0 && k > 0); 
    
    max_iter = 200;
    if (argc == 3) {
        assert(sscanf(argv[2], "%f", &rawMaxIter) == 1);
        max_iter = (int)rawMaxIter;
        assert(rawMaxIter - max_iter == 0 && max_iter > 0); 
    }
    
    readFile();
    initCentroids();
    
    clusters = (int **)calloc(k, numOfVectors*sizeof(int));
    assert(clusters != NULL);
    while ((counter<=max_iter) && (changes > 0)) {
        assignVectorToCluster();
        updateCentroidValue();
        counter += 1;
    }

    printResult();
    free(vectors);
    free(centroids);
    free(clusters);
    free(clustersSizes);
    */
    return 0;
}

/*
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
    "mykmeanssp",     // name of module exposed to Python
    "mykmeanssp Python wrapper for custom C extension library.", // module documentation
    -1,
    kmeansMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
*/