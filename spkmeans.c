#include <stdio.h>
#include <assert.h>
#include <math.h>

typedef struct eigenVector {
    double eigenVal;
    int columnIndex;
} eigenVector; 

int k, dimension, numOfVectors = 0, changes = 1, max_iter = 300;
float rawK, rawMaxIter;
double *eigenVals, *eigenGaps;
double **vectors, **centroids, **wam, **ddg, **lnorm, **V, **U;
int **clusters, *clustersSizes;
eigenVector *eigenVectors;

void *calloc(size_t nitems, size_t size);
void *malloc(size_t size);
void *realloc(void *ptr, size_t size);
void free(void *ptr);
char *strtok(char * str, const char *delim);
double atof(const char * str);
int strcmp (const char* str1, const char* str2);
void qsort(void *base, size_t nmemb, size_t size,
           int (*compar)(const void *, const void *));

void errorsAssert(int cond, int isInputError) {
    if (!cond) {
        if (isInputError==1) {
            printf("Invalid Input!");
        }
        else {
            printf("An Error Has Occured");
        }
        assert(0);
    }
}

int calcDimension(char buffer[]) {
    /*gets a buffer contains the first vector and calculates*/
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
    errorsAssert(vectorStr != NULL,0);
    vectors = (double **)malloc(1 * sizeof(*vectors));
    errorsAssert(vectors != NULL,0);
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
        errorsAssert(vectors[numOfVectors] != NULL,0);
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

void assignUToVectors() {
    /*put U vectors in vectors matrix for further calculations*/
    int i;
    free(vectors);
    vectors = (double **)calloc(numOfVectors, k * sizeof(double));
    errorsAssert(vectors != NULL,0);
    for (i = 0; i < k; i++) {
        vectors[i] = (double *)calloc(k, sizeof(double)); 
        errorsAssert(vectors[i] != NULL,0);
    }
    vectors = U;
}

void initCentroids() {
    /*Initialize the clusters and their centroids from the first K vectors*/
    int i,j;
    centroids = (double **)calloc(k, dimension*sizeof(double));
    errorsAssert(centroids != NULL,0);
    errorsAssert(k < numOfVectors,0);
    for (i = 0; i < k; i++) {
        centroids[i] = (double *)calloc(dimension, sizeof(double)); 
        errorsAssert(centroids[i] != NULL,0);
        for (j = 0; j < dimension; j++) {
            centroids[i][j] = vectors[i][j];
        }
    }
}

double distance(double *vector1, double *vector2) {
    /*Calculates the distance between two vectors*/
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
    errorsAssert(clustersSizes != NULL,0);

    /*Clearing all clusters (we do not want to remember what was here)*/
    for (i = 0; i < k; i++) { 
        clusters[i] = (int *)calloc(numOfVectors, sizeof(int));
        errorsAssert(clusters[i] != NULL,0);
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
    errorsAssert(sumVector != NULL,0);
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


void deepClone(double **a, double** b){
    /*clones b matrix values to a*/
    int i,j;
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < numOfVectors; j++) {
            a[i][j] = b[i][j];
        }
    }
}

void printMatrix(double** mat, int numOfRows, int numOfCols) {
    /*prints a matrix*/
    int i, j;
    for (i = 0; i < numOfRows; i++) {
        for (j = 0; j < numOfCols; j++) {
            printf("%.4f", mat[i][j]); /*format the floats precision to 4 digits*/
            if (j < numOfCols - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

double** matrixMultiplication(double** a, double** b){
    /*gets two matrixes and multiplies them*/
    int i,j,k;
    double** mul = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(mul != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        mul[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(mul[i] != NULL,0);
    }

    for(i = 0; i < numOfVectors; i++){    
        for(j = 0; j < numOfVectors; j++){    
            mul[i][j]=0;    
            for(k = 0; k < numOfVectors; k++){
                /*by matrixes multiplication rules*/    
                mul[i][j] += a[i][k] * b[k][j];     
            }    
        }    
    } 
    return mul;   
}

void squareMatrixTranspose(double **matrix, int numOfRows) {
    /*transpose a squared matrix*/
    int i,j,tmp;
    for (i = 1; i < numOfRows; i++) {
        for (j = 0; j < i; j++) {
            tmp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = tmp;
        }
    }
}

double calcWeightsForAdjacencyMatrix(double *vector1, double *vector2){
    /*gets two vectors and calculates wij for them*/
    double dis = 0;
    int i;
    for (i = 0; i < dimension; i++) { 
        /*calculates the euclidean distance between two vectors*/
        dis += pow(vector1[i]-vector2[i],2);
    }
    dis = sqrt(dis); /*added factors according to instructions*/
    dis = -0.5*dis; /*added factors according to instructions*/
    
    return exp(dis);
} 

double** weightedAdjacencyMatrix(){
    /*calculates weighted adjacency matrix after vectors matrix was set up*/
    int i, j;

    wam = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(wam != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        wam[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(wam[i] != NULL,0);
    }

    for (i = 0; i < numOfVectors; i++){
        double* vector1 = vectors[i]; /*gets vector i*/
        for (j = i+1; j < numOfVectors; j++){ 
            double* vector2 = vectors[j]; /*gets vector j*/
            wam[i][j] = calcWeightsForAdjacencyMatrix(vector1, vector2);
            wam[j][i] = wam[i][j]; /*wam is symetric*/
        }
    }

    return wam;
}

double** diagonalDegreeMatrix(int calcWam, int toPrint){
    /*/*calculates diagonal degree matrix*/
    int i,j;

    if (calcWam==1){ /*if wam wasn't calculated before, used in lnorm*/
        wam = weightedAdjacencyMatrix();
    }

    ddg = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(ddg != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        ddg[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(ddg[i] != NULL,0);
    }

    for (i = 0; i < numOfVectors; i++) {
        double sum = 0;
        for (j = 0; j < numOfVectors; j++){
            sum += wam[i][j]; /*sums the row*/
        }

        if (toPrint==1){ /*if was called for ddg goal, only need to be printed*/
            ddg[i][i] = sum;
        }
        else{
            ddg[i][i] = 1/sqrt(sum); /*if was called for further calculations*/
        }  
    }
        
    return ddg;
} 

double** laplacianNorm(){
    /*calculated the laplacian norm matrix*/
    int i,j;

    wam = weightedAdjacencyMatrix(); /*calling wam*/
    ddg = diagonalDegreeMatrix(0,0); /*calling ddg without the need to calculate wam*/

    lnorm = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(lnorm != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        lnorm[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(lnorm[i] != NULL,0);
    }

    /*multiplies according to formula*/
    lnorm =  matrixMultiplication(matrixMultiplication(ddg, wam),ddg); 

    for (i = 0; i < numOfVectors; i++){
        for (j = 0; j < numOfVectors; j++){
            if (i==j){
                lnorm[i][j] = 1-lnorm[i][j]; /*I - matrix*/
            }
            else{
                lnorm[i][j] = (-1)*lnorm[i][j]; /*I - matrix*/
            }
        }
    }
    return lnorm;
}

int* maxOffDiagonalValue(double** mat){
    /*calculates the indexes of max off-diagonal element in a matrix*/
    int i,j;
    int maxRow = 0;
    int maxCol = 1; /*initialized as the first off-diagonal element*/
    int* res;

    for (i = 0; i < numOfVectors; i++){
        for (j = i+1; j < numOfVectors; j++){
            /*finds the max off-diagonal element*/
            if (fabs(mat[i][j])>fabs(mat[maxRow][maxCol])){ 
                maxRow = i;
                maxCol = j;
            }
        }
    }
    
    res = malloc(2*sizeof(int));
    res[0] = maxRow;
    res[1] = maxCol;
    return res; /*returns as a tuple (i,j)*/
}

double calcTheta(double **matrix, int i, int j){
    /*calcs theta as part os jacobi computations*/
    return (matrix[j][j]-matrix[i][i])/(2*matrix[i][j]);
}

double calcT(double theta){
    /*calcs t as part os jacobi computations*/
    int sign = theta < 0 ? -1 : 1;
    double denom = fabs(theta)+sqrt(pow(theta,2)+1);
    return sign/denom;
}

double calcC(double t){
    /*calcs c as part os jacobi computations*/
    return 1 / (sqrt(pow(t,2)+1));
}

double calcS(double t, double c){
    /*calcs s as part os jacobi computations*/
    return t*c;
}

void createRotationMatrixP(double** P, int maxRow, int maxCol, double c, double s){
    /*calcs P rotation matrix as part of jacobi computations*/
    int i,j;
    for (i = 0; i < numOfVectors; i++) { /*I matrix*/
        for (j = 0; j < numOfVectors; j++) { 
            P[i][j] = i==j ? 1 : 0;
        }
    }

    P[maxRow][maxRow] = c;
    P[maxCol][maxCol] = c;
    P[maxRow][maxCol] = s;
    P[maxCol][maxRow] = (-1)*s;
}

void updateAPrime(double** A, double** APrime, int i, int j, double c, double s){
    /*gets A and A' and updates A' by five equations in instructions*/
    int r;
    for (r = 0; r < numOfVectors; r++){
        if ((r!=i) && (r!=j)){
            APrime[i][r] = APrime[r][i] = c*A[r][i]-s*A[r][j];
            APrime[j][r] = APrime[r][j] = c*A[r][j]+s*A[r][i];
        }
    }
    APrime[i][i] = pow(c,2)*A[i][i]+pow(s,2)*A[j][j]-2*s*c*A[i][j];
    APrime[j][j] = pow(s,2)*A[i][i]+pow(c,2)*A[j][j]+2*s*c*A[i][j];
    APrime[i][j] = 0;
    APrime[j][i] = 0;
}

double calcOffSquared(double** mat){
    /*gets a matrix and calculates the sum of off-diagonal elements squared*/
    int i,j;
    double sum = 0;
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
    /*gets A and A' matrixes and checks if their off values are closer than epsilon*/
    double epsilon = 0.001; /*constant from instructions*/
    double a = calcOffSquared(A);
    double ap = calcOffSquared(APrime);
    
    if ((a-ap)<=epsilon){
        return 1;
    }
    return 0;
}

void printJacobi(double **A, double **V) {
    /*gets A matrix (for eigenvalues) and V matrix (for eigenvectors) 
    and prints them according to instructions*/
    int i,j;
    for (i = 0; i < numOfVectors; i++) {
        printf("%.4f", A[i][i]); /*eigenvalues, Format to 4 digits*/
            if (i < numOfVectors - 1) {
                printf(",");
            }
    }
    printf("\n");
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < numOfVectors; j++) {
            printf("%.4f", V[j][i]); /*Transpose V, Format to 4 digits*/
            if (j < numOfVectors - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

double** jacobi(double **A, int toPrint){
    /*calculates jacobi iterations until convergence*/
    int i, maxRow, maxCol, count=0, isConverged=0;
    int* maxValInd;
    double theta, t, c, s;
    double **APrime, **P;

    APrime = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(APrime != NULL,0);
    P = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(P != NULL,0);
    V = (double **)calloc(numOfVectors, numOfVectors*sizeof(double));
    errorsAssert(V != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        APrime[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(APrime[i] != NULL,0);
        P[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(P[i] != NULL,0);
        V[i] = (double *)calloc(numOfVectors, sizeof(double));
        errorsAssert(V[i] != NULL,0);
    }

    for (i = 0; i < numOfVectors; i++) { 
        V[i][i] = 1; /*init V as I matrix for מeutrality to multiplication*/
    }
    deepClone(APrime, A);

    do {        
        maxValInd = maxOffDiagonalValue(A);        
        maxRow = maxValInd[0];
        maxCol = maxValInd[1];
        if (A[maxRow][maxCol] == 0) { /*matrix is already diagonal*/
            break;
        }
        theta = calcTheta(A, maxRow, maxCol);
        t = calcT(theta);
        c = calcC(t);
        s = calcS(t, c);

        createRotationMatrixP(P, maxRow, maxCol, c, s); /*creating P rotation matrix*/
        V = matrixMultiplication(V, P); /*updating eigenvectors matrix*/

        updateAPrime(A, APrime, maxRow, maxCol, c, s); /*updating A'*/
        isConverged = checkConvergence(A, APrime); /*checks convergence*/

        /* A = APrime, deep clone */
        deepClone(A,APrime);
        count++; /*iterations count*/
        }
    while ((isConverged==0)&&(count<100)); /*until convergence or 100 iterations*/

    if (toPrint==0) { /*if further calculations are necessary*/
        return A;
    }
    else { /*if goal was jacobi, only need to be printed*/
        printJacobi(A, V);
        return NULL;
    }
}  

int compareEigenVectors(const void *a, const void *b) {
    /*comperator of vectors for quicksort*/
    struct eigenVector *eva = (struct eigenVector *) a;
    struct eigenVector *evb = (struct eigenVector *) b;
    int diff = eva->eigenVal - evb->eigenVal; /*by eigen values*/
    return diff != 0 ? diff : eva->columnIndex - evb->columnIndex; /*then, by index*/
}

void sortEigenVectorsAndValues() {
    /*sorts eigen vecctors using quicksort 
    and sorts eigen values accordingly*/
    int i;
    eigenVectors = (eigenVector *)calloc(numOfVectors, numOfVectors*sizeof(eigenVector));
    errorsAssert(eigenVectors != NULL,0);
    for (i = 0; i < numOfVectors; i++) { /*sets eigenvector's attributes*/
        eigenVectors[i].columnIndex = i;
        eigenVectors[i].eigenVal = eigenVals[i];
    }
    
    /*sorting*/
    qsort(eigenVectors, numOfVectors, sizeof(eigenVector), compareEigenVectors);
    for (i = 0; i < numOfVectors; i++) {
        eigenVals[i] = eigenVectors[i].eigenVal;
    }
}

int eigengapHeuristic(){
    /*calculates eigengaps for eigengap heuristic and calculates k*/
    int i, limit, maxGap, k;
    double **A;
    eigenVals = (double *)calloc(numOfVectors, sizeof(double));
    A = jacobi(laplacianNorm(), 0); /*not for printing*/
    for (i = 0; i < numOfVectors; i++) {
        eigenVals[i] = A[i][i]; /*eigenvals are on the diagonal line*/
    }
    sortEigenVectorsAndValues(); /*sorting eigenvectors and eigenvals*/
    eigenGaps = (double *)calloc(numOfVectors - 1, sizeof(double));
    for (i = 0; i < numOfVectors - 1; i++) {
        /*calculates eigen gaps*/
        eigenGaps[i] = fabs(eigenVals[i]-eigenVals[i+1]);
    }
    limit = floor(numOfVectors / 2);
    for (i = 0; i < limit; i++) { /*finds k*/
        if (eigenGaps[i] > maxGap) {
            maxGap = eigenGaps[i];
            k = i;
        }
    }
    return k + 1; /*becuase count in intructions starts from 1*/
}

void normalizeUMatrix() {
    /*normalizes U matrix to T matrix according to formula*/
    int i,j;
    double sum;
    for (i = 0; i < numOfVectors; i++){
        for (j = 0; j < k; j++){
            sum += pow(U[i][j],2);
        }
        sum = sqrt(sum);
        for (j = 0; j < k; j++){
            U[i][j] = U[i][j] / sum;
        }
    }
}

void createUMatrix() {
    /*takes k-smallest-eigenvals vectors from V matrix*/
    int i,j;
    squareMatrixTranspose(V, numOfVectors); /*in order to place vectors as rows*/

    U = (double **)calloc(numOfVectors, k*sizeof(double));
    errorsAssert(U != NULL,0);
    for (i = 0; i < numOfVectors; i++) {
        U[i] = (double *)calloc(k, sizeof(double));
        errorsAssert(U[i] != NULL,0);
        for (j = 0; j < k; j++){
            /*takes relevant part of vectors and puts it as a column*/
            U[i][j] = V[eigenVectors[j].columnIndex][i]; 
        }
    }
    normalizeUMatrix();
}

void printVectors() {
    /*prints the vectors in vectors matrix*/
    int i, j;
    for (i = 0; i < numOfVectors; i++) {
        for (j = 0; j < dimension; j++) {
            printf("%.4f", vectors[i][j]); /*format the floats precision to 4 digits*/
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
    int counter = 1;
    
    errorsAssert(argc == 4,1); /*Checks if we have the right amount of args*/ 
    
    errorsAssert(sscanf(argv[1], "%f", &rawK) == 1,1);
    k = (int)rawK;
    errorsAssert(rawK - k == 0 && k >= 0,1); /*checks if k is a non-negative int*/
    
    file = fopen(argv[3],"r");
    readFile(file);

    goal = argv[2];

    if (strcmp(goal,"spk")==0){
        int calcK = eigengapHeuristic();
        if (k==0) {
            k = calcK;
        }
        createUMatrix();
        assignUToVectors();
        initCentroids();
        clusters = (int **)calloc(k, numOfVectors*sizeof(int));
        errorsAssert(clusters != NULL,0);
        while ((counter <= max_iter) && (changes > 0)) {
            assignVectorToCluster();
            updateCentroidValue();
            counter += 1;
        }
        printResult();
    } 
    else if (strcmp(goal,"wam")==0){
        printMatrix(weightedAdjacencyMatrix(),numOfVectors,numOfVectors);
    } 
    else if (strcmp(goal,"ddg")==0){
        printMatrix(diagonalDegreeMatrix(1,1),numOfVectors,numOfVectors);
    } 
    else if (strcmp(goal,"lnorm")==0){
        printMatrix(laplacianNorm(),numOfVectors,numOfVectors);
    } 
    else if (strcmp(goal,"jacobi")==0){
        jacobi(vectors, 1);
    } 
    else{
        errorsAssert(0==1,1); /*If the goal is unknown*/
    }
        
    free(vectors);
    free(centroids);
    free(clusters);
    free(clustersSizes);
    return 0;
}