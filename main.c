#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

size_t ROWS = 0, COL = 0;
int flag = 1;
double eps = 0.000001;
double M_PI = 3.14159265358979323846;

#define MAX_SIZE 100000

#define MATRIX "matrix.txt"

char* findSize(char *inputData) {

    FILE* file;
    int i=0;
    char ch;
    char buff[MAX_SIZE];

    fopen_s(& file, inputData, "r");
      
    while (fgets(buff, sizeof(buff), file)) {
        ROWS++; 
    }
   
    fclose(file);
   
    fopen_s(&file, inputData, "r");
    fgets(buff, sizeof(buff), file);

    while (buff[i] != '\0') {
        if (buff[i] == ',') {
            COL++;
        }
        i++;
    }
    COL++;

    fclose(file);

    printf("rows: %d \n", ROWS); 
    printf("cols: %d \n", COL); 

    return inputData;
}

/*Convert text file to a matrix in C*/
void readcsv(double* mat, char *inputData) {

    FILE* file;
    int i = 0;
    char *token;
    char* rest; 
    char buff[MAX_SIZE];

    fopen_s(&file, inputData, "r");
 
    while (fgets(buff, MAX_SIZE, file) != NULL) {
        token = strtok_s(buff, ",", &rest); // Get first token
        while (token != NULL) {
            mat[i] = atof(token);
            token = strtok_s(NULL, ",", &rest); // Get next token
            i++;
        }
    }
}
void writecsv(double* mat, char *inputfile, size_t ind, int type) {

    FILE* file;
    size_t i = 0, s = 0, b = 0;
    char str[1000];
    
    fopen_s(&file, inputfile, "w");
    if (type == 1) {
        b = ind * ind - ind;
        while (s < ind) {
            while (i < ind * ind) {
                sprintf_s(str, sizeof(str), "%f", mat[i]);
               
                if (i == b)
                    fprintf(file, "%s ", str);
                else
                    fprintf(file, "%s%s ", str, ", ");
                i = i + ind;
            }
            s++;
            i = s;
            b++;
            fprintf(file, "\n", "");
        }
    }
    else if (type == 2) { //segmaMatrix

        b = ind;
        for(i=0; i<b; i++){
            sprintf_s(str, sizeof(str), "%f",sqrt(fabs(mat[i])));
           
            if (i == b)
                fprintf(file, "%c", str);
            else
                fprintf(file, "%s%s ", str, ", ");
        }   
    }

    else {
        b = ind;
        while (s < ind) {
            while (i < b) {
                sprintf_s(str, sizeof(str), "%f", mat[i]);
        
                if (i == b)
                    fprintf(file, "%s", str);
                else
                    fprintf(file, "%s%s ", str, ", ");
                i++;
            }
            s++;
            b = b + ind;
            fprintf(file, "\n", "");
        }
    }

    fclose(file);
}
void sortArray(double *mat, int *arr, size_t n) {

    size_t i, j, ind;
    double x;

    for (i = 0; i < n; ++i){
        for (j = i + 1; j < n; ++j) {
            if (mat[i] < mat[j]) {
                x = mat[i];
                mat[i] = mat[j];
                mat[j] = x;

                ind = arr[i];
                arr[i] = arr[j];
                arr[j] = ind;
            }
        }
    }
}
void find_U_Matrix(double A[], double VT[], double segma[], double U[]) {

    double* tempA = (double*)calloc(ROWS * COL, sizeof(double));
    size_t j = 0;
    size_t z = 0;
    size_t g = 0;
    size_t ind1 = 1, ind2 = COL;
    size_t n = 0;
    size_t m = 0;
    size_t min = ROWS;
    size_t max = ROWS;


    if (ROWS > COL) {
        min = COL;
    }

    printf("---------------------------------------U-----------------------------------------------\n"); 
    for (size_t i = 0; i < min; i++) {
        if (segma[i] != 0) {
            double temp = 1 / sqrt(fabs(segma[i]));
            while (j < ROWS * COL) {
                tempA[j] = temp * A[j];
                j++;
            }
            j = 0;
            while (n < max) {
                while (m < min) {
                    U[z] = U[z] + tempA[j] * VT[g]; 
                    j++;
                    g++;
                    m++;
                }
                g = ind2 - COL;
                z = z + max;
                n++;
                m = 0;
            }
            z = ind1;
            g = ind2;
            ind2 = ind2 + COL;
            ind1++;
            j = 0;
            m = 0;
            n = 0;
        }
        else {
            z = ind1;
            g = ind2;
            ind2 = ind2 + COL;
            ind1++;
            j = 0;
            m = 0;
            n = 0;
        }
    }
    printf("-------------------------------------End U-------------------------------------------\n");
    printf("\n\n");
}

void printMartix(double VT[], double segma[], double U[]) {

    printf("VT temp matrix: \n");
    size_t j = 1, r = 0, stop = COL;
    size_t index = 1;

    for (int t = 0; t < COL; t++) {
        while (r < COL * COL) {
            printf("%f ", VT[r]);
            r = r + COL;
        }
        r = j;
        j++;
        stop = stop + COL;
        printf("\n");
    }

    printf("\n\n");

    r = 0;
    stop = COL;
    printf("segma Matrix: \n"); 
    for (int t = 0; t < COL; t++) {
        printf("%f  ", sqrt(fabs(segma[t])));
    }

    r = 0;
    j = 0;
    stop = ROWS;
    printf("\n\n");

    printf("U Matrix: \n"); 
    for (int t = 0; t < ROWS; t++) {
        while (j < stop) {
            printf("%f ", U[r++]);
            j++;
        }
        j = 0;
        printf("\n");
    }
}
void mat_identity(double v[]) {

    size_t i = 0;

    while (i < (COL * COL)) {

        v[i] = 1;
        i++;
        i = i + COL;
    }
    return;
}

void svd(double A[], double ATA[], double VT[], double segma[], double segmaTemp[], double U[], size_t sizeATA, int it_num, char *inputData) {

    size_t i, j, p1, p2, q1, q2, k = 0;
    double maxL_abs, maxL = -1;
    size_t b;
    size_t indI, indJ, indMax;
    size_t o;
    double* arrI1 = (double*)calloc(COL, sizeof(double));
    double theta, c, s, m, res;
    int flag = 1;

    mat_identity(VT);

    printf("The file is: %s \n", inputData);
    printf("The number of iterations: %d \n", it_num); 

    printf("\n\n");
    printf("------------------------------------VT Matrix------------------------------------------\n"); 

    for (o = 0; o < it_num; o++) {
        if (flag) {
            i = 1;
            j = 2;
            b = COL;
            maxL_abs = fabs(ATA[i]);
            maxL = ATA[i];

            while (i < COL * COL) { //find max from L
                while (i < b) {
                    if (fabs(ATA[i]) >= maxL_abs) {
                        maxL_abs = fabs(ATA[i]);
                        maxL = ATA[i];
                        indI = i / COL;
                        indJ = i % COL;
                        indMax = i;
                    }
                    i++;
                }
                i = i + j;
                j++;
                b = b + COL;
            }

            if (fabs(maxL) <= eps) {
                flag = 0;
            }

            size_t ii = indI * COL + indI;
            size_t jj = indJ * COL + indJ;

            if (fabs(ATA[ii] - ATA[jj]) < eps) {
                theta = (M_PI / 4);
            }
            else {
                m = (ATA[jj] - ATA[ii]);
                res = (2 * ATA[indMax]) / m;
                theta = atan(res) * 0.5;
            }

            c = cos(theta);
            s = sin(theta);

            size_t tempindI = indJ;
            size_t tempindJ = indI;
            size_t tempindMax = tempindI * COL + tempindJ;

            double tempI = ATA[ii];//Value in diagonal ii
            double tempJ = ATA[jj];//Value in diagonal jj
            double tempMaxIJ = ATA[indMax];//upper max
            double tempMaxJI = ATA[tempindMax];//lower max


            double m = tempI - tempJ;
            if (fabs(m) < eps) {
                m = 0;
            }

            ATA[ii] = (c * c * tempI) - (2 * s * c * tempMaxIJ) + (s * s * tempJ);
            ATA[jj] = (s * s * tempI) + (2 * s * c * tempMaxIJ) + (c * c * tempJ);
            //ATA[indMax] = ATA[tempindMax] = (c * c - s * s) * (tempMaxIJ) + (s * c * (tempI - tempJ));
            ATA[indMax] = ATA[tempindMax] = (c * c - s * s) * (tempMaxIJ)+(s * c * (m));


            p1 = COL * indI;
            p2 = indI;
            q1 = indJ;
            q2 = COL * indJ;

            size_t x = 0, y = 0;

            size_t itr = COL * (indI + 1);

            while (p1 < itr) {

                if (p1 != ii && p1 != indMax) {
                    arrI1[x++] = ATA[p1];
                    ATA[p1] = ATA[p2] = c * ATA[p1] - s * ATA[q1];
                }
                p1++;
                p2 = p2 + ROWS;
                q1 = q1 + ROWS;
            }

            p1 = COL * indI;
            p2 = indI;
            q1 = indJ;
            q2 = COL * indJ;

            x = 0;

            while (p1 < itr) {
                if (q2 != jj && q2 != tempindMax) {
                    ATA[q1] = ATA[q2] = s * arrI1[x] + c * ATA[q1];
                }
                p1++;
                x++;
                q1 = q1 + COL;
                q2++;
            }

            size_t VindI = indI;
            size_t VindJ = indJ;
            size_t z = 0;

            while (z < ROWS) {
                double vip = VT[VindI];
                double viq = VT[VindJ];
                VT[VindI] = c * vip - s * viq;
                VT[VindJ] = c * viq + s * vip;
                VindI = VindI + COL;
                VindJ = VindJ + COL;
                z++;
            }
            o++;

            k = 0;
            b = COL + 1;
            i = 0;
            size_t n = COL;

            while (k < COL * COL) {
                while (k < n) {
                    if (fabs(ATA[k]) < eps && k != b && k != 0) {
                        if (ATA[k] != 0)
                            ATA[k] = 0;
                    }
                    k++;
                }
                i++;
                b = b + COL + i;
                n = n + COL;
            }
        }
    }
    int* arrindexSort = (int*)calloc(COL, sizeof(int));
    double* tempVecsort = (double*)calloc(sizeATA, sizeof(double));

    i = 0;
    j = 0;
    while (i < COL * COL) {
        segmaTemp[j] = ATA[i];
        arrindexSort[j] = j;
        i++;
        i = i + COL;
        j++;
    }

    sortArray(segmaTemp, arrindexSort, COL);

    size_t g = 0, h=0;
    for (size_t i = 0; i < COL; i++) {
        j = arrindexSort[i];
        while (g < COL) {
            tempVecsort[h] = VT[j];
            g++;
            j = j + COL;
            h++;
        }
        g = 0;
    }

    for (size_t i = 0; i < COL*COL; i++) {
        VT[i] = tempVecsort[i];    
    }

    printf("------------------------------------end VT-----------------------------------------\n");
    printf("\n\n");

    printf("-----------------------------------Segma----------------------------------------------\n"); 

    printf("----------------------------------End Segma--------------------------------------------\n"); 
    printf("\n\n");

    find_U_Matrix(A, VT, segmaTemp, U);
}


/***********main*****************/
int main() {

    clock_t start, end;
    double cpu_time_used;
    char *inputFile;

    start = clock();//Time is running out

    inputFile=findSize(MATRIX);//

    size_t row1 = ROWS, col1 = COL;// size of matrix A
    size_t row2 = COL, col2 = ROWS;//size of matrix AT

    double* A = (double*)calloc(row1*col1, sizeof(double)); 

    readcsv(A, inputFile);//

    double* AT = (double*)calloc(row2 * col2, sizeof(double));

    size_t k = 0;
    size_t i = 0;
    size_t j = 1;
    size_t d = 0;

    while (i < ROWS * COL) {
        while (k < ROWS * COL) {
            AT[i++] = A[k];
            k = k + COL;
        }
        k = j;
        j++;
    }
    
    size_t sizeATA = row2 * col1;

    double* ATA = (double*)calloc(sizeATA, sizeof(double));

    i = 0;
    k = 0;
    j = 0;
    d = 0;

    size_t stop1 = col1, stop2 = col2;
    size_t indJ = 1, indK = 0;

    while (i < sizeATA) {
        while (i < stop1) {
            while (d < stop2) {
                ATA[i] = ATA[i] + AT[k] * A[j];
                k++;
                j = j + col1;
                d++;
            }
            i++;
            k = indK;
            j = indJ;
            indJ++;
            d = 0;

        }
        stop1 = stop1 + col1;
        indK = indK + col2;
        k = indK;
        j = 0;
        indJ = 1;
    }
 
    size_t sizeS = col1;
    if (row1 > col1)
        sizeS = row1;

    double* segma = (double*)calloc(sizeATA, sizeof(double));
    double* segmaTemp = (double*)calloc(sizeATA / col1, sizeof(double));
    double* VT = (double*)calloc(sizeATA, sizeof(double));
    double* U = (double*)calloc(row1*col2, sizeof(double));

    svd(A, ATA, VT, segma, segmaTemp, U, sizeATA, 2000, inputFile);
    //printMartix(VT, segmaTemp, U);//print result

    writecsv(VT, "VT_res.txt", COL, 1);
    writecsv(segmaTemp, "SEGMA_res.txt", COL, 2);
    writecsv(U, "U_res.txt", ROWS, 3);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("\n\nThe total running time is: %f", cpu_time_used);

    printf("\n\n");

    return 0;
}


/* r = 0;
    stop = COL;
    printf("ATA:\n");
    for (int t = 0; t < row2; t++) {
        while (r < stop) {
            printf("%f ", ATA[r++]);
        }
        stop = stop + col2;
        printf("\n");
}*/
/*  int r = 0;
    size_t stop = col2;
    printf("AT:\n"); 
    for (int t = 0; t < row2; t++) {
        while (r < stop) {
            printf("%f ", AT[r++]); 
        }
        stop = stop + col2;
       printf("\n");
}*/