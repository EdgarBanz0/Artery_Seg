/*Implementation of morphological operations
    *erosion
    *dilation
    *opening
    *closing
    *gradient
    *tophat
    *skeletonization

    Biomedical Image Processing
    Edgar Aguilera Hernández
    02/01/2025
*/

#include <string>

using namespace std;
/*Allocate a row*column matrix initiallized at value */
int** createMatrix(int rows, int cols, int fill_n);

/*Allocate a double row*column matrix initiallized at value */
double** createDoubleMatrix(int rows, int cols, int fill_n);

/*create a flat structuring element of specified shape and size
    allowed shapes: square, cross, disk, line, diamond
*/
int** createStructuringElement(string shape, int fill_n, int rows=0, int cols=0, int radius = 0, int angle = 0);

/*Copy image into new matrix*/
int** copyImage(int **img, int rows, int cols);

/*Convoluttion operation with variable operator*/
int** convolution(int **img, int **kernel, int rows, int cols, int k_radius, int op, int** mask = nullptr);

/*Erosion morphological operation*/
int** erosion(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*Dilation morphological operation*/
int** dilation(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*opening morphological operation*/
int** opening(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*closing morphological operation*/
int** closing(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*gradient morphological operation*/
int** gradient(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*top-hat morphological operation*/
int** top_hat(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);

/*black-hat morphological operation*/
int** black_hat(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr);