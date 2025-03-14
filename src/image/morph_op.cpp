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
#include <iostream>

using namespace std;
/*Allocate a row*column matrix initiallized with some value*/
int** createMatrix(int rows, int cols, int fill_n){
    int **matrix = new int*[rows];
    for(int i = 0; i < rows; i++){
        matrix[i] = new int[cols];
        for(int j = 0; j < cols; j++){
            matrix[i][j] = fill_n;
        }
    }
    return matrix;
}

/*Allocate a double row*column matrix initiallized with some value*/
double** createDoubleMatrix(int rows, int cols, int fill_n){
    double **matrix = new double*[rows];
    for(int i = 0; i < rows; i++){
        matrix[i] = new double[cols];
        for(int j = 0; j < cols; j++){
            matrix[i][j] = fill_n;
        }
    }
    return matrix;
}

/*Copy image into new matrix*/
int** copyImage(int **img, int rows, int cols){
    int **copy = createMatrix(rows, cols, 0);

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            copy[i][j] = img[i][j];
        }
    }
    return copy;
}

/*create a flat structuring element of specified shape and size
    allowed shapes: square, cross, disk, line, diamond
*/
int** createStructuringElement(string shape, int fill_n,int rows=0, int cols=0, int radius = 0, int angle = 0){
    int **strel;
    if(radius == 0){
        strel = createMatrix(rows, cols, 0);
    }
    else{
        strel = createMatrix(2*radius+1, 2*radius+1, 0);
        rows = 2*radius+1;
        cols = 2*radius+1;
    }

    if (shape == "square"){
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                strel[i][j] = fill_n;
            }
        }
    }
    else if (shape == "cross"){
        for(int i = 0; i < rows; i++){
            strel[i][cols/2] = fill_n;
        }
        for(int j = 0; j < cols; j++){
            strel[rows/2][j] = fill_n;
        }
    }
    else if (shape == "disk" && radius > 0){
        for(int i = 0; i < radius*2+1; i++){
            for(int j = 0; j < radius*2+1; j++){
                if((i - radius)*(i - radius) + (j - radius)*(j - radius) <= radius*radius){
                    strel[i][j] = fill_n;
                }
            }
        }
    }
    else if(shape == "line" && angle >= 0){
        if(angle == 0){
            for(int i = 0; i < rows; i++){
                strel[i][cols/2] = fill_n;
            }
        }
        else if(angle == 90){
            for(int j = 0; j < cols; j++){
                strel[rows/2][j] = fill_n;
            }
        }
        else if(angle == 135){
            for(int i = 0; i < rows; i++){
                strel[i][i] = fill_n;
            }
        }
        else if(angle == 45){
            for(int i = 0; i < rows; i++){
                strel[i][cols - (i+1)] = fill_n;
            }
        }
    }
    else if(shape == "diamond" && radius > 0){
        for(int i = 0; i < radius*2+1; i++){
            for(int j = 0; j < radius*2+1; j++){
                if(abs(i - radius) + abs(j - radius) <= radius){
                    strel[i][j] = fill_n;
                }
            }
        }
    }
    return strel;
}

/*Convoluttion operation with variable operator*/
int** convolution(int **img, int **kernel, int rows, int cols, int k_radius, int op, int** mask = nullptr){
    int **result;
    
    if(op == 1){
        //initiallize with lowest value
        result = createMatrix(rows,cols,0);
    }else{
        //initiallize with highest value
        result = createMatrix(rows,cols,255);
    }

    //image coordinates
    int x = 0, y = 0;

    //max min values for difference operator
    int max = 0, min = 255;

    //image index
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            //kernel index
            if(mask != nullptr && mask[i][j] == 0)
                continue;

            for(int k = 0; k < k_radius*2+1; k++){
                for(int l = 0; l < k_radius*2+1; l++){
                    //operate only over structuring element
                    if(kernel[k][l] != 0){
                        //calculate image-kernel window
                        x = i + k - k_radius;
                        y = j + l - k_radius;
                        //exclude out-of-boundary pixels
                        if(x >= 0 && x < rows && y >= 0 && y < cols){
                            //keep maximum value (dilation)
                            if(op == 1){
                                if(img[x][y]*kernel[k][l] > result[i][j]){
                                    result[i][j] = img[x][y]*kernel[k][l];
                                }
                            }
                            //keep minimum value (erosion)
                            else if(op == 2){
                                if(img[x][y]*kernel[k][l] < result[i][j]){
                                    result[i][j] = img[x][y]*kernel[k][l];
                                }
                            }
                            //keep difference between maximum and minimum value (gradient)
                            else if(op == 3){
                                if(img[x][y]*kernel[k][l] > max){
                                    max = img[x][y]*kernel[k][l];
                                }
                                if(img[x][y]*kernel[k][l] < min){
                                    min = img[x][y]*kernel[k][l];
                                }
                                result[i][j] = max - min;
                            }
                        }
                    }
                }
            } 
        }
    }

    return result;
}

/*Erosion morphological operation*/
int** erosion(int **img, int **kernel, int rows, int cols, int k_radius,int** mask = nullptr){
    return convolution(img, kernel, rows, cols, k_radius, 2, mask);
}

/*Dilation morphological operation*/
int** dilation(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    return convolution(img, kernel, rows, cols, k_radius, 1, mask);
}

/*opening morphological operation*/
int** opening(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    //apply erosion
    int** temp =  convolution(img, kernel, rows, cols, k_radius, 2, mask);
    //followed by dilation
    int** result = convolution(temp, kernel, rows, cols, k_radius, 1, mask);

    //free temp matrix
    delete temp[0];
    delete temp;

    return result;
}

/*closing morphological operation*/
int** closing(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    //apply dilation
    int** temp =  convolution(img, kernel, rows, cols, k_radius, 1,mask);
    //followed by erosion
    int** result = convolution(temp, kernel, rows, cols, k_radius, 2, mask);

    //free temp matrix
    delete temp[0];
    delete temp;

    return result;
}

/*gradient morphological operation*/
int** gradient(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    return convolution(img, kernel, rows, cols, k_radius, 3, mask);
}

/*top-hat morphological operation*/
int** top_hat(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    //apply opening
    int** result = opening(img,kernel,rows,cols,k_radius, mask);

    //subtract result morph operation image from original image
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            result[i][j] = img[i][j] - result[i][j];
            if(result[i][j] < 0)
                result[i][j] = 0;
        }
    }

    return result;
}

/*black-hat morphological operation*/
int** black_hat(int **img, int **kernel, int rows, int cols, int k_radius, int** mask = nullptr){
    //apply opening
    int** result = closing(img,kernel,rows,cols,k_radius, mask);

    //subtract result morph operation image from original image
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            result[i][j] -= img[i][j];
            if(result[i][j] < 0)
                result[i][j] = 0;
        }
    }

    return result;

}