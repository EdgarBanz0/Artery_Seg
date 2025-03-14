/*Artery segmentation using morphological operations

    Biomedical Image Processing
    Edgar Aguilera Hernández
    02/01/2025
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <bits/stdc++.h>
#include "image/morph_op.hpp"

//number of elements in dataset
int db_size;
//first number in dataset
int db_init;
//root path for dataset
string db_path;
//path to store enhanced images
string save_path_enhance;
//path to store segmented images
string save_path_segment;
//path for masks
string mask_path;
//path for groundtruth
string gt_path;

using namespace std;
class Image{
    private:
        int** img;
        int rows;
        int cols;
        int const MAXVALUE = 255;

    public:
        Image(){
            rows = 0;
            cols = 0;
        }

        ~Image(){
            if(rows != 0){
                delete img[0];
                delete img;
            }
        }

        void setImage(int** image,int new_rows, int new_cols, bool erase_prev = true, bool copy = false){
            
            if(rows != 0 && erase_prev){
                delete img[0];
                delete img;
            }
            
            if(copy)
                img = copyImage(image, new_rows,new_cols);
            else
                img = image;
            rows = new_rows;
            cols = new_cols;
        }

        int **getImage(){
            return img;
        }

        int **getCopyImage(){
            return copyImage(img,rows,cols);
        }

        int getRows(){
            return rows;
        }

        int getCols(){
            return cols;
        }
        
        /*read a pgm file in format P2(ASCII) and P5(Binary)*/
        int pgmRead(string fileName){

            ifstream file;
            string line; /* for character input from file */
            int maximumValue = 0; /* max value from header */
            int binary;           /* flag to indicate if file is binary (P5) or ascii (P2)*/
            long numberRead = 0;  /* counter for number of pixels read */
            long i,j;             /* (i,j) for loops */
            int test,temp;        /* for detecting EOF(test) and temp storage */

            //check for previous stored image
            if(rows != 0){
                delete img[0];
                delete img;
            }

            /* Open the file, return an error if necessary. */
            file.open(fileName, ios_base::binary);
            if (!(file.is_open())) {
                printf ("ERROR: cannot open file to read\n\n");
                file.close();
                return (0);
            }
            
            /* Check the file signature ("Magic Numbers" P2 and P5); skip comments
                and blank lines (CR with no spaces before it).*/
            getline(file,line);
            while (line[0]=='#' || line[0]=='\n') 
                getline(file,line);

            if (line[0]=='P' && (line[1]=='2')) 
                binary = 0;
            
            else if (line[0]=='P' && (line[1]=='3')) 
                binary = 1;
            else if (line[0]=='P' && (line[1]=='5')) 
                binary = 1;
            
            else {
                printf ("ERROR: incorrect file format\n\n");
                file.close();
                return (0);
            } 

            /* Input the width, height and maximum value, skip comments and blank lines. */
            getline(file,line);
            while (line[0]=='#' || line[0]=='\n') 
                getline(file,line);
            sscanf (line.c_str(),"%d %d",&cols,&rows);

            getline(file,line);
            while (line[0]=='#' || line[0]=='\n') 
                getline(file,line);
            sscanf (line.c_str(),"%d",&maximumValue);

            if (cols<1 || rows<1 || maximumValue<0 || maximumValue>MAXVALUE){
                printf ("ERROR: invalid file specifications (cols/rows/max value)\n\n");
                file.close();
                return (0);
            }

            // creating memory for the input image
            img = createMatrix(rows,cols,0);

            /* Read in the data for binary (P5) or ascii (P2) PGM formats*/
            if (binary) {
                for (i = 0; i < rows; i++) {
                    for (j = 0; j < cols; j++) {
                        file.read((char*)&img[i][j],1); //read byte
                        numberRead++;
                    }
                }
            }
            else {
                for (i= 0; i < rows; i++) {
                    for (j =0; j < cols; j++) {
                        file >> temp;                    
                        img[i][j] = (int)temp;
                        numberRead++;
                    }
                    test = file.eof() ? EOF : 1;
                    if (test == EOF) break;
                }
            }

            /* close the file and return 1 indicating success */
            file.close();
            return (1);
        }



        /*create a pgm file from image matrix*/
        int pgmWrite(string fileName, string comment_string, int** image = nullptr, int i_rows = 0, int i_cols = 0)
        {
            ofstream file;
            long nwritten = 0; /* counter for the number of pixels written */
            long i;          /* for loop counters */
            
            int **img_write;
            if(image == nullptr){
                img_write = img;
            }else{
                img_write = image;
            }

            if(i_rows == 0){
                i_rows = rows;
                i_cols = cols;
            }
                

            /* open the file; write header and comments specified by the user. */
            file.open(fileName);
            if (!(file.is_open())) {
                printf ("ERROR: cannot open file to write\n\n");
                file.close();
                return (0);
            }
            
            //file headers
            file << "P2\n";
            if (comment_string.size() != 0) 
                file << "#" << comment_string <<endl;
            file << i_cols <<" "<< i_rows <<endl;//file size
            file << "255\n";//max value

            /* Write data */
            for (i = 0;i < i_rows;i++) {
                for(int j = 0; j < i_cols; j++){
                    file << img_write[i][j]<<endl;
                }
                nwritten ++;
            }
            //cout << "Image saved :" <<nwritten<<endl;
                
            file.close();
            return(1);
        }


        /*apply morphological operation
            allowed operations: erosion, dilation, opening, tophat, gradient
        */
        int **morphOp(string op, int** strel,int strel_radius,bool inplace = false, int** mask = nullptr){

            int **temp;

            if(op == "erosion"){
                temp = erosion(img, strel, rows, cols, strel_radius, mask);
            }
            else if(op == "dilation"){
                temp = dilation(img, strel, rows, cols, strel_radius, mask);
            }
            else if(op == "opening"){
                temp = opening(img, strel, rows, cols, strel_radius, mask);
            }
            else if(op == "gradient"){
                temp = gradient(img, strel, rows, cols, strel_radius, mask);
            }
            else if(op == "tophat"){
                temp = top_hat(img, strel, rows, cols, strel_radius, mask);
            }
            else if(op == "blackhat"){
                temp = black_hat(img, strel, rows, cols, strel_radius, mask);
            }
            else{
                cout<< "Operacion morfologica no valida\n";
            }

            if(inplace){
                delete img[0];
                delete img;
                img = temp;
            }

            return temp;
        }

        /*apply convolution operation with a given kernel*/
        int** convolution(int** kernel,int k_rows, int k_cols, int inplace = true){
            int **img_write;
            if(inplace){
                img_write = img;
            }else{
                img_write = createMatrix(rows,cols,0);
            }
            
            //filter matrix center
            int center_i = k_rows/2;
            int center_j = k_cols/2;

            //window mutiply accumulator and division coefficient
            int aux;
            int div_c = 0;
            
            //calculate division coeficient
            for(int i= 0; i < k_rows; i++){
                for(int j= 0; j < k_cols; j++){
                    div_c += kernel[i][j];
                }
            }
            if(div_c == 0)
                div_c = 1;
            
            //iterator over patches pixels
            for ( int y = 0; y < rows; y++ ){
                for ( int x = 0; x < cols; x++ ){
                    //reset accumulator values
                    aux = 0;
                    //filter window coordinates
                    for(int i = y-center_i; i <= y+center_i; i++){
                        for(int j = x-center_j; j <= x+center_j; j++){
                            
                            //exclude non fitting filter pixels
                            if((i<0) || (j<0) || (i>=rows) || (j>=cols)){
                                continue;
                            }
                            //multiply pixel by kernel cell
                            aux += img[i][j] * kernel[center_i-(y-i)][center_j-(x-j)];
                        }
                    }

                    img_write[y][x] = aux / abs(div_c);
                }
            }

            return img_write;
        }

        /*calculate difference between original image and img_subtract*/
        int** diffImage(int **img_subtract, bool inplace= true){
            int **img_write;
            if(inplace){
                img_write = img;
            }else{
                img_write = copyImage(img,rows,cols);
            }

            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    img_write[i][j] -= img_subtract[i][j];
                    if(img_write[i][j] < 0)
                        img_write[i][j] = 0;
                }
            }

            return img_write;
        } 

        /*add img_add values to original image*/
        void addImage(int **img_add,int **mask = nullptr){
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    if(mask != nullptr && mask[i][j] == 0)
                        continue;

                    img[i][j] += img_add[i][j];
                    if(img[i][j] > 255)
                        img[i][j] = 255;
                }
            }
        }

        /*add img_Add values where there are no previous values in the original image*/
        void addCountour(int **img_add){

            //cover left half
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    
                    if(img[i][j] < 50){
                        img[i][j] = img_add[i][j];
                        if(img[i][j] > 255)
                            img[i][j] = 255;
                    }
                    
                }
            }

        }

        /*fill outside circular given mask*/
        void fillCountour(int **mask){
            int** img_contour = createMatrix(rows,cols,0);
            //cover left half
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< (cols/2)+1; j++){

                    img_contour[i][j] = img[i][j];

                    if(mask[i][j-1] > 100)
                        break;
                }
            }

            //cover right  half
            for(int i = 0; i< rows; i++){
                for(int j = cols-1; j> cols/2; j--){

                    img_contour[i][j] = img[i][j];

                    if(mask[i][j+1] > 100)
                        break;
                }
            }
            //keep cut image
            delete img[0];
            delete img;

            img = img_contour;
        }

        /*difference from original image and img_diff values outside circular given mask*/
        void diffCountour(int **img_add,int **mask){

            //cover left half
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols/2+1; j++){
                    if(mask[i][j] != 0)
                        break;

                    img[i][j] -= img_add[i][j];
                    if(img[i][j] < 0)
                        img[i][j] = 0;
                }
            }

            //cover right  half
            for(int i = 0; i< rows; i++){
                for(int j = cols-1; j> cols/2; j--){
                    if(mask[i][j] != 0)
                        break;

                    img[i][j] -= img_add[i][j];
                    if(img[i][j] < 0)
                        img[i][j] = 0;
                }
            }
        }

        /*invert image*/
        void invertImage(int** mask = nullptr){
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    if(mask != nullptr && mask[i][j] == 0)
                        img[i][j] = 0;
                    else
                        img[i][j] = 255 - img[i][j];
                }
            }
        }

        /*Normalize int image into 0-250 values*/
        int** normalize(int** image = nullptr, int** mask = nullptr){
            int **img_input;
            if(image == nullptr){
                img_input = img;
            }else{
                img_input = image;
            }

            int** img_write = createMatrix(rows,cols,0);

            int max = 0;
            int min = 100000;

            //get max value
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    if(mask != nullptr){
                        if(mask[i][j] == 0)
                            continue;
                    }
                    if(max < img_input[i][j])
                        max = img_input[i][j];
                    if(min > img_input[i][j])
                        min = img_input[i][j];

                }
            }

            //apply normalization
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    if(mask != nullptr){
                        if(mask[i][j] == 0)
                            continue;
                    }
                    img_write[i][j] = round((img_input[i][j] - min)*(255)/ (max-min));
                }
            }

            //change un-normalized image
            delete img_input[0];
            delete img_input;

            if(image == nullptr){
                img = img_write;
            }

            return img_write;
        }

        /*Normalize double image into range of values, returning the int rounded version*/
        int** normalizeDouble(double** image_d, int rows, int cols, int min_range, int max_range){
            int **img_write = createMatrix(rows,cols,0);

            double max = 0;
            double min = 100000;

            //get max value
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    if(max < image_d[i][j])
                        max = image_d[i][j];
                    if(min > image_d[i][j])
                        min = image_d[i][j];

                }
            }
    
            //apply normalization
            for(int i = 0; i< rows; i++){
                for(int j = 0; j< cols; j++){
                    img_write[i][j] = round((image_d[i][j] - min)*(max_range-min_range)/ (max-min) +min_range);
                }
            }
            
            return img_write;
        }

        /*Compute mean from whole image*/
        int mean(int** image = nullptr, int** mask = nullptr){

            int **img_input;
            if(image == nullptr){
                img_input = img;
            }else{
                img_input = image;
            }
            
            //pixel accumulator
            int aux = 0;
            int n_pixels = 0;

            //iterator over patches pixels
            for ( int y = 0; y < rows; y++ ){
                for ( int x = 0; x < cols; x++ ){
                    if(mask != nullptr){
                        if(mask[y][x] != 0){
                            n_pixels++;
                            aux += img_input[y][x];
                        }
                            
                    }else{
                        n_pixels++;
                        aux += img_input[y][x];
                    }
                        
                    
                }
            }

            return round(aux / n_pixels);
        }

        /*find max value coordinates */
        void maxCoordinates(int yx_center[]){
            int max = 0;

            //iterator over image
            for ( int i = 0; i < rows; i++ ){
                for ( int j = 0; j < cols; j++ ){

                    if(img[i][j] > max){
                        max = img[i][j];
                        yx_center[0] = i;
                        yx_center[1] = j;
                    }
                        
                    
                }
            }

        }

        /*Compute mean from a 5x5 window*/
        void meanWindow(int** image = nullptr){
            int **img_write;
            int **img_input;
            if(image == nullptr){
                img_write = img;
            }else{
                img_write = image;
            }
            img_input = copyImage(img_write,rows,cols);

            int k_len = 5;

            //Mean window
            int kernel[5][5] = {{1,1,1,1,1},
                                {1,1,1,1,1},
                                {1,1,1,1,1},
                                {1,1,1,1,1},
                                {1,1,1,1,1}};
            
            
            //filter matrix center
            int center_i = k_len/2;
            int center_j = k_len/2;

            //window accumulator
            int aux;

            //iterator over patches pixels
            for ( int y = 0; y < rows; y++ ){
                for ( int x = 0; x < cols; x++ ){
                    //reset accumulator values
                    aux = 0;
                    //filter window coordinates
                    for(int i = y-center_i; i <= y+center_i; i++){
                        for(int j = x-center_j; j <= x+center_j; j++){
                            
                            //exclude non fitting filter pixels
                            if((i<0) || (j<0) || (i>=rows) || (j>=cols)){
                                continue;
                            }
                            //multiply pixel by kernel cell
                            aux += img_input[i][j] * kernel[center_i-(y-i)][center_j-(x-j)];
                        }
                    }

                    img_write[y][x] = aux / (k_len*k_len);
                }
            }
            delete img_input[0];
            delete img_input;
        }

        /*Apply gaussian filter to smooth image*/
        int** gauss_filter(bool inplace = false){
            int k_len = 5;

            //Gauss Kernel window
            int kernel[5][5] = {{1,4,7,4,1},
                                {4,16,26,16,4},
                                {7,26,41,26,7},
                                {4,16,26,16,4},
                                {1,4,7,4,1}};
            
            
            //filter matrix center
            int center_i = k_len/2;
            int center_j = k_len/2;

            //window mutiply accumulator and division coefficient
            int aux;
            int div_c = 0;

            // allocate memory if new image is expected
            int **img_write = createMatrix(rows,cols,0);

            //calculate division coeficient
            for(int i= 0; i < k_len*k_len; i++){
                div_c += *(kernel[0] + i);
            }

            //iterator over patches pixels
            for ( int y = 0; y < rows; y++ ){
                for ( int x = 0; x < cols; x++ ){
                    //reset accumulator values
                    aux = 0;
                    //filter window coordinates
                    for(int i = y-center_i; i <= y+center_i; i++){
                        for(int j = x-center_j; j <= x+center_j; j++){
                            
                            //exclude non fitting filter pixels
                            if((i<0) || (j<0) || (i>=rows) || (j>=cols)){
                                break;
                            }
                            //multiply pixel by kernel cell
                            aux += img[i][j] * kernel[center_i-(y-i)][center_j-(x-j)];
                        }
                    }

                    img_write[y][x] = aux / div_c;
                }
            }

            if(inplace){
                delete img[0];
                delete img;
                img = img_write;
            }
            
            return  img_write;
        }

        /*Apply Sobel filter for border detection */
        int **scharr_gradient(bool inplace = false, int** angle_gradient = nullptr){
            //Scharr kernel
            int Gx[3][3] = { {3, 0, -3}, {10, 0, -10}, {3, 0, -3} };
            int Gy[3][3] = { {3, 10, 3}, {0, 0, 0}, {-3, -10, -3} };
            //kernel center
            int center_i = 1;
            int center_j = 1;

            //window mutiply accumulator
            int aux_gx, aux_gy;

            // allocate memory
            int **img_write = createMatrix(rows,cols,0);

            //center pixel of image
            for(int y = 0; y < rows; y++){
                for(int x = 0; x < cols; x++){
                    aux_gx = 0;
                    aux_gy = 0;

                    //kernel coordinates
                    for(int i= y-center_i; i <= y+center_i; i++){
                        for(int j = x-center_j; j <= x+center_j; j++){
                            //exclude non fitting filter pixels
                            if((i<0) || (j<0) || (i>=rows) || (j>=cols)){
                                break;
                            }
                            
                            //take any of the color channel value
                            aux_gx += img[i][j] * Gx[center_i-(y-i)][center_j-(x-j)];
                            aux_gy += img[i][j] * Gy[center_i-(y-i)][center_j-(x-j)];
                        }
                    }
                    
                    img_write[y][x] = (int)(sqrt(aux_gx*aux_gx + aux_gy*aux_gy));
                    
                    //compute angle in degrees if matrix is given as parameter
                    if(angle_gradient != nullptr){
                        if(aux_gx == 0){
                            angle_gradient[y][x] = 0;
                        }else{
                            angle_gradient[y][x] = round(atan(aux_gy/aux_gx)*180/M_PI);
                        }
                    }
                        
                }
            }

            if(inplace){
                delete img[0];
                delete img;
                img = img_write;
            }

            return img_write;
        }

        /*compute edge using Canny algorithm*/
        void cannyEdge(){

            //1. smooth image
            gauss_filter(true);

            //2. compute gradient magnitude and direction matrix
            int** grad_angle = createMatrix(rows,cols,0);
            scharr_gradient(true,grad_angle);  
            
            //3. Non-maximum supression
            for(int i = 1; i< rows-1; i++){
                for(int j = 1; j< cols-1; j++){
                    //check for direction of gradient
                    if(abs(grad_angle[i][j]) >= 68){
                        if(img[i][j] < img[i-1][j] || img[i][j] < img[i+1][j])
                            img[i][j] = 0;

                    }else if(grad_angle[i][j] < 68 && grad_angle[i][j] >= 22){
                        if(img[i][j] < img[i-1][j+1] || img[i][j] < img[i+1][j-1])
                            img[i][j] = 0;

                    }else if(abs(grad_angle[i][j]) < 22){
                        if(img[i][j] < img[i][j+1] || img[i][j] < img[i][j-1])
                            img[i][j] = 0;

                    }else if(grad_angle[i][j] <= -22 && grad_angle[i][j] > -68){
                        if(img[i][j] < img[i-1][j-1] || img[i][j] < img[i+1][j+1] )
                            img[i][j] = 0;
                    }
                    
                    
                }
            }
            
            //set gradient into 0 - 255 level
            //normalize();

            //4. Double thrheshold filtering
            int max_t = 2000;
            int min_t = 1900;

            //weak pixel matrix
            int** weak_edges = copyImage(img,rows,cols);

            for(int i = 1; i< rows-1; i++){
                for(int j = 1; j< cols-1; j++){
                    //strong gradient
                    if(img[i][j] >= max_t){
                        weak_edges[i][j] = 0;

                    }else if(img[i][j] >= min_t){
                        img[i][j] = 0; //weak gradient

                    } else if(img[i][j] < min_t){
                        img[i][j] = 0;       //supressing gradients
                        weak_edges[i][j] = 0;
                    }   
                }
            }

            //5. Evaluate strong, weak pixels by hysteresis
            for(int i = 1; i< rows-1; i++){
                for(int j = 1; j< cols-1; j++){
                    //check for weak gradient
                    if(weak_edges[i][j] == 0)
                        continue;

                    //window around weak pixel
                    for(int k = i - 1; k <= i+1; k++){
                        for(int l = j - 1; l <= j+1; l++){                         
                            if(img[k][l] != 0){
                                img[i][j] = weak_edges[i][j];
                                break;
                            }   
                        }
                    }

                }
            }
            
            
            delete grad_angle[0];
            delete grad_angle;
            delete weak_edges[0];
            delete weak_edges;
            
        }

        /*Measure distance from a skeletonized image to the edge of its structure in a given image_matrix*/
        int** radialEdgeSearch(int** edge_image, int** mask, int radial_t, int scale = 1){

            int radius;
            int** radii_matrix = createMatrix(rows,cols,0);

            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    if(mask == 0 || img[i][j] == 0)
                        continue;

                    radius = 1;
                    while(radius <= radial_t ){

                        //window around central skeleton
                        for(int k = i - radius; k <= i+radius; k++){
                            for(int l = j - radius; l <= j+radius; l++){  
                                //exclude non fitting or pixel width found
                                if((k<0) || (l<0) || (k >= rows) || (l >= cols) || (radii_matrix[i][j] != 0))
                                    continue;
                                
                                if((k == i && l == j) && edge_image[k][l] > 0)
                                    radii_matrix[i][j] = -5;  //decimal part representation
                                else if(edge_image[k][l] > 0)
                                    radii_matrix[i][j] = radius*scale;
                                
                            }
                        }

                        //edge reached
                        if(radii_matrix[i][j] != 0)
                            break;

                        radius++;
                    }

                }
            }

            return radii_matrix;

        }

        /*get original image values that result from applying a threshold; foreground = true (if threshold < img), foreground = false (if threshold > img)*/
        int** getImageFromMask(int threshold, bool foreground = true){
            int **img_output = createMatrix(rows,cols,0);

            //iterator over original image
            for (int y = 0; y < rows; y++ ){
                for ( int x = 0; x < cols; x++ ){
                    if(foreground && img[y][x] > threshold)
                        img_output[y][x] = img[y][x];
                    if(!foreground && img[y][x] < threshold)
                        img_output[y][x] = img[y][x];
                }
            }
            return img_output;
        }

        void printImage(){
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    cout << img[i][j] << " ";
                }
                cout << endl;
            }
        }

        void printStrel(int** strel,int strel_radius){
            for(int i = 0; i < strel_radius*2+1; i++){
                for(int j = 0; j < strel_radius*2+1; j++){
                    cout << strel[i][j] << " ";
                }
                cout << endl;
            }
        }
};

class ROC{
    private:
        vector<Image*> elements;
        vector<Image*> groundtruth;
        vector<Image*> mask;
        string strel;
        double area;
        int strel_param[5];
        double confusion[4][256];   //TP, TN, FP, FN x 255 tresholds
        double sens[256];
        double spec[256];
    
    public:
        ROC(){
            strel = "none";
            area = 0;

        }

        ROC(string struct_element,int newstrel_param[]){

            strel = struct_element;
            strel_param[0] = newstrel_param[0];
            strel_param[1] = newstrel_param[1];
            strel_param[2] = newstrel_param[2];
            strel_param[3] = newstrel_param[3];
            strel_param[4] = newstrel_param[4];

            for(int i = 0; i < 256; i++){
                confusion[0][i] = 0;
                confusion[1][i] = 0;
                confusion[2][i] = 0;
                confusion[3][i] = 0;
            }
        }

        ~ROC(){
            //free allocated memory for every enhanced and groundthruth image
            for(int i= 0; i<(int)elements.size(); i++){
                delete elements[i];
            }

            if((int)groundtruth.size() > 0){
                for(int i= 0; i<(int)groundtruth.size(); i++){
                    delete groundtruth[i];
                }
            }
            if((int)mask.size() > 0){
                for(int i= 0; i<(int)mask.size(); i++){
                    delete mask[i];
                }
            }
        }

        /*read mask images from dataset into an array*/
        void buildMaskArray(string mask_path, int n_images, int init_ref=1){

            //check if mask array is already filled
            if((int)mask.size() == n_images){
                return;
            }

            Image *ptr;
            for(int i = init_ref; i < init_ref + n_images; i++){
                ptr = new Image();
                ptr->pgmRead(mask_path + to_string(i) +"_training_mask.pgm");
                // add image object to vector
                mask.push_back(ptr);
            }

            //cout<<"--------------Masks procesadas: "<< mask.size()<<endl;
        }

        /*read groundtruth images from dataset into an array*/
        void buildGroundthruthArray(string gt_path, int n_images, int init_ref=1){

            //check if groundtruth array is already filled
            if((int)groundtruth.size() == n_images){
                return;
            }

            Image *ptr;
            for(int i = init_ref; i < init_ref + n_images; i++){
                ptr = new Image();
                ptr->pgmRead(gt_path + to_string(i) +"_manual1.pgm");
                //ptr->pgmRead(gt_path + to_string(i) +"_gt.pgm");

                // add image object to vector
                groundtruth.push_back(ptr);
            }

            //cout<<"--------------Grounthruth procesadas: "<< groundtruth.size()<<endl;

        }

        /*read images from dataset into an array*/
        void buildImageArray(string dataset_path, int n_images, int init_ref, string new_strel = "test"){
            //update datset strel
            strel = new_strel;

            Image *ptr;
            for(int i = init_ref; i < init_ref + n_images; i++){
                ptr = new Image();
                cout<<dataset_path + to_string(i) +"_enhance.pgm"<<endl;
                ptr->pgmRead(dataset_path + to_string(i) +"_enhance.pgm");
                // add image object to vector
                elements.push_back(ptr);
            }

            cout<<"--------------Imagenes cargadas: "<< elements.size()<<endl;
        }

        /*insert image into array*/
        void insertEnhanceImage(Image *newImage){
            if((int)elements.size() > db_size){
                cout<<"Enhanced images array is full"<<endl;
            }
            elements.push_back(newImage);
        }

        /*reset enhanced image array*/
        void clearEnhanceArray(){
            if((int)elements.size() > 0 ){
                for(int i= 0; i<(int)elements.size(); i++){
                    delete elements[i];
                }
                elements.clear();
            }
        }

        /*calculate confusion matrix*/
        void calculateConfusionMatrix(bool available_mask = false){
            int **image;
            int **gt;
            int **mk;

            for(int i = 0; i < 256; i++){
                confusion[0][i] = 0;
                confusion[1][i] = 0;
                confusion[2][i] = 0;
                confusion[3][i] = 0;
            }

            //evalute per pixel (m x n x i_images)
            for(int i = 0; i < (int)elements.size(); i++){
                image = elements[i]->getImage();
                gt = groundtruth[i]->getImage();
                if(available_mask)
                    mk = mask[i]->getImage();  

                //apply calculation per threshold
                //***using inverted image as we need vessels appear brighter
                for(int j = 0; j< 256; j++){

                    for(int m = 0; m < elements[i]->getRows(); m++){
                        for(int n = 0; n < elements[i]->getCols(); n++){
                            if(!available_mask || mk[m][n] != 0){
                                //TP
                                if(image[m][n] >= j && gt[m][n] != 0){
                                    confusion[0][j] += 1; 
                                }
                                //TN
                                if(image[m][n] < j && gt[m][n] == 0){
                                    confusion[1][j] += 1; 
                                }
                                //FP
                                if(image[m][n] >= j && gt[m][n] == 0){
                                    confusion[2][j] += 1; 
                                }
                                //FN
                                if(image[m][n] < j && gt[m][n] != 0){
                                    confusion[3][j] += 1; 
                                }
                            }
                        }
                    }
                }
                
            }
        }

        /*Calculate Sensitivity and Specificity*/
        void calculateSensSpec(){
            for(int i = 0; i < 256; i++){
                //TP/(TP + FN)
                sens[i] = confusion[0][i] / (confusion[0][i] + confusion[3][i]);
                //TN/(TP + FP)
                spec[i] = confusion[1][i] / (confusion[1][i] + confusion[2][i]);
            }

            //sort ascending by sensitivity with selection algorithm
            double aux1,aux2;
            for(int i= 0; i < 256; i++){
                for(int j = i; j < 256; j++){
                    if(sens[j] < sens[i]){
                        aux1 = sens[j];
                        sens[j] = sens[i];
                        sens[i] = aux1;
                        aux2 = spec[j];
                        spec[j] = spec[i];
                        spec[i] = aux2;
                    }
                    
                }
            }

        }

        /*Calculate area under the curve, using Riemman sum with mid value*/
        void calculateAUC(){
            area = 0;
            //step size
            for(int i = 1; i < 256; i++){
                area += abs((1-spec[i]) - (1-spec[i-1])) * ((sens[i]+sens[i-1])/2);
            }
        }

        /*write file with set of image,groundthruth and mask*/
        void saveImageSet(int image_id, int init_ref){
            elements[image_id-init_ref]->pgmWrite(to_string(image_id) + "_writeTest.pgm","Prueba");
            groundtruth[image_id-init_ref]->pgmWrite(to_string(image_id) + "_GT_writeTest.pgm","Prueba");
            mask[image_id-init_ref]->pgmWrite(to_string(image_id) + "_Mask_writeTest.pgm","Prueba");
        }

        /*get one of the images lodaded at given index (type 1: image, 2: mask, 3: groundtruth)*/
        Image* getImage(int index,int type){
            if(type == 1)
                return elements[index];
            if(type == 2)
                return mask[index];
            else
                return groundtruth[index];
        }

        double getArea(){
            return area;
        }
 
        int getArraySize(){
            return elements.size();
        }

        void printROCData(){
            cout << setw(20) << left << strel << setw(10) << left << strel_param[3] << setw(10) << area << endl;
        }
};

//------------------------------------------Segmentation method based on adaptive local threshold
class Segment{
    private:
        vector<Image*> image;
        vector<Image*> groundtruth;
        vector<Image*> mask;
        vector<Image*> segmented;
        double confusion[4];   //TP, TN, FP, FN 

    public:

    Segment(){
        confusion[0] = 0;
        confusion[1] = 0;
        confusion[2] = 0;
        confusion[3] = 0;
        
    }

    ~Segment(){
        //free allocated memory for every enhanced and groundthruth image
        for(int i= 0; i<(int)image.size(); i++){
            delete image[i];
        }

        /*if((int)groundtruth.size() > 0){
            for(int i= 0; i<(int)groundtruth.size(); i++){
                delete groundtruth[i];
            }
        }
        if((int)mask.size() > 0){
            for(int i= 0; i<(int)mask.size(); i++){
                delete mask[i];
            }
        }
        if((int)segmented.size() > 0){
            for(int i= 0; i<(int)segmented.size(); i++){
                delete segmented[i];
            }
        }*/
    }

    /*read mask images from dataset into an array*/
    void buildMaskArray(string mask_path, int n_images, int init_ref=1){

        //check if mask array is already filled
        if((int)mask.size() == n_images){
            return;
        }

        Image *ptr;
        for(int i = init_ref; i < init_ref + n_images; i++){
            ptr = new Image();
            ptr->pgmRead(mask_path + to_string(i) +"_training_mask.pgm");
            // add image object to vector
            mask.push_back(ptr);
        }

        cout<<"--------------Masks leidas: "<< mask.size()<<endl;
    }

    /*read groundtruth images from dataset into an array*/
    void buildGroundthruthArray(string gt_path, int n_images, int init_ref=1){

        //check if mask array is already filled
        if((int)groundtruth.size() == n_images){
            return;
        }

        Image *ptr;
        for(int i = init_ref; i < init_ref + n_images; i++){
            ptr = new Image();
            ptr->pgmRead(gt_path + to_string(i) +"_manual1.pgm");
            // add image object to vector
            groundtruth.push_back(ptr);
        }

        cout<<"--------------Grountruth leidas: "<< groundtruth.size()<<endl;

    }

    /*read images from dataset into an array; array_type = 1 for image array, array_type = 2 for segmented array*/
    void buildImageArray(string dataset_path, int n_images, int init_ref=1, int array_type = 1, int img_type = 1){

        //type of image to load
        string img_appendix;
        if(img_type == 1){
            img_appendix = "_enhance.pgm";
        }
        if(img_type == 2){
            img_appendix = "_segmented.pgm";
        }

        Image *ptr;
        for(int i = init_ref; i < init_ref + n_images; i++){
            ptr = new Image();
            ptr->pgmRead(dataset_path + to_string(i) + img_appendix);
            // add image object to vector
            if(array_type == 1)
                image.push_back(ptr);
            else
                segmented.push_back(ptr);
        }

        // add image object to vector
        if(array_type == 1)
            cout<<"--------------Imagenes cargadas: "<< image.size()<<endl;
        else
            cout<<"--------------Imagenes cargadas: "<< segmented.size()<<endl;
    }

    /*Segment a series of images using Yanowitz thresholding surface method*/
    void yanowitz_method(int n_images, string save_path,int maxima_t, int connected_thresh){
        
        //--------------------------------------------------Image segmentation workflow
        Image *ptr;
        bool VisB;
        cout<<"Vessel is Black? : ";
        cin>>VisB;

        cout<<"Segmentando...";
        for(int i = 0; i< n_images; i++){
            int** original_img = image[i]->getCopyImage();
            int** img_mask = mask[i]->getCopyImage();
            int rows = image[i]->getRows();
            int cols = image[i]->getCols();
            
            //1. smooth
            image[i]->gauss_filter(true);

            //2. gradient
            image[i]->scharr_gradient(true);
            image[i]->normalize(nullptr,mask[i]->getImage());

            //3. local maxima
            
            int** max_grad_mask = localMaxima(image[i]->getImage(),rows,cols,20,maxima_t);
            image[i]->pgmWrite(save_path + to_string(i+db_init)+"_gradient_max.pgm","max local gradient image",max_grad_mask);

            //4. Get original gray levels on local maxima
            int** eval_max_mask = evaluateMaxima(original_img,max_grad_mask,rows,cols);
            image[i]->pgmWrite(save_path + to_string(i+db_init)+"_potential.pgm","potential threshold points",eval_max_mask);

            //5. interpolate with SOR over laplace derivative
            int** threshold_surface = interpolatePoints(eval_max_mask,rows,cols,1.5,3,1000);
            image[i]->pgmWrite(save_path + to_string(i+db_init)+"_thresh_surf.pgm","threshold surface",threshold_surface);

            //6. Apply threshold surface 
            int** segmented_img = segmentImage(original_img,threshold_surface,img_mask,rows,cols,VisB);

            //7. Apply connected elements algorithm to keep objects > threshold
            connected_BFS(segmented_img,rows,cols,connected_thresh);
            cout<<save_path + to_string(i)+"_segmented.pgm"<<endl;
            image[i]->pgmWrite(save_path + to_string(i+db_init)+"_segmented.pgm","post processed image segmented with Yanowitz threshold surface",segmented_img,rows,cols);

            //store segmented image
            ptr = new Image();
            ptr->setImage(segmented_img,rows,cols);
            // add image object to vector
            segmented.push_back(ptr);
            
            //clear allocated images
            delete max_grad_mask[0];
            delete max_grad_mask;
            delete eval_max_mask[0];
            delete eval_max_mask;
            delete threshold_surface[0];
            delete threshold_surface;
            delete original_img[0];
            delete original_img;
            delete img_mask[0];
            delete img_mask;

        }
        cout<<"\nProceso finalizado\n";
    }

    /*Segment a series of images using iterative thresholding method*/
    void iterative_method(int n_images, string save_path, int c_thresh){
        
        //--------------------------------------------------Image segmentation workflow
        Image *ptr;
        int** img_background;
        int** img_foreground;
        int threshold = 0,threshold_new;
        int rows = image[0]->getRows();
        int cols = image[0]->getCols();

        bool VisB;
        cout<<"Vessel is Black? : ";
        cin>>VisB;

        cout<<"Segmentando...";

        for(int i = 0; i< n_images; i++){
            //normalize image
            image[i]->normalize(nullptr, mask[i]->getImage());
            //1. initial threshold from image mean
            threshold_new = image[i]->mean(nullptr, mask[i]->getImage());
            
            while(abs(threshold - threshold_new) > 1){
                //update threshold
                threshold = threshold_new;

                //2. apply threshold to get background and foreground
                img_foreground = image[i]->getImageFromMask(threshold,true);
                image[i]->pgmWrite(save_path + to_string(i+db_init)+"_foreground.pgm","image segmented with iterative threshold method",img_foreground,rows,cols);

                img_background = image[i]->getImageFromMask(threshold,false);
                image[i]->pgmWrite(save_path + to_string(i+db_init)+"_background.pgm","image segmented with iterative threshold method",img_foreground,rows,cols);

                //3. Compute new threshold from mean of foreground and background images
                threshold_new = (image[i]->mean(img_foreground, mask[i]->getImage()) + image[i]->mean(img_background, mask[i]->getImage())) / 2;

                delete img_background[i];
                delete img_background;
                delete img_foreground[i];
                delete img_foreground;
            }
            //segment using best estimated threshold
            img_foreground = segmentImage(image[i]->getImage(),threshold_new,mask[i]->getImage(),rows,cols,VisB);


            //7. Apply connected elements algorithm
            connected_BFS(img_foreground,rows,cols,c_thresh);
            cout<<save_path + to_string(i)+"_segmented.pgm"<<endl;
            image[i]->pgmWrite(save_path + to_string(i+db_init)+"_segmented.pgm","image segmented with iterative threshold method",img_foreground,rows,cols);

            //store segmented image
            ptr = new Image();
            ptr->setImage(img_foreground,rows,cols);
            // add image object to vector
            segmented.push_back(ptr);
        

        }
        cout<<"\nProceso finalizado\n";
    }

    /*calculate confusion matrix*/
    void calculateConfusionMatrix(){
        int **img;
        int **gt;
        int **mk;
        
        //evalute per pixel (m x n x i_images)
        for(int i = 0; i < (int)segmented.size(); i++){
            img = segmented[i]->getImage();
            gt = groundtruth[i]->getImage();
            mk = mask[i]->getImage();

            
            for(int m = 0; m < segmented[i]->getRows(); m++){
                for(int n = 0; n < segmented[i]->getCols(); n++){
                    
                    if(mk[m][n] != 0){
                        //TP
                        if((img[m][n] == gt[m][n]) && (gt[m][n] != 0)){
                            confusion[0] += 1; 
                        }
                        //TN
                        if((img[m][n] == gt[m][n]) && (gt[m][n] == 0)){
                            confusion[1] += 1; 
                        }
                        //FP
                        if((img[m][n] != gt[m][n]) && (gt[m][n] == 0)){
                            confusion[2] += 1; 
                        }
                        //FN
                        if((img[m][n] != gt[m][n]) && (gt[m][n] != 0)){
                            confusion[3] += 1; 
                        }
                    }
                }
            }
        }
    }

    /*Create binary mask of local maxima over graddient image*/
    int** localMaxima(int** grad_image,int rows, int cols, int window_size, int threshold){
        int** grad_max = createMatrix(rows,cols,0);
        int max_yx[2] = {0,0};

        //iterator over window 
        for ( int y = 0; y < rows; y++ ){
            for ( int x = 0; x < cols; x++ ){
                //reset max values
                max_yx[0] = y;
                max_yx[1] = x;
                //filter window coordinates
                for(int i = y; i <= y + window_size; i++){
                    for(int j = x; j <= x + window_size; j++){
                        //exclude non fitting filter pixels
                        if((i>=rows) || (j>=cols)){
                            break;
                        }
                        //evaluate max value in local window
                        if(grad_image[max_yx[0]][max_yx[1]] < grad_image[i][j]){
                            max_yx[0] = i;
                            max_yx[1] = j;
                        }
                    }
                }

                if(grad_image[max_yx[0]][max_yx[1]]  >= threshold)
                    grad_max[max_yx[0]][max_yx[1]] = 255;
            }
        }

        return grad_max;
    }

    /*apply gradient local maxima mask over original image to map gary levels into a new image*/
    int** evaluateMaxima(int** image,int** mask, int rows, int cols){
        int** eval_image = createMatrix(rows,cols,0);

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                if(mask[i][j] > 0)
                    eval_image[i][j] = image[i][j];
            }
        }
        return eval_image;
    }

    /*Interpolate potential surfaces solving laplacian local derivatives with SOR method*/
    int** interpolatePoints(int **img,int rows, int cols,int beta, int eps, int iter_max){
        int** laplace_grad = copyImage(img,rows,cols);
        int** surface_thresh = copyImage(img,rows,cols);

        int iterations = 1;
        int max_residual = 125;

        //SOR method to reduce residuals
        while ((iterations < iter_max) && (max_residual > eps)){
            max_residual = 0;
            //residuals surface to minimize
            for(int i = 1; i < rows-1; i++){
                for(int j = 1; j < cols-1; j++){
                    laplace_grad[i][j] = surface_thresh[i][j+1] + surface_thresh[i][j-1] +surface_thresh[i+1][j] + surface_thresh[i-1][j+1] - 4*surface_thresh[i][j];
                }
            }

            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    //keep interpolation points unchange
                    if(img[i][j] ==  0){
                        surface_thresh[i][j] += (beta * laplace_grad[i][j]) / 4;
                    
                        //keep max residual for convergence indicator
                        if(max_residual < laplace_grad[i][j])
                            max_residual = laplace_grad[i][j];
                    }
                        
                }
            }
            iterations++;
        }
        return surface_thresh;
    }

    /*apply threshold surface to segment image*/
    int** segmentImage(int** img,int** thresh_surface,int** mask,  int rows, int cols, bool VisB = true){
        int** segmented_img = createMatrix(rows,cols,0);

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                if((VisB && img[i][j] <= thresh_surface[i][j]) && mask[i][j] != 0)
                    segmented_img[i][j] = 255;  
                
                if((!VisB && img[i][j] >= thresh_surface[i][j]) && mask[i][j] != 0)
                    segmented_img[i][j] = 255; 
            }
        }

        return segmented_img;
    }

    /*apply threshold to segment image*/
    int** segmentImage(int** img, int threshold, int** mask, int rows, int cols, bool VisB = true){
        int** segmented_img = createMatrix(rows,cols,0);

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                if((VisB && img[i][j] <= threshold) && mask[i][j] != 0){
                        segmented_img[i][j] = 255; 
                }
                if((!VisB && img[i][j] >= threshold) && mask[i][j] != 0){
                    segmented_img[i][j] = 255; 
                }   
                     
            }
        }

        return segmented_img;
    }

    /*Apply Breadth first search for connected elements detection*/
    void connected_BFS(int** img, int rows, int cols, int size_threshold){
        //matrix of labeled pixels
        int** connected_map = createMatrix(rows,cols,0);
        //element label id
        int label = 1;
        //array of connected elements
        queue<vector<int>> queue_connect;
        //count of connected elements
        vector<int> vector_connect;
        int c_x,c_y;
        
        //labeling pixels
        for ( int y = 0; y < rows; y++ ){
            for ( int x = 0; x < cols; x++ ){
                //check for already visited pixels or background
                if(connected_map[y][x] != 0 || img[y][x]== 0){
                    continue;
                }

                //enqueue unvisited pixel
                queue_connect.push({y,x});

                //set new label
                connected_map[y][x] = label;
                //store new label in vector count
                vector_connect.push_back(1);
                label++;

                //evaluate neighborhood from queue
                while(!queue_connect.empty()){
                    c_y = queue_connect.front()[0];
                    c_x = queue_connect.front()[1];
                    queue_connect.pop();

                    for(int i= c_y-1; i <= c_y+1; i++){
                        for(int j = c_x-1; j <= c_x+1; j++){
                            //exclude non fitting filter pixels, already visited or background
                            if((i<0) || (j<0) || (i>=rows) || (j>=cols) || (connected_map[i][j]!=0) || (img[i][j]==0)){
                                continue;
                            }
                            
                            //set same label as center
                            connected_map[i][j] = connected_map[c_y][c_x];
                            //add neighbor to queue
                            queue_connect.push({i,j});
                            vector_connect[label-2]++; 
                        }
                    }
                }
            }
        }

        //keep connected elements with length > threshold
        for( int y = 0; y < rows; y++){
            for ( int x = 0; x < cols; x++){
                if(connected_map[y][x] != 0){
                    if((vector_connect[connected_map[y][x]-1] >=  size_threshold))
                        img[y][x] = 255;
                    else
                        img[y][x] = 0;
                }

            }
        }

        delete connected_map[0];
        delete connected_map;
    }

    /*Thin white elements from binary images using Zhang-Suen algorithm*/
    void zhangSuenSkeletonization(){
        
        string save_path_skeleton = "src/db_coronary/skeletonized/";
        bool change_flag;   //changes made during any of both phases
        bool BtoW_change;   //transition from black too white
        int** image_ptr;
        int count_b,count_a;

        for(int n = 0; n < (int)image.size(); n++){
            image_ptr = image[n]->getImage();
            change_flag = true;

            while(change_flag){

                change_flag = false;
                /* order of window elements:
                    9 2 3
                    8 1 4
                    7 6 5
                */
               //phase 1
                for(int i = 1; i < image[n]->getRows() - 1; i++){
                    for(int j = 1; j < image[n]->getCols() - 1; j++){
                        if(image_ptr[i][j] == 255){
                            
                            int window[9][2] = {{i-1,j},{i-1,j+1},{i,j+1},{i+1,j+1},{i+1,j},{i+1,j-1},{i,j-1},{i-1,j-1},{i-1,j}};
                            
                            //1. Count number of white neighbour pixels
                            count_b = 0;
                            for(int k=0; k<8; k++){
                                if(image_ptr[window[k][0]][window[k][1]] == 255)
                                    count_b++; 
                            }
                            if(count_b < 2 || count_b > 6)
                                continue;
                            
                            //2. Count number of black -> white transitions (closing circle in window 2,3,...,2)
                            count_a = 0;
                            for(int k=0; k<9; k++){
                                if(image_ptr[window[k][0]][window[k][1]] == 0)
                                    BtoW_change = true;
                                else if(BtoW_change && image_ptr[window[k][0]][window[k][1]] == 255){
                                    BtoW_change = false;
                                    count_a++;
                                }
                                     
                            }
                            if(count_a != 1)
                                continue;
    
                            //3. check for at least 2,4,6 are black pixels
                            if(image_ptr[window[0][0]][window[0][1]] *image_ptr[window[2][0]][window[2][1]] * image_ptr[window[4][0]][window[4][1]] != 0)
                                continue;
                            
                            //4. check for at least 4,6,8 are black pixels
                            if(image_ptr[window[6][0]][window[6][1]] *image_ptr[window[2][0]][window[2][1]] * image_ptr[window[4][0]][window[4][1]] != 0)
                                continue;
                            
                            //set 0 in current pixel
                            image_ptr[i][j] = 0;   
                            change_flag = true;
                        }
                    }
                }
    
                //phase 2
                for(int i = 1; i < image[n]->getRows()-1; i++){
                    for(int j = 1; j < image[n]->getCols()-1; j++){
                        if(image_ptr[i][j] == 255){
                            int window[9][2] = {{i-1,j},{i-1,j+1},{i,j+1},{i+1,j+1},{i+1,j},{i+1,j-1},{i,j-1},{i-1,j-1},{i-1,j}};
                            
                            //1. Count number of white neighbour pixels
                            count_b = 0;
                            for(int k=0; k<8; k++){
                                if(image_ptr[window[k][0]][window[k][1]] == 255)
                                    count_b++; 
                            }
                            if(count_b < 2 || count_b > 6)
                                continue;
                            
                            //2. Count number of black -> white transitions (closing circle in window 2,3,...,2)
                            count_a = 0;
                            for(int k=0; k<9; k++){
                                if(image_ptr[window[k][0]][window[k][1]] == 0)
                                    BtoW_change = true;
                                else if(BtoW_change && image_ptr[window[k][0]][window[k][1]] == 255){
                                    BtoW_change = false;
                                    count_a++;
                                }
                            }
                            if(count_a != 1)
                                continue;
    
                            //3. check for at least 2,4,8 are black pixels
                            if(image_ptr[window[0][0]][window[0][1]] *image_ptr[window[2][0]][window[2][1]] * image_ptr[window[6][0]][window[6][1]] != 0)
                                continue;
                            
                            //4. check for at least 2,6,8 are black pixels
                            if(image_ptr[window[0][0]][window[0][1]] *image_ptr[window[4][0]][window[4][1]] * image_ptr[window[6][0]][window[6][1]] != 0)
                                continue;
                            
                            //set 0 in current pixel
                            image_ptr[i][j] = 0;   
                            change_flag = true;
                        }
                    }
                }
            }

            //save into folder
            cout<<save_path_skeleton + to_string(n + db_init) +"_skeleton.pgm"<<endl;
            image[n]->pgmWrite(save_path_skeleton + to_string(n+db_init)+"_skeleton.pgm","Skeletonized image with Shang-Suen algorithm");
        }
        
        
    }

    /*reset enhanced image array and confusion matrix values*/
    void clearArray(int array_type){
        if(array_type == 1 && (int)image.size() > 0 ){
            for(int i= 0; i<(int)image.size(); i++){
                delete image[i];
            }
            image.clear();
        }

        if(array_type == 2 && (int)mask.size() > 0 ){
            for(int i= 0; i<(int)mask.size(); i++){
                delete mask[i];
            }
            mask.clear();
        }

        if(array_type == 3 && (int)groundtruth.size() > 0 ){
            for(int i= 0; i<(int)groundtruth.size(); i++){
                delete groundtruth[i];
            }
            groundtruth.clear();
        }

        if(array_type == 4 && (int)segmented.size() > 0 ){
            for(int i= 0; i<(int)segmented.size(); i++){
                delete segmented[i];
            }
            segmented.clear();
        }


        confusion[0] = 0;
        confusion[1] = 0;
        confusion[2] = 0;
        confusion[3] = 0;
    }

    /*print metrics F1-score, Jaccard, Accuracy and Recall*/
    void metrics(){
        float dice = 2*confusion[0] / (2*confusion[0] + confusion[2] + confusion[3]);
        float jaccard = confusion[0] / (confusion[0] + confusion[2] + confusion[3]);
        float accuracy = (confusion[0] + confusion[1]) / (confusion[0] + confusion[1] + confusion[2] + confusion[3]);
        float recall = (confusion[0]) / (confusion[0] + confusion[3]); 

        cout<<"|---------------Results-----------------|\n";
        cout<<"F1-score         "<<dice<<endl;
        cout<<"Jaccard          "<<jaccard<<endl;
        cout<<"Accuracy         "<<accuracy<<endl;
        cout<<"Recall           "<<recall<<endl;
    }

    /*print formated data from ROC  testing*/
    void printROCData(string method, string params){
        float dice = 2*confusion[0] / (2*confusion[0] + confusion[2] + confusion[3]);

        cout << setw(20) << left << method << setw(10) << left << params << setw(10) << dice << endl;
    }
        
};



/*create structuring element
    allowed shapes: square, cross, disk, line, diamond
*/
int** createStrel(string shape, int fill_n, int rows=0, int cols=0, int radius=0, int angle=0){

    int** strel;

    //create structuring element
    if(radius > 0){
        strel = createStructuringElement(shape, fill_n, 0, 0, radius);
    }
    else if(rows > 0 && cols > 0){
        radius = rows/2;
        if(angle == 0)
            strel = createStructuringElement(shape, fill_n,rows, cols);
        else{
            strel = createStructuringElement(shape, fill_n, rows, cols, 0, angle);
            
        }
    }

    return strel;
}

/*rotate kernel*/
int** rotatekernel(int** kernel,int rows, int cols, int angle, int inplace = false){
    int** kernel_rotated;
    if(inplace)
        kernel_rotated = kernel;
    else 
        kernel_rotated = createMatrix(rows,cols,0);

    int center_i = rows/2;
    int y_center;
    int center_j = cols/2;
    int x_center;

    //rotation kernel
    double rotation[2][2] = {{cos(angle*M_PI/180),-sin(angle*M_PI/180)},{sin(angle*M_PI/180),cos(angle*M_PI/180)}};

    //new rotateed coordinates
    int i;
    int j;
    
    //iterator over patches pixels
    for ( int y = 0; y < rows; y++ ){
        for ( int x = 0; x < cols; x++ ){
            //center coordinates
            y_center = y - center_i;
            x_center = x - center_j;
            //rotation coordinates
            i = round(y_center*rotation[0][0] + x_center*rotation[0][1]);
            j = round(y_center*rotation[1][0] + x_center*rotation[1][1]);

            //center new coordinates

            if(i+center_i >= 0 && i+center_i < rows && j+center_j >= 0 && j+center_j < cols)
                kernel_rotated[y][x] = kernel[i+center_i][j+center_j];
        }
    }

    return kernel_rotated;
}

/*create gaussian matching filter from parameters
    [0]-sigma, [1]-L_len, [2]-T_len
*/
int** createGMFkernel(int* gmf_params, int extra_L = 4, int extra_T = 6){
    
    double** kernel = createDoubleMatrix(gmf_params[1],gmf_params[2],0);

    int center_j = round(gmf_params[2]/2);

    //sum of elements
    double sum = 0;
    double scale = 10;


    for(int i = 0; i < gmf_params[1]; i++){
        for(int j = 0; j < center_j+1; j++){
            //gaussian distribution function
            kernel[i][center_j + j] = - exp(-((double)(j*j)/(double)(2*(gmf_params[0])*(gmf_params[0]))));
            kernel[i][center_j - j] = kernel[i][center_j + j];
            if(j != 0)
                sum += 2*kernel[i][center_j + j];
            else
                sum += kernel[i][center_j + j];
        }
    }
    double mean = sum / (gmf_params[1]*gmf_params[2]);

    for(int i = 0; i < gmf_params[1]; i++){
        for(int j = 0; j < center_j+1; j++){
            //gaussian distribution function
            kernel[i][center_j + j] =  (kernel[i][center_j + j] - mean);
            kernel[i][center_j - j] = kernel[i][center_j + j];
        }
    }

    //fit in bigger canvas
    int rows_extend = gmf_params[1] + extra_L;
    int cols_extend = gmf_params[2] + extra_T;

    int** kernel_fit = createMatrix(rows_extend,cols_extend,0);

    for(int i = 0; i < gmf_params[1]; i++){
        for(int j = 0; j < gmf_params[2]; j++){
            kernel_fit[(extra_L/2) + i][(extra_T/2) + j] = round(scale*kernel[i][j]);
        }
    }
    

    //clear memory
    delete kernel[0];
    delete kernel;
    return kernel_fit;
}

/*create binary descriptor from initial structure and change percentage*/
int** randomBinaryDescriptor(int** initial_strel, int change_percent, int radius){

    //create binary descriptor from initial structure
    int** strel = copyImage(initial_strel,radius*2+1,radius*2+1);
    //initiallize random generator
    srand(time(0));
    //pixel coordinates to change
    int x,y;
    //number of pixel changes
    int n_changes = (int)((change_percent/100.0) * (radius*2+1)*(radius*2+1));
    for(int i = 0; i < n_changes; i++){
        x = rand() % (radius*2+1);
        y = rand() % (radius*2+1);
        strel[x][y] = rand() % 2;
    }

    return strel;
}


/*evaluate strel with whole dataset, returning its AUC and save images*/
void addStrelEvaluation(ROC* roc_curve, int** strel, string strel_name, int strel_param[], int enhancetype){
    //image operation result
    int **img_bright;
    int **img_dim;

    //image object pointer
    Image *img;

    //pgm descriptor
    string pgm_desc;
    if(enhancetype == 1)
        pgm_desc = "image - topkhat, ";
    else
        pgm_desc = "image + tophat - blackhat, ";


    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        img = new Image();
        //build path to image 
        img->pgmRead(db_path + to_string(i) +"_training.pgm");

        //Apply morphological operations
        if(enhancetype == 1){
            //enhance original image by decreasing light (image - blackhat)
            img_bright = img->morphOp("tophat",strel,strel_param[3]);
            img->diffImage(img_bright);

        }else{
            //enhance original image by increasing contrast (image + tophat - blackhat)
            img_bright = img->morphOp("tophat",strel,strel_param[3]);
            img_dim = img->morphOp("blackhat",strel,strel_param[3]); 
            img->addImage(img_bright);
            img->diffImage(img_dim);
            delete img_dim[0];
            delete img_dim;
        }
        
        // Save the result
        //img->pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm", pgm_desc + strel_name);
        // Invert image
        img->invertImage();
        // Store into ROC array object
        roc_curve->insertEnhanceImage(img);
        //clear memory
        delete img_bright[0];
        delete img_bright;
        
    }

}

/*enhance whole dataset with morphological kernels, saving images*/
void enhanceDataset(int** strel, string strel_name, int strel_param[], int enhancetype, int ref_path){
    //image operation result
    int **img_bright;
    int **img_dim;

    //image object pointer
    Image *img;

    //pgm descriptor
    string pgm_desc;
    if(enhancetype == 1)
        pgm_desc = "image - topkhat, ";
    else
        pgm_desc = "image + tophat - blackhat, ";


    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        img = new Image();
        //path to image to enhance
        if(ref_path == 1) 
            img->pgmRead(db_path + to_string(i) +"_training.pgm");
        else
            img->pgmRead(save_path_enhance + to_string(i) + "_enhance.pgm");

        //Apply morphological operations
        if(enhancetype == 1){
            //enhance original image by decreasing light (image - blackhat)
            img_bright = img->morphOp("tophat",strel,strel_param[3]);
            img->diffImage(img_bright);

        }else{
            //enhance original image by increasing contrast (image + tophat - blackhat)
            img_bright = img->morphOp("tophat",strel,strel_param[3]);
            img_dim = img->morphOp("blackhat",strel,strel_param[3]); 
            img->addImage(img_bright);
            img->diffImage(img_dim);
            delete img_dim[0];
            delete img_dim;
        }
        
        // Save the result
        img->pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm", pgm_desc + strel_name);

        //clear memory
        delete img_bright[0];
        delete img_bright;
        
    }

}



/*Local Search algorithm to improve strel response*/
void localSearch(ROC* roc_best, int** strel,int strel_params[], int change_percent, int radius, int iterations){

    //strel to impove
    int** strel_aux;

    //temporal roc pointer
    ROC *roc_temp;

    //initial best ROC object
    addStrelEvaluation(roc_best,strel,"binary descriptor",strel_params,2);
    //confusion matrix
    roc_best->calculateConfusionMatrix(true);
    //sensitivity and specificity
    roc_best->calculateSensSpec();
    //AUC
    roc_best->calculateAUC();
    cout<<"Initial local AUC: "<<roc_best->getArea()<<endl;

    //auxiliary ROC object
    ROC *roc_aux = new ROC("binary_desc",strel_params);
    //add training mask
    roc_aux->buildMaskArray(mask_path,db_size,db_init);
    //add training groundthruth
    roc_aux->buildGroundthruthArray(gt_path,db_size,db_init);

    for(int i = 0; i< iterations; i++){
        //perturbate strel at "change percent" of pixels
        strel_aux = randomBinaryDescriptor(strel,change_percent,radius);

        //evaluate strel
        addStrelEvaluation(roc_aux,strel_aux,"binary descriptor",strel_params,2);
        roc_aux->calculateConfusionMatrix(true);
        roc_aux->calculateSensSpec();
        roc_aux->calculateAUC();
        cout<<"local search AUC: "<<roc_aux->getArea()<<endl;

        if(roc_aux->getArea() > roc_best->getArea()){
            
            roc_temp = roc_best;
            roc_best = roc_aux;
            roc_aux = roc_temp;

            delete strel[0];
            delete strel;
            strel = strel_aux;
        }else{
  
            delete strel_aux[0];
            delete strel_aux;
        }
        roc_aux->clearEnhanceArray();
    }


    //clear memory
    delete roc_aux;
}


/*Apply iterated local search to improve initial strel response*/
ROC *iteratedLocalSearch(int** init_strel, int radius, int iterations){

    //strel params
    int** strel_aux;
    int strel_param[5] = {1,0,0,radius,0};

    //temporal roc pointer
    ROC *roc_temp;

    //initial best ROC object
    ROC *roc_best = new ROC("binary_desc",strel_param);
    addStrelEvaluation(roc_best,init_strel,"binary descriptor",strel_param,2);
    //add training mask
    roc_best->buildMaskArray(mask_path,db_size,db_init);
    //add training groundthruth
    roc_best->buildGroundthruthArray(gt_path,db_size,db_init);
    //confusion matrix
    roc_best->calculateConfusionMatrix(true);
    //sensitivity and specificity
    roc_best->calculateSensSpec();
    //AUC
    roc_best->calculateAUC();
    cout<<"Initial AUC: "<<roc_best->getArea()<<endl;

    //auxiliary ROC object
    ROC *roc_aux = new ROC("binary_desc",strel_param);
    //add training mask
    roc_aux->buildMaskArray(mask_path,db_size,db_init);
    //add training groundthruth
    roc_aux->buildGroundthruthArray(gt_path,db_size,db_init);

    for(int i = 0; i< iterations; i++){
        //perturbate strel at 75% of pixels
        strel_aux = randomBinaryDescriptor(init_strel,75,radius);
        //local search to improve strel respoonse
        localSearch(roc_aux,strel_aux,strel_param,15,radius,iterations);
        cout<<"best local search AUC: "<<roc_best->getArea()<<endl;
        //evaluate strel
        if(roc_aux->getArea() > roc_best->getArea()){
            roc_temp = roc_best;
            roc_best = roc_aux;
            roc_aux = roc_temp;

            delete init_strel[0];
            delete init_strel;
            init_strel = strel_aux;
        }else{
            delete strel_aux[0];
            delete strel_aux;
        }
        roc_aux->clearEnhanceArray();
        
    }

    //clear memory
    delete roc_aux;
    return roc_best;
}


/*Evaluate a single strel params with the whole dataset returning the ROC as an object*/
ROC *strelType(string strel_name, int strel_param[], int enhancetype){

    ROC* roc_curve;
    int** strel;
    //create strel evaluation from parameters
    if (strel_name != "binary descriptor"){
        strel = createStrel(strel_name,strel_param[0],strel_param[1],strel_param[2],strel_param[3],strel_param[4]);
        roc_curve = new ROC(strel_name,strel_param);
        addStrelEvaluation(roc_curve,strel,strel_name,strel_param,enhancetype);

    }else{
        strel = createStrel("diamond",strel_param[0],strel_param[1],strel_param[2],strel_param[3],strel_param[4]);
        roc_curve = iteratedLocalSearch(strel,strel_param[3],10);
        //print best strel
        Image *img = new Image();
        img->setImage(strel, strel_param[3]*2+1,strel_param[3]*2+1);
        img->pgmWrite("best_strel.pgm","best strel from ILS algorithm");
        delete img;
        
    }
    
    delete strel[0];
    delete strel;
    return roc_curve;
}

/*Test different parameters for image enhancement, returning the best set of strel and size*/
void testEnhanceParams(string* strel_names, int* strel_params, int* strel_radii, int n_tests, int enhancetype){

    //array of ROC evaluation per strel combination
    ROC* roc_array[n_tests];

    //evalute every strel name and size combination
    cout << setw(20) << left << "|Strel" << setw(10) << left << "|Radio" << setw(10) << "|AUC" << endl;

    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 4; j++){
            strel_params[1] = strel_radii[j]*2 + 1;
            strel_params[2] = strel_radii[j]*2 + 1;
            strel_params[3] = strel_radii[j];
            //enhance images
            roc_array[j + 6*i] = strelType(strel_names[i],strel_params,enhancetype);
            //add training mask
            roc_array[j + 6*i]->buildMaskArray(mask_path,db_size,db_init);
            //add training groundthruth
            roc_array[j + 6*i]->buildGroundthruthArray(gt_path,db_size,db_init);
            //confusion matrix
            roc_array[j + 6*i]->calculateConfusionMatrix(true);
            //sensitivity and specificity
            roc_array[j + 6*i]->calculateSensSpec();
            //AUC
            roc_array[j + 6*i]->calculateAUC();
            roc_array[j + 6*i]->printROCData();
        }
    }
}

/*Enhance images applying traditional symetric structuring element*/
void enhanceSymetricStrel(string strel_name, int strel_params[], int enhancetype,int ref_path){
    
    int** strel = createStrel(strel_name,strel_params[0],strel_params[1],strel_params[2],strel_params[3],strel_params[4]);
    enhanceDataset(strel,strel_name,strel_params,enhancetype,ref_path);

    delete strel[0];
    delete strel;
}

/*Enhance images applying the best binary descriptor from ILS algorithm*/
void enhanceBinaryDescriptor(int* strel_params){
    strel_params[3] = 8;
    int** strel = createStrel("diamond",strel_params[0],strel_params[1],strel_params[2],strel_params[3],strel_params[4]);
    ROC* roc_curve = iteratedLocalSearch(strel,strel_params[3],5);
    roc_curve->printROCData();

    //print best strel
    Image *img = new Image();
    img->setImage(strel, strel_params[3]*2+1,strel_params[3]*2+1);
    img->pgmWrite("best_strel.pgm","best strel from ILS algorithm");

    delete img;
    delete roc_curve;
    delete strel[0];
    delete strel;
}

/*enhance images with matching gaussian filter*/
void gaussianMatchingFilter(int* gmf_params, int ref_path){

    //create gaussian matching filter
    int extraL = gmf_params[1]/2;
    int extraT = gmf_params[2]/2;
    int** gmf_kernel = createGMFkernel(gmf_params,extraL,extraT);
    int** gmf_rotated[12];
    //create image matrix
    int** img_matrix = createMatrix(584,565,0);
    int** img_aux;

    Image *img;

    //compute rotational kernels
    Image *kernel = new Image();
    for(int j = 0; j < 12; j++){
        //rotate kernel
        gmf_rotated[j] = rotatekernel(gmf_kernel,gmf_params[1]+extraL,gmf_params[2]+extraT,15*(j),false);
        //print kernel
        kernel->setImage(gmf_rotated[j],gmf_params[1]+extraL,gmf_params[2]+extraT,true,true);
        kernel->normalize();
        kernel->pgmWrite("kernel" + to_string(15*(j)) + ".pgm","kernel rotated");
        
    }

    //apply filter to datset
    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        img = new Image();

        //path to image to enhance
        if(ref_path == 1) 
            img->pgmRead(db_path + to_string(i) +"_training.pgm");
        else
            img->pgmRead(save_path_enhance + to_string(i) + "_enhance.pgm");

        for(int j = 0; j < 12; j++){

            //apply filter
            img_aux = img->convolution(gmf_rotated[j],gmf_params[1]+extraL,gmf_params[2]+extraT,false);

            //selecting max response angle for every pixel
            for(int k=0; k < img->getRows(); k++){
                for(int l = 0; l < img->getCols(); l++){
                    if(img_matrix[k][l] < img_aux[k][l]){
                        img_matrix[k][l] = img_aux[k][l];
                    }
                }
            }
            delete img_aux[0];
            delete img_aux;
        }
        
        // Set resulting image
        img_matrix = img->normalize(img_matrix);
        img->pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm","image enhanced with gaussian matching filter",img_matrix);
        delete img;
        
    }

    delete kernel;
}

/*soft the hiighest valued gradient edge of the set of images*/
void ROI(int threshold){
    int** img_matrix;
    int** mask;
    int** strel = createStrel("square",1,0,0,2,0);

    Image img[db_size];
    Image mask_img;

    //dataset 
    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        img[i-db_init].pgmRead(db_path + to_string(i) +"_training.pgm");

        //apply sharr edge detection
        img_matrix = img[i-db_init].scharr_gradient(false);
        img_matrix = img[i-db_init].normalize(img_matrix);
        
        //selecting gradient values above threshold
        mask = createMatrix(img[i-db_init].getRows(),img[i-db_init].getCols(),0);
        mask_img.setImage(mask,img[i-db_init].getRows(),img[i-db_init].getCols(),false);
        for(int k=0; k < img[i-db_init].getRows(); k++){
            for(int l = 0; l < img[i-db_init].getCols(); l++){
                
                if(img_matrix[k][l] > threshold){
                        mask[k][l] = threshold;
                    }
            }
        }
        //dilate over edge
        for(int k = 0; k < 50; k++){
            mask_img.morphOp("dilation",strel,2,true);
            //mean
            mask_img.meanWindow();
        }
        
        //cut inner growth dilation
        mask_img.fillCountour(img_matrix);

        //add to original image
        img[i-db_init].addCountour(mask_img.getImage());

        img[i-db_init].pgmWrite("dilated_result.pgm","edge dilation result",mask_img.getImage());

        //clear memory
        delete img_matrix[0];
        delete img_matrix;

        // Set resulting image
        img[i-db_init].pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm","image enhanced with ROI to soft edge");
    }
}

/*Smooth set of images with gaussian filtering*/
void smoothImages(int ref_path){
    Image img;
    //mask;

    //apply filter to dataset
    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        //path to image to enhance
        if(ref_path == 1) 
            img.pgmRead(db_path + to_string(i) +"_training.pgm");
        else
            img.pgmRead(save_path_enhance + to_string(i) + "_enhance.pgm");

        //read mask
        //mask.pgmRead(mask_path + to_string(i) +"_training_mask.pgm");
        img.gauss_filter(true);
        
        // Set resulting image
        img.pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm","image with inverted values");
    }
}

/*invert set of images */
void invertImages(int ref_path){
    Image img,mask;

    //apply filter to dataset
    for(int i = db_init; i < db_init + db_size; i++){
        // Read image
        //path to image to enhance
        if(ref_path == 1) 
            img.pgmRead(db_path + to_string(i) +"_training.pgm");
        else
            img.pgmRead(save_path_enhance + to_string(i) + "_enhance.pgm");

        //read mask
        mask.pgmRead(mask_path + to_string(i) +"_training_mask.pgm");
        img.invertImage(mask.getImage());
        
        // Set resulting image
        img.pgmWrite(save_path_enhance + to_string(i) + "_enhance.pgm","image with inverted values");
    }
}

/*Test different parameters for image segmentation, returning the best set of params*/
void testSegmentParams(){

    int strel_params[5] = {1,0,0,8,0};
    //testing segmentation parameters
    int connected_thresh[] = {50, 60, 70, 80};

    //testing max gradient threshold (Yanowitz method)
    int max_thresh[] = {20, 40, 60, 80};

    //enhance pipeline
    //smooth original images
    smoothImages(1);
    //increase contrast whitehat - blackhat
    enhanceSymetricStrel("diamond",strel_params,2,2);
    enhanceSymetricStrel("diamond",strel_params,2,2);
    enhanceSymetricStrel("diamond",strel_params,2,2);

    //evalute every strel name and size combination
    cout << setw(20) << left << "|Method" << setw(10) << left << "|max-thresh; connect-thresh" << setw(10) << "|F-1 Score" << endl;

    Segment drive_training;
    drive_training.buildImageArray(save_path_enhance,db_size,db_init);
    drive_training.buildMaskArray(mask_path,db_size,db_init);
    drive_training.buildGroundthruthArray(gt_path,db_size,db_init);

    //connect threshold 
    for(int i = 0; i < 4; i++){

        //max gradient threshold
        for(int j = 0; j < 4; j++){
            //threshold surface method
            drive_training.yanowitz_method(db_size,save_path_segment,max_thresh[j], connected_thresh[i]);
            drive_training.calculateConfusionMatrix();
            drive_training.printROCData("Surface",to_string(max_thresh[j]) + "," + to_string(connected_thresh[i]));
            drive_training.clearArray(4);
        }


        //iterative threshold method
        drive_training.iterative_method(db_size,save_path_segment,connected_thresh[i]);
        drive_training.calculateConfusionMatrix();
        drive_training.printROCData("Iterative",";" + to_string(connected_thresh[i]));
        drive_training.clearArray(4);
        
    }
}

/*compute metrics with enhanced images without segmentation*/
void computeMetrics(){
    
    Segment drive_training;
    drive_training.buildImageArray(save_path_enhance,db_size,db_init,2);
    drive_training.buildMaskArray(mask_path,db_size,db_init);
    drive_training.buildGroundthruthArray(gt_path,db_size,db_init);
    drive_training.calculateConfusionMatrix();
    drive_training.metrics();
    string temp;
    cin>>temp;
}

/*Segment using surface of images*/
void segmentSurfaceYanowitz(int maxima_t, int threshold){
    Segment drive_training;
    drive_training.buildImageArray(save_path_enhance,db_size,db_init);
    drive_training.buildMaskArray(mask_path,db_size,db_init);
    drive_training.buildGroundthruthArray(gt_path,db_size,db_init);
    drive_training.yanowitz_method(db_size,save_path_segment,maxima_t,threshold);
    drive_training.calculateConfusionMatrix();
    drive_training.metrics();
    cout<<">>Segmentation process finished"<<endl;
    string temp;
    cin>>temp;
}

/*Segment using iterative thresholding computation*/
void segmentSurfaceIterative(int threshold){
    Segment drive_training;
    drive_training.buildImageArray(save_path_enhance,db_size,db_init);
    drive_training.buildMaskArray(mask_path,db_size,db_init);
    drive_training.buildGroundthruthArray(gt_path,db_size,db_init);
    drive_training.iterative_method(db_size,save_path_segment,threshold);
    drive_training.calculateConfusionMatrix();
    drive_training.metrics();
    cout<<">>Segmentation process finished"<<endl;
    string temp;
    cin>>temp;
}

/*skeletonize from segmented images*/
void skeletonization(){
    Segment image_segmented;
    image_segmented.buildImageArray(save_path_segment,db_size,db_init,1,2);
    image_segmented.buildMaskArray(mask_path,db_size,db_init);
    image_segmented.zhangSuenSkeletonization();
    cout<<">>Skeletonization process finished"<<endl;
    string temp;
    cin>>temp;
}

/*calculate vessel width */
void vesselWidth(){

    string save_path_width =  "src/db_coronary/width/";
    string save_path_skeleton = "src/db_coronary/skeletonized/";
    int vessel_t = 10;  // vessel max threshold
    int scale = 50;     // scaling factor to visualize radii image
    Image edge, skeleton, mask;
    int** radii_image;

    //stats variables
    double max = 0, min = 1000, avg = 0, n_pixels = 0;

    cout << setw(13) << left << "|Image ID" << setw(14) << left << "|Min (pixel)" << setw(14) << "|Max (pixel)"  << setw(14) << left << "|Avg (pixel)" <<endl;
    
    for(int i = db_init; i < db_init + db_size; i++){
        //compute Canny edge detection
        edge.pgmRead(save_path_segment + to_string(i) + "_segmented.pgm");
        edge.cannyEdge();
        edge.pgmWrite(save_path_width + to_string(i)+"_edge.pgm","Canny edge detection");
        
        //load ROI from mask
        mask.pgmRead(mask_path + to_string(i) + "_training_mask.pgm");
        
        //compare Canny edge with skeletonized vessel
        skeleton.pgmRead(save_path_skeleton + to_string(i) + "_skeleton.pgm");
        radii_image = skeleton.radialEdgeSearch(edge.getImage(),mask.getImage(),vessel_t,scale);

        //print radii map
        skeleton.pgmWrite(save_path_width + to_string(i)+"_radii.pgm","Radii map",radii_image);

        //compute vessel stats
        for(int k = 0; k < skeleton.getRows(); k++){
            for(int l = 0; l < skeleton.getCols(); l++){
                if(radii_image[k][l] == 0)
                    continue;

                if(max < radii_image[k][l]/scale)
                    max = radii_image[k][l] / scale;
                if(min > radii_image[k][l]/scale){
                    if(radii_image[k][l] < 0)
                        min = 0.5;  //possible fraction represented as unique negative value
                    else
                        min = radii_image[k][l] / scale;
                }
                    
                n_pixels++;
                if(radii_image[k][l] < 0)
                    avg += 0.5;
                else
                    avg += (radii_image[k][l] / scale);
            }
        }

        cout << setw(13) << left << "  "+to_string(i) << setw(14) << left << "  " + to_string(min) << setw(14) << "  " + to_string(max)  << setw(14) << left << avg/n_pixels <<endl;
        max = 0;
        min = 1000;
        avg = 0; 
        n_pixels = 0;

        delete radii_image[0];
        delete radii_image;
    }

    cout<<">>Vessel width process finished"<<endl;
    string temp;
    cin>>temp;
}

/*solve equation system with Cramer's determinants*/
void solveEquationCramer(int A[3][3], int b[3], double coef[3]){
    //determinants
    double det_s = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    double det_A = b[0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1]*(b[1]*A[2][2] - A[1][2]*b[2]) + A[0][2]*(b[1]*A[2][1] - A[1][1]*b[2]);
    double det_B = A[0][0]*(b[1]*A[2][2] - A[1][2]*b[2]) - b[0]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*b[2] - b[1]*A[2][0]);
    double det_C = A[0][0]*(A[1][1]*b[2] - b[1]*A[2][1]) - A[0][1]*(A[1][0]*b[2] - b[1]*A[2][0]) + b[0]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    coef[0] = det_A/det_s;
    coef[1] = det_B/det_s;
    coef[2] = det_C/det_s;
}

/*Select 2 random symetrci points on each side of the optical disk*/
void randomVesselPoint(int** img_matrix, int rows, int cols, int A[3][3],int b[3]){
    //initiallize random generator
    //srand(time(0));
    //select random row and column until it find any vessel pixel
    int range_x, range_y, column = 0, row = 0;

    if(b[0] < cols/2)
        range_x = cols - (b[0]) + 1;    //range from optical disk to the image width
    else
        range_x = (b[0]) - 0 + 1;
                
    range_y = (A[0][1]) - 0 + 1;        //range from optical disk to the image heigth
    
    while(img_matrix[row][column] == 0){
        if(b[0] < cols/2)
            column = rand() % range_x + (b[0]);
        else
            column = rand() % range_x;
        row = rand() % range_y + 0; 
    }

    A[1][0] = row*row;
    A[1][1] = row;
    A[1][2] = 1;
    b[1] = column;

    //(second coordinate) select random row and column until it find any vessel pixel
    column = 0, row = 0;
    range_y =  rows - (A[0][1]) + 1;    //range from optical disk to the image heigth
    
    while(img_matrix[row][column] == 0){
        if(b[0] < cols/2)
            column = rand() % range_x + (b[0]);
        else
            column = rand() % range_x;
        row = rand() % range_y + (A[0][1]); 
    }
    
    A[2][0] = row*row;
    A[2][1] = row;
    A[2][2] = 1;
    b[2] = column;
}

/*select best fitting parabola from 3 points using RANSAC algorithm*/
void RANSACparabola(int** mask,int** segmented, int rows, int cols, int yx_opticdisk[2],  double* best_parabola, int max_iterations, int eps){
    double parabola_params[3];
    int iteration = 0, inliers = 0, best_inliers = 0;
    int min_y, max_y, x;
    int A[3][3],b[3];   //parabola parameters equation system

    //build matrix with first equation as optic disk center
    A[0][0] = yx_opticdisk[0]*yx_opticdisk[0];
    A[0][1] = yx_opticdisk[0];
    A[0][2] = 1;
    b[0] = yx_opticdisk[1];

    while(iteration < max_iterations){
        inliers = 0;
        //select 2 random additional points
        
        randomVesselPoint(segmented,rows,cols, A, b);
        cout<<iteration<<": "<<b[1]<<","<<A[1][1]<<" - "<<b[0]<<","<<A[0][1]<<" - "<<b[2]<<","<<A[2][1] <<endl;
        
        //Compute parabola coefficients Ay^2 + By + C = x
        solveEquationCramer(A,b,parabola_params);
        
        //determine range of ROI
        for(int i = 1; i< rows; i++){
            if(mask[i][cols/2] != 0 && min_y == 0)
                min_y = i;
            if(mask[i][cols/2] == 0 && mask[i-1][cols/2] != 0 && max_y == 0)
                max_y = i;
        }
        
        //count inliers from ROI around epsilon distance of proposed parabola
        for(int y = min_y; y <= max_y; y++){
            //compute Ay^2 + By + C = x
            x = round(parabola_params[0]*y*y + parabola_params[1]*y + parabola_params[2]);
            if(x > cols || x < 0)
                continue;
            
            //check epsilon window  for inliers
            for(int k = y-eps; k <= y + eps; k++){
                for(int l = x-eps; l <= x + eps; l++){
                    if(segmented[k][l] != 0)
                        inliers++;
                }
            }
        }
        
        if(inliers > best_inliers){
            best_parabola[0] = parabola_params[0];
            best_parabola[1] = parabola_params[1];
            best_parabola[2] = parabola_params[2];
            best_inliers = inliers;
        }
        
        iteration++;
    }

}

/*Compute best parabolic model of Mayor Temporal Arcade (MTA) */
void MTAModeling(){
    
    string save_path_MTA =  "src/db_coronary/MTA_model/";
    int yx_opticdisk[2] = {0};
    double parabola_params[3];
    int x;
    
    Image image, mask, segmented;
    
    for(int i = db_init; i < db_init + db_size; i++){
        //load ROI from mask
        mask.pgmRead(mask_path + to_string(i) + "_training_mask.pgm");
        
        //load segmented image
        segmented.pgmRead(save_path_segment + to_string(i) + "_segmented.pgm");
        //Find optic disk center
        image.pgmRead(db_path + to_string(i) + "_training.pgm");
        image.maxCoordinates(yx_opticdisk);
        
        //best parabola selection
        RANSACparabola(mask.getImage(),segmented.getImage(),mask.getRows(),mask.getCols(), yx_opticdisk, parabola_params,10, 5);

        //plotting parabola
        int** parabola_img = image.getCopyImage();//createMatrix(image.getRows(),image.getCols(),0);
        for(int k = 0; k < image.getRows(); k++){
            for(int l = 0; l < image.getCols(); l++){
                //(compute Ay^2 + By + C = x)
                x = round(parabola_params[0]*k*k + parabola_params[1]*k + parabola_params[2]);
                if(x > image.getCols() || x < 0)
                    continue;

                parabola_img[k][x] = 255;
            }
        }
        
        cout<<save_path_MTA + to_string(i)+"_parabola.pgm"<<endl;
        image.pgmWrite(save_path_MTA + to_string(i)+"_parabola.pgm","Best adjusted parabola ",parabola_img);
        delete parabola_img[0];
        delete parabola_img;
    }

    cout<<">>MTA modeling process finished"<<endl;
    string temp;
    cin>>temp;
    
}


/*Define dataset subpaths basedon root path*/
void setDatasetPaths(string db_path){
    //set dataset paths

    //path to store enhanced images
    save_path_enhance = db_path + "training_enhance/";
    //path to store segmented images
    save_path_segment = db_path + "training_segmentation/";
    //path for masks
    mask_path = db_path + "training/mask/";
    //path for groundtruth
    gt_path = db_path + "training/groundtruth/";

}

void morphological_interface(int enhancetype, int ref_path){
    int option = 100;
    string strel_names[]= {"diamond","disk"};
    int strel_radii[] = {2,4,6,8};
    int strel_params[5] = {1,0,0,0,0};
    int n_tests = 2*4;

    while (option != 0){
        /*morphological submenu */
        cout<<"|----------------------------------|\n";
        cout<<"| 1. ROC evaluation                |\n";
        cout<<"| 2. Diamond                       |\n";
        cout<<"| 3. Round                         |\n";
        cout<<"| 4. Binary descriptor             |\n";
        cout<<"| 0. Return                        |\n";
        cout<<"|----------------------------------|\n";
        cout<<"Select a structuring element: ";
        cin>>option;
        
        switch (option){
        case 1:
            //------------------------------------------------test for differents enhance params
            //structuring element parameters: weight, rows, cols, radius, angle=0
            testEnhanceParams(strel_names,strel_params,strel_radii,n_tests,enhancetype-1);
            break;
        case 2:
            strel_params[3] = 8;
            cout<<"Procesando imagenes..."<<endl;
            enhanceSymetricStrel(strel_names[0],strel_params,enhancetype-1,ref_path);
            break;
        case 3:
            strel_params[3] = 8;
            enhanceSymetricStrel(strel_names[1],strel_params,enhancetype-1,ref_path);
            break;
        case 4:
            enhanceBinaryDescriptor(strel_params);
            break;
        case 0:
            break;
        default:
            cout<<"Invalid option"<<endl;
            break;
        }
        cout<< u8"\033[2J\033[1;1H"; //clear console
    }
}

void enhance_interface(){
    int option = 100;
    int threshold = 200;
    int ref_path = 1;
    int gmf_params[3] = {2,9,13};

    while (option != 0){
        /*Enhance submenu */
        cout<<"|----------------------------------|\n";
        cout<<"| 1. ROI                           |\n";
        cout<<"| 2. I - Tophat                    |\n";
        cout<<"| 3. I + Tophat - Blackhat         |\n";
        cout<<"| 4. Gaussian smooth filter        |\n";
        cout<<"| 5. Gaussian Matching Filter      |\n";
        cout<<"| 6. Invert values                 |\n";
        
        cout<<"\n\n| 8. Set original image path       |\n";
        cout<<"| 9. Set last enhance as new image |\n";
        cout<<"| 0. Return                        |\n";

        cout<<"|----------------------------------|\n";
        cout<<"Select an option: ";
        cin>>option;
        
        switch (option)
        {
        case 0:
            break;
        case 1:
            cout<<"Select a threshold: ";
            cin>>threshold;
            ROI(threshold);
            
            break;
        case 2:
            morphological_interface(option,ref_path);
            break;
        case 3:
            morphological_interface(option,ref_path);
            break;
        case 4:
            smoothImages(ref_path);
            break;
        case 5:
            gaussianMatchingFilter(gmf_params,ref_path);
            break;
        case 6:
            invertImages(ref_path);
            break;
        case 7:

            break;
        case 8:
            ref_path = 1;
            break;
        case 9:
            ref_path = 2;
            break;
        default:
            cout<<"Invalid option"<<endl;
            break;
        }
        cout<< u8"\033[2J\033[1;1H"; //clear console
    }
}

void segment_interface(){
    int option = 100;
    int threshold;
    int maxima_t;

    while (option != 0){
        /*Segment submenu */
        cout<<"|----------------------------------|\n";
        cout<<"| 1. Thresholding surface          |\n";
        cout<<"| 2. Convex Hull                   |\n";
        cout<<"| 3. Iterative thresholding        |\n";
        cout<<"| 4. Evaluate methods              |\n";
        cout<<"| 0. Return                        |\n";
        cout<<"|----------------------------------|\n";
        cout<<"Select an option: ";
        cin>>option;
        
        switch (option)
        {
        case 1:
            cout<<"Connecting element threshold: ";
            cin>>threshold;
            cout<<"threshold for local maxima: ";
            cin>>maxima_t;
            segmentSurfaceYanowitz(maxima_t,threshold);

            break;
        case 2:
            
            break;
        case 3:
            cout<<"Connecting element threshold: ";
            cin>>threshold;
            segmentSurfaceIterative(threshold);
            break;
        case 4:
            testSegmentParams();
            break;
        case 0:
            break;
        default:
            cout<<"Invalid option"<<endl;
            break;
        }
        cout<< u8"\033[2J\033[1;1H"; //clear console
    }
}

void interface(){
    int option = 100;
    while (option != 0)
    {
        cout<<"|----------------------------------|\n";
        cout<<"| 1. Enhance images                |\n";
        cout<<"| 2. Segment images                |\n";
        cout<<"| 3. Skeletonization               |\n";
        cout<<"| 4. Vessel radius                 |\n";
        cout<<"| 5. MTA model                     |\n";
        cout<<"| 0. Exit                          |\n";
        cout<<"|----------------------------------|\n";
        cout<<"Select an option: ";
        
        cin>>option;
        switch (option)
        {
        case 1:
            enhance_interface();
            break;
        case 2:
            segment_interface();
            break;
        case 3:
            skeletonization();
            break;
        case 4:
            vesselWidth();
            break;
        case 5:
            MTAModeling();
        case 0:
            break;
        default:
            cout<<"Invalid option"<<endl;
            break;
        }
        cout<< u8"\033[2J\033[1;1H";
    }
    
     
}



int main(int argc, char **argv){
    if(argc != 4){
        cout << "Error, params: 1. db_path, 2.db_size, 3.db_init" << endl;
        return 1;
    }

    //dataset path
    db_path = argv[1] + string("training/");
    db_size = atoi(argv[2]);
    db_init = atoi(argv[3]);
    //set dataset paths
    setDatasetPaths(argv[1]);

    //init interface
    interface();
    
    return 0;
}