#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <cmath>
#define _USE_MATH_DEFINES
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"



// DEBUG FUNCTIONS

void doesItMakeIt(int iter = 0) {
    std::cout << "yep" << std::endl;
    if (iter > 0) {
        std::cout << "yep, " << iter << std::endl;
    }
}

void intPrint(int print) {
    std::cout << print << std::endl;
}

void floatPrint(float print) {
    std::cout << print << std::endl;
}

void stringPrint(std::string print) {
    std::cout << print << std::endl;
}


// FUNCTION DEFINITIONS

float gaussianRoller(int n) {
    // Produces a pseudorandom number between 0 and 1 from a pseudo Gaussian centered at 0.5
    //   More n produces more Gaussian!

    float output = 0;
    srand(time(0));

    for (int roll = 0; roll < n; roll++) {
        output += rand() % 6;
    }
    output = output/(n*5);
   
    return output;
}

int getPixel(uint8_t *img,int pixel_pos[], int dims[]) {
    // This function takes in an image loaded as a 1D array and a desired pixel position, then returns the grayscale value of that pixel
    //    255 = white, 1 (or 0?) = black

    unsigned int pixel, width,x,y;
    x = pixel_pos[0];
    y = pixel_pos[1];
    width = dims[0];
    pixel = img[y * width + x];
    return pixel;
}

std::vector<std::array<int, 2>> getLineFull(uint8_t *img, int dims[]) {
    // Produces all points along an input line of black points in {(x1,y1),(x2,y2),...} format

    int height,width,hits;
    width = dims[0];
    height = dims[1];
    hits = 0;

    std::vector<std::array<int,2>> line_points;
    std::array<int, 2> position = {0,0};

    for (int j = 0; j < height; j++) {
        for (int i = width-1; i >= 0; i--) {
            if (img[j * width + i] < 128) {
                position = {i,j};
                line_points.push_back(position);
                 
            }
        }
    }
    return line_points;
}

std::vector<std::array<int, 2>> getLinePoints(uint8_t *img, int dims[], int sample) {
    // Produces a sample of points along an input line of black points in {(x1,y1),(x2,y2),...} format

    int height,width,hits;
    width = dims[0];
    height = dims[1];
    hits = 0;

    std::vector<std::array<int,2>> line_points;
    std::array<int, 2> position = {0,0};

    for (int j = 0; j < height; j++) {
        for (int i = width-1; i >= 0; i--) {
            if (img[j * width + i] < 128) {
                
                hits++;
                if (hits % sample == 0) {
                    position = {i,j};
                    line_points.push_back(position);
                } 
            }
        }
    }
    return line_points;
}

std::vector<std::array<int, 4>> getPointNormals(std::vector<std::array<int, 2>> line_points) {
    
    int x1,x2, y1,y2, dx,dy;

    std::vector<std::array<int, 4>> point_normals;
    for (int s = 0; s < line_points.size()-1; s++) {
        int e;
        x1 = line_points[s][0];
        y1 = line_points[s][1];

        x2 = line_points[s+1][0];
        y2 = line_points[s+1][1];

        dx = x2 - x1; 
        dy = y2 - y1;

        // The distance from each normal point to the line is equal to the distance between point 1 and point 2
        point_normals.push_back({x1, y1, -dy + x1, dx + y1});

    }
    return point_normals;
}



std::vector<std::array<int, 4>> flipNormalsVertical (std::vector<std::array<int, 4>> nrmals, int dims[]) {
    // Takes the normals = {x1,y1, n1,n2} format and flips the vertical coordinates around the middle of the total image

    int height = dims[1]-1;

    for (int s = 0; s < nrmals.size(); s++) {
        nrmals[s][1] = height - nrmals[s][1];
        nrmals[s][3] = height - nrmals[s][3];
    }

    return nrmals;
}

std::array<std::array<int, 2>, 4> flipBuilding (std::array<std::array<int, 2>, 4> building_corners, int dims[]) {
    // Takes in an array of 4 (x,y) points and flips all y components around the middle of the total image

    int height = dims[1]-1;

    for (int s = 0; s < 4; s++) {
        building_corners[s][1] = height - building_corners[s][1];
    }
    return building_corners;
}


std::array<std::array<int, 2>, 4> drawBuilding(int corner_index, std::array<int,2> ftprint, std::vector<std::array<int, 4>> normals_list, int dims[]) {
    // Pixels get counted from top to bottom, so all of my math is upside down. I'm just going to flip it here and then flip it back at the end:
    std::vector<std::array<int, 4>> nrmals;
    nrmals = flipNormalsVertical(normals_list, dims);
    
    std::array<std::array<int, 2>, 4> building_corners;
    // Points are listed counterclockwise
    int length = ftprint[0];
    int depth = ftprint[1];
    //std::cout << "nrmals[corner_index][2] = " << nrmals[corner_index][2] << std::endl;

    // Bottom left corner, cornerstone:
    int x1 = nrmals[corner_index][2];
    int y1 = nrmals[corner_index][3];
    building_corners[0] = {x1, y1};
    
    // Can't actually use x1, y1 for calculations unfortunately because it may need to be changed if the building doesn't fit the road.

    std::array<int, 2> test_normal;
    int i = 0;
    float dist = 0;
    int dx,dy;
    float test_func;

    while (dist < length) {  // Finds another normal to draw the front wall along
        i++;
        test_normal = {nrmals[corner_index-i][2], nrmals[corner_index-i][3]};

        dist = sqrt(pow(building_corners[0][0]-test_normal[0],2) + pow(building_corners[0][1]-test_normal[1],2)); // Prolly should make this into a function
        dx = test_normal[0] - building_corners[0][0];
        dy = test_normal[1] - building_corners[0][1];
        for(int k = 0; k < i; k++) {
            test_func = (dy/dx)*(nrmals[corner_index-k][0] - test_normal[0]) + test_normal[1];
            if (nrmals[corner_index - k][1] > test_func) {
                stringPrint("fails");
                corner_index--;
                building_corners[0] = {nrmals[corner_index][2], nrmals[corner_index][3]};
                i = 0, dist = 0;
            }
        }
    }

    float theta = atan(dy/dx);
    std::cout << "theta = " << theta << std::endl;
    int x2 = building_corners[0][0] + ceil(length*cos(theta));
    int y2 = building_corners[0][1] + ceil(length*sin(theta));
    building_corners[1] = {(int)x2, (int)y2};

    int x3 = round(x2 - depth*cos(M_PI_2 - theta));
    int y3 = round(y2 + depth*sin(M_PI_2 - theta));
    building_corners[2] = {(int)x3, (int)y3};

    int x4 = round(building_corners[0][0] - depth*cos(M_PI_2 - theta));
    int y4 = round(building_corners[0][1] + depth*sin(M_PI_2 - theta));
    building_corners[3] = {(int)x4, (int)y4};
    
    // Need to: 
    //   - Still draw building along line to neighboring normal point if the next point doesn't intersect (function test(s))

    building_corners = flipBuilding(building_corners, dims);
    return building_corners;
}

std::vector<std::array<int, 2>> drawManyBuildings(int n_buildings, std::array<int,2> ftprint, std::vector<std::array<int, 4>> normals_list, int dims[]) {

    std::vector<std::array<int, 2>> building_data;
    std::array<std::array<int, 2>, 4> building;
    int random_index;
    srand(time(0));
    for (int n = 0; n < n_buildings; n++) {
        //random_index = round(normals_list.size()*gaussianRoller(10));
        random_index = rand() % normals_list.size();
        building = drawBuilding(random_index, ftprint, normals_list, dims);
        for (int s = 0; s < 4; s++) {
            building_data.push_back(building[s]);
        }
    }

    return building_data;
}

// PRINT FORMATTING/VISUALIZATION FUNCTIONS



void getPixelPrint(uint8_t *img,int pixel_pos[], int dims[]) {
    // Function that prints the grayscale value of a desired pixel to the terminal (see getPixel)

    std::cout << "Value of pixel at (" << pixel_pos[0]<<","<<pixel_pos[1]  << ") = " << getPixel(img, pixel_pos, dims) << std::endl;
}

void printAllPixels(uint8_t *img, int dims[]) {
    // Function that prints all grayscale pixel values in a loaded image

    int height, width;
    width = dims[0];
    height = dims[1];
    unsigned int pixel;

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            pixel = img[j * width + i];
            std::cout << pixel << std::endl;
        }
    }
}

void printPoints(std::vector<std::array<int,2>> points) {
    // Function that prints each point pair in a list of (x,y) points

    for (int s = 0; s < points.size(); s++) {
        std::cout << "("
        << points[s][0] << ","
        << points[s][1] << ")" << std::endl;
        
    }
}

uint8_t *arrayToImage(uint8_t *img, std::vector<std::array<int,2>> points, int dims[], const char *filename, int chnnls, int qual, bool blank) {
    // Follows std_image_write formatting to produce a white .jpg file with each sample point colored in black

    int height, width;
    width = dims[0];
    height = dims[1];

    if (blank == true) {
        for (int j = 0; j< height; j++) {
            for (int i = 0; i<width; i++) {
                img[j * width + i] = 255;
            }
        }
    }

    int length = points.size();
    int i,j;
    for (int s = 0; s< length; s++) {
        i = points[s][0];
        j = points[s][1];
        img[j * width + i] = 1;
    }

    stbi_write_jpg(filename, width, height, chnnls, img, qual);
    return img;
}



//---------------------------



int main()
{
    //std::string filename = "village_line.jpg";
    int width, height, channels, ok, dims;
    //ok = stbi_info("village_line.jpg", &w, &h, &comp);
    uint8_t *image = stbi_load("village_line.jpg", &width, &height, &channels, 1);
    if(image == NULL) {
        std::cout << "Error loading the image\n" << std::endl;
        exit(1);
    }

    int dimensions[2];
    dimensions[0]=width;
    dimensions[1]=height;
    

    // Writing an image of test pixels
    int sample_rate = 9; // Plots 1 out of every <sample_rate> points 
    std::vector<std::array<int, 2>> data = getLinePoints(image, dimensions, sample_rate);


    int quality = 100; // out of 100
    bool blank = true;
    uint8_t *sampled_image = arrayToImage(image, data, dimensions, "test.jpg", 1, quality, blank);
    blank = false;
    std::vector<std::array<int, 4>> normals = getPointNormals(data);
    std::vector<std::array<int, 2>> normals_data;
    for (int s = 0; s < normals.size(); s++) {
        normals_data.push_back({normals[s][2], normals[s][3]});
    }
    // printPoints(normals_data);
    uint8_t *normals_image = arrayToImage(sampled_image, normals_data, dimensions, "testNormals.jpg", 1, quality, blank);


    // Probably should read in a text file with several preset building footprints and possibly an assoicated priority int
    std::array<int,2> footprint = {20,10};
    
    /*int random_index = round(normals.size()*gaussianRoller(10));
    std::array<std::array<int, 2>, 4> building = drawBuilding(random_index, footprint, normals, dimensions);

    std::vector<std::array<int, 2>> building_data;
    for (int s = 0; s < 4; s++) {
        building_data.push_back(building[s]);
    }
    printPoints(building_data);*/

    std::vector<std::array<int, 2>> building_data;
    int number_of_buildings = 8;
    building_data = drawManyBuildings(number_of_buildings, footprint, normals, dimensions);
    blank = true;
    arrayToImage(normals_image, building_data, dimensions, "testBuilding.jpg", 1, quality, blank);
    

    image = stbi_load("village_line.jpg", &width, &height, &channels, 1);
    std::vector<std::array<int, 2>> original_line = getLineFull(image, dimensions);
    std::vector<std::array<int, 2>> total_data;
    for (int s = 0; s < original_line.size(); s++) {
        total_data.push_back(original_line[s]);
    }
    for (int s = 0; s < building_data.size(); s++) {
        total_data.push_back(building_data[s]);
    }

    blank = false;
    arrayToImage(image, total_data, dimensions, "testManyBuildings.jpg", 1, quality, blank);
    //blank = false;
    //arrayToImage(image, building_data, dimensions, "testManyBuildings.jpg", 1, quality, blank);


    return 0;
}

// Graveyard:
    // Prints particular image pixel's value:
    /*int pixel_position[2] = {223,631};
    getPixelPrint(image, pixel_position, dimensions)*/


    /*uint8_t *processImage(uint8_t *img, int dims[]) {
        int height,width;
        width = dims[0];
        height = dims[1];

        for (int j = 0; j< height; j++) {
            for (int i = 0; i<width; i++) {
                if (img[j * width + i] < 128) {
                    img[j * width + i] = 1;
                } else {
                    img[j * width + i] = 0;
                }
            }
        }
        return img;
    }*/


/*int getHits(uint8_t *img, int dims[]) {
    // Unnecessary function, but useful if you have to initialize a vector with an explicit length

    int height,width,hits;
    width = dims[0];
    height = dims[1];
    hits = 0;

    for (int j = 0; j< height; j++) {
        for (int i = 0; i<width; i++) {
            if (img[j * width + i] < 128) {
                hits++;
            }
        }
    }
    return hits;
}*/
