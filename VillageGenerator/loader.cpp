#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"


// FUNCTION DEFINITIONS

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

std::vector<std::array<int, 2>> getLinePoints(uint8_t *img, int dims[], int sample) {
    // Produces a sample of points along an input line of black points in {(x1,y1),(x2,y2),...} format

    int height,width,hits;
    width = dims[0];
    height = dims[1];
    hits = 0;

    std::vector<std::array<int,2>> line_points;
    std::array<int, 2> position = {0,0};

    for (int j = 0; j< height; j++) {
        for (int i = 0; i<width; i++) {
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

void arrayToImage(uint8_t *img, std::vector<std::array<int,2>> points, int dims[], const char *filename, int chnnls, int qual) {
    // Follows std_image_write formatting to produce a white .jpg file with each sample point colored in black

    int height, width;
    width = dims[0];
    height = dims[1];

    for (int j = 0; j< height; j++) {
        for (int i = 0; i<width; i++) {
            img[j * width + i] = 255;
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
}

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

//---------------------------

// Graveyard:
    // Prints particular image pixel's value:
    /*int pixel_position[2] = {223,631};
    getPixelPrint(image, pixel_position, dimensions)*/



int main()
{
    //std::string filename = "village_line.jpg";
    int width, height, channels, ok, dims;
    //ok = stbi_info("village_line.jpg", &w, &h, &comp);
    uint8_t *image = stbi_load("village_line.jpg", &width, &height, &channels, 1);
    if(image == NULL) {
        std::cout << "Error in loading the image\n" << std::endl;
        exit(1);
    }

    int dimensions[2];
    dimensions[0]=width;
    dimensions[1]=height;
    

    // Writing an image of test pixels
    int sample_rate = 9; // Plots 1 out of every <sample_rate> points 
    std::vector<std::array<int, 2>> data = getLinePoints(image, dimensions, sample_rate);


    int quality = 100; // out of 100
    arrayToImage(image, data, dimensions, "test.jpg", 1, quality);


    return 0;
}



