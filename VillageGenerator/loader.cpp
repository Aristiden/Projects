#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

// DOCUMENTATION of stb_image.h
//
// Limitations:
//    - no 12-bit-per-channel JPEG
//    - no JPEGs with arithmetic coding
//    - GIF always returns *comp=4
//
// Basic usage (see HDR discussion below for HDR usage):
//    int x,y,n;
//    unsigned char *data = stbi_load(filename, &x, &y, &n, 0);
//    // ... process data if not NULL ...
//    // ... x = width, y = height, n = # 8-bit components per pixel ...
//    // ... replace '0' with '1'..'4' to force that many components per pixel
//    // ... but 'n' will always be the number that it would have been if you said 0
//    stbi_image_free(data);
//
// Standard parameters:
//    int *x                 -- outputs image width in pixels
//    int *y                 -- outputs image height in pixels
//    int *channels_in_file  -- outputs # of image components in image file
//    int desired_channels   -- if non-zero, # of image components requested in result

// To query the width, height and component count of an image without having to
// decode the full file, you can use the stbi_info family of functions:
//
//   int x,y,n,ok;
//   ok = stbi_info(filename, &x, &y, &n);
//   // returns ok=1 and sets x, y, n if image is a supported format,
//   // 0 otherwise.

int getPixel(uint8_t *img,int pixel_pos[], int dims[]) {
    unsigned int pixel, width,x,y;
    x = pixel_pos[0];
    y = pixel_pos[1];
    width = dims[0];
    pixel = img[y * width + x];
    return pixel;
}

void getPixelPrint(uint8_t *img,int pixel_pos[], int dims[]) {
    std::cout << "Value of pixel at (" << pixel_pos[0]<<","<<pixel_pos[1]  << ") = " << getPixel(img, pixel_pos, dims) << std::endl;
}

std::vector<std::array<int, 2>, unsigned int> getLinePoints(uint8_t *img, int dims[], int sample) {
    int height,width,hits;
    width = dims[0];
    height = dims[1];
    hits = 0;

    std::vector<std::array<int,2>, unsigned int> line_points;
    std::array<int, 2> position = {0,0};

    for (int j = 0; j< height; j++) {
        for (int i = 0; i<width; i++) {
            if (img[j * width + i] < 128) {
                
                hits++;
                if (hits % sample == 0) {
                    position = {i,j};
                    line_points.push_back(position);
                } 

            /*    img[j * width + i] = 1;
            } else {
                img[j * width + i] = 0;
            }*/
            }
        }
    }
    return line_points;
}

uint8_t *arrayToImage(std::vector<std::array<int,2>, unsigned int> points, int dims[]) {
    uint8_t *img;
    int height,width;
    width = dims[0];
    height = dims[1];

    int length = points.size();
    int i,j;
    for (int s = 0; s< length; s++) {
        i = points[s][0];
        j = points[s][1];
        img[j * width + i] = 1;
    }

    return img;
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


    // Prints each image pixel value in a list:
    /*unsigned int pixel;
    for (int j = 0; j< height; j++){
        for (int i = 0; i<width; i++)
            pixel = image[j * width + i];
            std::cout << pixel << std::endl;
    }*/
    


    // Prints particular image pixel's value:
    /*int pixel_position[2] = {223,631};
    getPixelPrint(image, pixel_position, dimensions)*/
    


    // Writing an image of test pixels
    int sample_rate = 4; // Plots 1 out of every <sample_rate> points 
    std::vector<std::array<int, 2>, unsigned int> data = getLinePoints(image, dimensions, sample_rate);
    uint8_t *output = arrayToImage(data, dimensions);
    int quality = 80; // out of 100
    stbi_write_jpg("test.jpg", width, height, channels, output, quality);

    return 0;
}



