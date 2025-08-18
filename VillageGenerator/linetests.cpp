#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>

int main()
{
    //std::vector<std::array<int, 2>, unsigned int> points;

    //std::vector<std::array<int, 2>, unsigned int> line_points;
    std::vector<std::vector<int>, unsigned int> line_points;
    std::vector<int> position = {0, 0};

    position = {1, 2};
    line_points.push_back(position);

    position = {2, 8};
    line_points.push_back(position);


    for (int n : line_points) {
        std::cout << n << std::endl;
    }

    return 0;
}
