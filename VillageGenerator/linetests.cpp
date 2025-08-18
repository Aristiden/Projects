#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>

int main()
{
    //std::vector<std::array<int, 2>, unsigned int> points;

    std::vector<std::array<int, 2>> line_points;
    std::array<int, 2> position{0, 0};

    //std::vector<std::vector<int>, int> line_points;
    //std::vector<int> position{0, 0};

    std::vector<int> test{0, 1};
    test.push_back(3);
    for (int n : test) {
        //std::cout << n << std::endl;
    }

    //std::vector<std::vector<int>> line_points;

    position = {1, 2};
    line_points.push_back(position);

    position = {2, 8};
    line_points.push_back(position);
    for (int s = 0; s < line_points.size(); s++) {
        for (int i = 0; i < 2; i++) {
            std::cout << line_points[s][i] << std::endl;
        }
    }
    
    return 0;
}
