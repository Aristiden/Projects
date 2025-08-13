#include "test.h"
#include "test.cpp"

int main()
{
    Review review = Review( 5, "Excellent", "Exceeded my expectations");
    review.displayDetails();

    return 0;
}