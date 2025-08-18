#include "test.h"
#include "test.cpp"

int main()
{
    Review review = Review( 5, "Excellent", "Exceeded my expectations");
    review.displayDetails();


    review.setRating(4);
    review.setTitle("Blood Meridian");
    review.setText("The width x length of America will destroy us all, and our private gods with it.");
    review.displayDetails();


    return 0;
}
