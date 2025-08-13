#include "test.h"

// This is test code based off of the OOP LinkedIn Learning course. Look, I haven't used C++ in like 8 years.


Review::Review(unsigned int r, const std::string& t, const std::string& txt)
{
    rating = r;
    title = t;
    text = txt;
}

Review::~Review()
{
    std::cout << "Review object is deleted." << std::endl;
}

void Review::displayDetails() const
{
    std::cout << "Rating: " << rating << "/5\nTitle: " << title << "\nText: " << text << std::endl;
}



