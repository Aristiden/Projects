#ifndef TEST_H
#define TEST_H

//Header guard to avoid double defining class stuff.

#include <iostream>

// This is test code based off of the OOP LinkedIn Learning course. Look, I haven't used C++ in like 8 years.

class Review
{
public:
    Review(unsigned int r, const std::string& t, const std::string& txt);
    ~Review();
    void displayDetails() const;
    

private:
    unsigned int rating;
    std::string title;
    std::string text;
};

#endif

