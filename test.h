#ifndef TEST_H
#define TEST_H

// Header guard to avoid double defining class stuff.

#include <iostream>

// This is test code based off of the OOP LinkedIn Learning course. Look, I haven't used C++ in like 8 years.

class Review
{
public:
    Review(unsigned int r, const std::string &t, const std::string &txt);
    ~Review();
    void displayDetails() const;

    unsigned int getRating() const { return rating; }
    std::string getTitle() const { return title; }
    std::string getText() const { return text; }

    void setRating(unsigned int r);
    void setTitle(const std::string &t);
    void setText(const std::string &text);

private:
    unsigned int rating;
    std::string title;
    std::string text;

    static const unsigned int min_rating = 1;
    static const unsigned int max_rating = 5;
    static const unsigned int max_title_length = 128;
    static const unsigned int max_text_length = 2048;
    std::string validateAndTrim(const std::string &str,
                                unsigned int maxLength,
                                const std::string &fieldName) const;
};

#endif
