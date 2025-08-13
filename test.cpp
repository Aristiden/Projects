#include "test.h"

// This is test code based off of the OOP LinkedIn Learning course. Look, I haven't used C++ in like 8 years.

Review::Review(unsigned int r, const std::string &t, const std::string &txt)
{
    setRating(r);
    setTitle(t);
    setText(txt);
}

Review::~Review()
{
    std::cout << "Review object is deleted." << std::endl;
}

void Review::displayDetails() const
{
    std::cout << "Rating: " << rating << "/5\nTitle: " << title << "\nText: " << text << std::endl;
}

void Review::setRating(unsigned int r)
{
    if (r < min_rating || r > max_rating)
    {
        throw std::invalid_argument("Rating must be between 1 and 5");
    }
    rating = r;
}

void Review::setTitle(const std::string &t)
{
    title = validateAndTrim(t, max_title_length, "Title");
}

void Review::setText(const std::string &txt)
{
    text = validateAndTrim(txt, max_text_length, "Review text");
}

std::string Review::validateAndTrim(const std::string &str,
                            unsigned int maxLength,
                            const std::string &fieldName) const
{
    if (str.empty())
    {
        throw std::invalid_argument(fieldName + " cannot be empty");
    }

    return str.length() > maxLength ? str.substr(0, maxLength) : str;
}
