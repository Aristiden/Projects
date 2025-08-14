#include "test2.h"

// Inheritance code from test.h class

SongReview::SongReview(unsigned int r, const std::string &t, const std::string &txt, const std::string &time, const std::string &gnr) : Review(r, t, txt)
{
    setSongLength(time);
    setGenre(gnr);
}

SongReview::~SongReview()
{
    std::cout << "Song review object is deleted." << std::endl;
}

void SongReview::displayDetails() const
{
    Review::displayDetails();
    std::cout << "Song Length: " << songLength << "\n"
    << "Genre: " << genre << std::endl;

}



void SongReview::setSongLength(const std::string &time)
{
    songLength = Review::validateAndTrim(time, 6, "Song length");
}

void SongReview::setGenre(const std::string &gnr)
{
    genre = Review::validateAndTrim(gnr, 128, "Song genre");
}


