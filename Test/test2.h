#ifndef SONG_REVIEW_H
#define SONG_REVIEW_H

#include <iostream>
#include "test.h"

// Multiple inheritance works by adding comma, then the next class you want to the line below
class SongReview : public Review
{
public:
    SongReview(unsigned int r,
               const std::string &t,
               const std::string &txt,
               const std::string &time,
               const std::string &gnr);

    ~SongReview();
    void displayDetails() const;

    std::string getSongLength() const { return songLength; }
    std::string getGenre() const { return genre; }

    void setSongLength(const std::string &time);
    void setGenre(const std::string &gnr);

private:
    // Song specific parameters
    std::string songLength;
    std::string genre;


    // Composition:
    // Could also inherit features of a desired class by calling it here. "This class contains features of another class" type beat.    == Sociable social;
    //   If so, you have to define functions for this class that use functions from the called class                                    == void like() { social.like(); }
};

#endif
