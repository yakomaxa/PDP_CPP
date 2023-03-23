#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <iostream>

class Segment {
public:
    // Constructors
    Segment();
    Segment(int from, int to, double score);

    // Getters and setters
    int getFrom() const;
    void setFrom(int from);
    int getTo() const;
    void setTo(int to);
    double getScore() const;
    void setScore(double score);

    // Overloaded operators
    bool operator==(const Segment& other) const;
    bool operator!=(const Segment& other) const;
    bool operator<(const Segment& other) const;
    bool operator>(const Segment& other) const;
    bool operator<=(const Segment& other) const;
    bool operator>=(const Segment& other) const;

    // Print the segment
    friend std::ostream& operator<<(std::ostream& os, const Segment& segment);

private:
    int from_;
    int to_;
    double score_;
};

#endif  // SEGMENT_HPP

