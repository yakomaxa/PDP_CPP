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

    int getFromOrg() const;
    void setFromOrg(int from);
    int getToOrg() const;
    void setToOrg(int to);
    std::string getChain() const;
    void setChain(std::string chain);

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
  int from_org_;
  int to_org_;
  std::string chain_;

};

#endif  // SEGMENT_HPP

