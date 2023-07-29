#include "Segment.hpp"

// Constructors
Segment::Segment() : from_(0), to_(0), score_(0.0) {}
Segment::Segment(int from, int to, double score) : from_(from), to_(to), score_(score) {}

// Getters and setters
int Segment::getFrom() const { return from_; }
void Segment::setFrom(int from) { from_ = from; }
int Segment::getTo() const { return to_; }
void Segment::setTo(int to) { to_ = to; }

int Segment::getFromOrg() const { return from_org_; }
void Segment::setFromOrg(int from) { from_org_ = from; }
int Segment::getToOrg() const { return to_org_; }
void Segment::setToOrg(int to) { to_org_ = to; }

double Segment::getScore() const { return score_; }
void Segment::setScore(double score) { score_ = score; }

std::string Segment::getChain() const { return chain_; }
void Segment::setChain(std::string chain) { chain_ = chain; }


// Overloaded operators
bool Segment::operator==(const Segment& other) const {
    return (from_ == other.from_) && (to_ == other.to_) && (score_ == other.score_);
}
bool Segment::operator!=(const Segment& other) const { return !(*this == other); }
bool Segment::operator<(const Segment& other) const {
    if (from_ != other.from_) {
        return from_ < other.from_;
    }
    return to_ < other.to_;
}
bool Segment::operator>(const Segment& other) const { return other < *this; }
bool Segment::operator<=(const Segment& other) const { return !(other < *this); }
bool Segment::operator>=(const Segment& other) const { return !(*this < other); }

std::ostream& operator<<(std::ostream& os, const Segment& segment) {
  //os << "Segment [from=" << segment.from_ << ", to=" << segment.to_ << ", score=" << segment.score_ << "]";
  //  os <<  segment.from_ + 1 << "-" << segment.to_ + 1;
  os << "(resi " << segment.from_org_ << "-" << segment.to_org_ << " and chain " << segment.chain_ + ")";
  return os;
}
