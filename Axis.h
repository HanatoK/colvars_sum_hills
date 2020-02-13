#ifndef AXIS_H
#define AXIS_H
#include <cstddef>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>

using std::string;
using std::vector;

class Axis
{
public:
    Axis();
    Axis(double lowerBound, double upperBound, size_t bin, bool periodic = false);
    bool isInitialized() const;
    double width() const;
    size_t bin() const;
    bool isInBoundary(double value) const;
    bool isPeriodic() const;
    bool wrap(double& x) const;
    bool index(double x, size_t& idx) const;
    size_t index(double x) const;
    double distance(double start, double end) const;
    string infoHeader() const;
    vector<double> middlePoint() const;
    friend bool operator==(const Axis& lhs, const Axis& rhs);
private:
    double          mLowerBound;
    double          mUpperBound;
    size_t          mBin;
    double          mWidth;
    bool            mPeriodic;
};

Axis::Axis():
mLowerBound(0.0), mUpperBound(0.0), mBin(0), mWidth(0.0), mPeriodic(false) {}

Axis::Axis(double lowerBound, double upperBound, size_t bin, bool periodic):
mLowerBound(lowerBound), mUpperBound(upperBound), mBin(bin), mWidth((upperBound - lowerBound) / static_cast<double>(bin)), mPeriodic(periodic) {}

bool Axis::isInitialized() const {
    return static_cast<bool>(mBin);
}

double Axis::width() const {
    return mWidth;
}

size_t Axis::bin() const {
    return mBin;
}

bool Axis::isInBoundary(double value) const {
    if (value < mLowerBound || value > mUpperBound) {
        return false;
    } else {
        return true;
    }
}

bool Axis::isPeriodic() const {
    return mPeriodic;
}

bool Axis::wrap(double& x) const {
    if (!isPeriodic()) {
        if (!isInBoundary(x)) return false;
        else return true;
    } else {
        const double period = mUpperBound - mLowerBound;
        if (x > mUpperBound) {
            double dist_to_bound = x - mUpperBound;
            x = x - (std::floor(dist_to_bound / period) + 1) * period;
            return true;
        } else if (x < mLowerBound) {
            double dist_to_bound = mLowerBound - x;
            x = x + (std::floor(dist_to_bound / period) + 1) * period;
            return true;
        } else {
            return true;
        }
    }
}

bool Axis::index(double x, size_t& idx) const {
    if (!wrap(x)) return false;
    idx = std::floor((x - mLowerBound) / mWidth);
    if (idx == mBin) {
        --idx;
    }
    return true;
}

size_t Axis::index(double x) const {
    size_t idx = std::floor((x - mLowerBound) / mWidth);
    idx = (idx == mBin) ? (mBin - 1) : idx;
    return idx;
}

double Axis::distance(double start, double end) const {
    if (!isPeriodic()) {
        return (end - start);
    } else {
        wrap(start);
        wrap(end);
        const double period = mUpperBound - mLowerBound;
        const double dist = end - start;
        if (std::abs(dist) > (period * 0.5)) {
            if (start > end) {
                return (dist + period);
            } else if (start < end) {
                return (dist - period);
            } else {
                return dist;
            }
        } else {
            return dist;
        }
    }
}

string Axis::infoHeader() const {
    std::stringstream ss;
    ss << "# " << std::fixed << std::setprecision(9)
       << mLowerBound << ' ' << mWidth << ' ' << mBin << ' ' << short(mPeriodic);
    return ss.str();
}

vector<double> Axis::middlePoint() const {
    double tmp = mLowerBound - 0.5 * mWidth;
    vector<double> result(mBin, 0.0);
    for(auto &i : result) {
        tmp += mWidth;
        i = tmp;
    }
    return result;
}

bool operator==(const Axis& lhs, const Axis& rhs) {
    if (lhs.mLowerBound != rhs.mLowerBound) {
        return false;
    }
    if (lhs.mUpperBound != rhs.mUpperBound) {
        return false;
    }
    if (lhs.mBin != rhs.mBin) {
        return false;
    }
    return true;
}

#endif
