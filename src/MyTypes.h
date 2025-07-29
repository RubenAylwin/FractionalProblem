#ifndef MY_TYPES
#define MY_TYPES
#include <complex>
#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <chrono>
#include <initializer_list>
#include <list>
#include <iostream>
#include <Point.h>

//////////////////////////////////////////
// File with custom types and typedefs. //
//////////////////////////////////////////

class EmpiricalInterpolation;

namespace BEM {;
    class Interval1D;
    using Complex = std::complex<double>;
    using Support1D = std::pair<double, double>;
    using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixP = Eigen::Matrix<std::complex<long double>, Eigen::Dynamic, Eigen::Dynamic>;
    using ColVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;
    using OneDimFunction = std::function<BEM::Complex (double t)>;
    using TwoDimFunction = std::function<BEM::Complex (double t, double s)>;
    using Transformation = std::function<Point2D (Point2D P)>;
    using UVector = std::vector<unsigned>;
    using DVector = std::vector<double>;
    using CVector = std::vector<BEM::Complex>;
    using TimePoint = std::chrono::time_point<std::chrono::system_clock>;
    extern Complex I;
    using EmpiricalInterpolationPtr = std::shared_ptr<EmpiricalInterpolation>;

    /**
     * @brief: Auxiliary class for interval.
     */
    class Interval1D
    {
    public:
        Interval1D(double first_, double second_) : first{first_}, second{second_} {};
        Interval1D operator+(double shift) const {return Interval1D(first+shift, second+shift);}
        double first = 0.0;
        double second = 0.0;
    };

    /**
     * @brief: Auxiliary class for list of intervals.
     */
    class Support1DL : public std::list<Interval1D>
    {
    public:
        Support1DL(void) = default;
        Support1DL(std::initializer_list<Interval1D> list) {
            this->insert(this->begin(), list.begin(), list.end());
        };
        Support1DL(Interval1D interval) : std::list<Interval1D>() {
            this->push_back(interval);
        };
        Support1DL operator+(double shift) const {
            Support1DL other{};
            for (const auto &interval : *this) {
                other.push_back(interval + shift);
            }
            return other;
        };
    };
}

/**
 * @brief: To print an interval.
 */
std::ostream& operator<<(std::ostream &os, const BEM::Interval1D& interval);

/**
 * @brief: To substract two vectors. TODO: Use concepts.
 */
template <typename T>
std::vector<T> operator-(const std::vector<T> &v1, const std::vector<T> &v2)
{
    assert(v1.size() == v2.size());
    auto result = v1;
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] -= v2[i];
    }
    return result;
}

/**
 * @brief: To multiply a vector by a scalar. TODO: Use concepts.
 */
template <typename T>
std::vector<T> operator*(const std::vector<T> &v1, const BEM::Complex &c)
{
    auto result = v1;
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i]*c;
    }
    return result;
}

#endif
