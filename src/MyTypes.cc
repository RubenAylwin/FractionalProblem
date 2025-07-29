#include <MyTypes.h>

/**
 * @brief: Imaginary unit.
 */
BEM::Complex BEM::I(0.,1.);

/**
 * @brief: Implementation for printing intervals.
 */
std::ostream& operator<<(std::ostream &os, const BEM::Interval1D& interval)
{
    os << "(" << interval.first << ", " << interval.second << ")";
    return os;
}
