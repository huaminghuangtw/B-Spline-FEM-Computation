#include "utilities.hpp"
#include <stdexcept>

namespace cie
{
namespace splinekernel
{

void runtime_check(bool result, const char message[])
{
    if( !result )
    {
        throw std::runtime_error{ message };
    }
}

} // namespace splinekernel
} // namespace cie