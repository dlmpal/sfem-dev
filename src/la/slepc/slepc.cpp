#ifdef SFEM_HAS_SLEPC

#include "slepc.hpp"

namespace sfem::la::slepc
{
    //=============================================================================
    void initialize(int argc, char *argv[])
    {
        SlepcInitialize(&argc, &argv, nullptr, nullptr);
    }
    //=============================================================================
    void finalize()
    {
        SlepcFinalize();
    }
}

#endif // SFEM_HAS_SLEPC