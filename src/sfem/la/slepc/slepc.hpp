#pragma once

#ifdef SFEM_HAS_SLEPC

#include <slepc.h>

namespace sfem::la::slepc
{
    // Convenience functions for initializing and finalizing SLEPc
    void initialize(int argc, char *argv[]);
    void finalize();
}

#endif // SFEM_HAS_SLEPC