#include "init.hpp"

namespace sfem
{
    //=============================================================================
    Application &initialize(int argc, char *argv[],
                            const std::string &name,
                            const std::filesystem::path &log_filename)
    {
        return Application::instance(argc, argv, name, log_filename);
    }
}