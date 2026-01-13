#include "init.hpp"

namespace sfem
{
    //=============================================================================
    Application &initialize(int argc, char *argv[],
                            bool write_log_file,
                            const std::string &name)
    {
        return Application::instance(argc, argv, name, write_log_file);
    }
}