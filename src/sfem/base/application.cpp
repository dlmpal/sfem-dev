#include "application.hpp"
#include <sfem/parallel/mpi.hpp>
#include <sfem/la/petsc/petsc.hpp>
#include <sfem/la/slepc/slepc.hpp>
#include <iostream>
#include <format>

namespace sfem
{
    //=============================================================================
    Application::Application(int argc, char *argv[],
                             const std::string &name,
                             const std::filesystem::path &log_filename)
        : name_(name),
          log_level_(LogLevel::info),
          log_file_(log_filename)
    {
        if (log_file_.is_open() == false and log_filename.empty() == false)
        {
            log_message(std::format("Could not create log file at {}, falling back to standard streams\n", log_filename.string()),
                        LogLevel::warning);
        }

        // MPI
        mpi::initialize(argc, argv);

// PETSc
#ifdef SFEM_HAS_PETSC
        la::petsc::initialize(argc, argv);
#endif

// SLEPc
#ifdef SFEM_HAS_SLEPC
        la::slepc::initialize(argc, argv);
#endif
    }
    //=============================================================================
    Application::~Application()
    {
// SLEPc
#ifdef SFEM_HAS_PETSC
        PetscBool is_slepc_init;
        SlepcInitialized(&is_slepc_init);
        if (is_slepc_init)
            la::slepc::finalize();
#endif

// PETSc
#ifdef SFEM_HAS_PETSC
        PetscBool is_petsc_init;
        PetscInitialized(&is_petsc_init);
        if (is_petsc_init)
            la::petsc::finalize();
#endif

        // MPI
        mpi::finalize();
    }
    //=============================================================================
    void Application::log_message(const std::string &msg, LogLevel level) const
    {
        if (level < log_level_)
        {
            return;
        }

        if (log_file_.is_open())
        {
            log_file_ << msg;
        }
        else
        {
            if (level == LogLevel::error)
            {
                std::cerr << msg;
                mpi::abort();
            }
            else
            {
                std::cout << msg;
            }
        }
    }
    //=============================================================================
    std::string Application::name() const
    {
        return name_;
    }
    //=============================================================================
    void Application::set_log_level(LogLevel level)
    {
        log_level_ = level;
    }
    //=============================================================================
    Application &Application::instance(int argc, char *argv[],
                                       const std::string &name,
                                       const std::filesystem::path &log_filename)
    {
        static Application app(argc, argv, name, log_filename);
        return app;
    }
}
