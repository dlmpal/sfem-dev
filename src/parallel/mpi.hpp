#pragma once

#include <vector>
#include <span>
#include <source_location>

#define SFEM_MPI_CHECK_ERROR(error_code)  \
    {                                     \
        if (error_code != MPI_SUCCESS)    \
            sfem::mpi::error(error_code); \
    }

namespace sfem::mpi
{
    // Convenience functions for initializing and finalizing MPI
    void initialize(int argc, char *argv[]);
    void finalize();

    /// @brief Get the root process rank
    int root();

    ///@brief Get the MPI process rank
    int rank();

    /// @brief Get the number of MPI processes
    int n_procs();

    /// @brief Abort the application
    void abort(int error = -1);

    /// @brief Handle a nonzero code returned by an MPI function
    void error(int code, std::source_location location = std::source_location::current());

    /// @brief Send data from all processes to all processes
    /// @param data Data
    /// @param dest Destination processes
    /// @return Data receiced from all processes,
    /// along with the received size per process and the displacements
    template <typename T>
    std::tuple<std::vector<T>, std::vector<int>, std::vector<int>>
    send_to_dest(const std::span<const T> data, const std::span<const int> dest);

    /// @brief Distribute (broadcast) data from the root process to all processes
    /// @param data Data
    /// @param dest Destination proceses
    /// @return Data belonging to this process
    template <typename T>
    std::vector<T> distribute(const std::span<const T> data, const std::span<const int> dest);
}