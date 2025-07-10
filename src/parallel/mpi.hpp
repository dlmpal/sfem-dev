#pragma once

#include <vector>
#include <span>
#include <source_location>

/// @brief MPI-related functionality
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

    enum class ReduceOperation
    {
        min = 0,
        max = 1,
        sum = 2,
        prod = 3
    };

    /// @brief Perform a reduce operation across all processes
    /// @param value This process' value
    /// @param op Operation to be performed
    /// @return Reduced value
    template <typename T>
    T reduce(T value, ReduceOperation op);

    /// @brief Send data from all processes to all processes
    /// @param data Data
    /// @param dest Destination processes
    /// @return Data receiced from all processes,
    /// along with the received size per process and the displacements
    template <typename T>
    std::tuple<std::vector<T>, std::vector<int>, std::vector<int>>
    send_to_dest(const std::span<const T> data, const std::span<const int> dest, int block_size = 1);

    /// @brief Distribute (broadcast) data from the root process to all processes
    /// @param data Data
    /// @param dest Destination proceses
    /// @return Data belonging to this process
    template <typename T>
    std::vector<T> distribute(const std::span<const T> data, const std::span<const int> dest);
}