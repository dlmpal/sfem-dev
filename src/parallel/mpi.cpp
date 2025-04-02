#include "mpi.hpp"
#include "../base/error.hpp"
#include "../base/config.hpp"
#include <mpi.h>
#include <vector>
#include <numeric>
#include <format>

namespace sfem::mpi
{
    //=============================================================================
    void initialize(int argc, char *argv[])
    {
        MPI_Init(&argc, &argv);
    }
    //=============================================================================
    void finalize()
    {
        MPI_Finalize();
    }
    //=============================================================================
    int root()
    {
        return 0;
    }
    //=============================================================================
    int rank()
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }
    //=============================================================================
    int n_procs()
    {
        int num;
        MPI_Comm_size(MPI_COMM_WORLD, &num);
        return num;
    }
    //=============================================================================
    void abort(int error)
    {
        MPI_Abort(MPI_COMM_WORLD, error);
    }
    //=============================================================================
    void error(int error_code, std::source_location location)
    {
        char *error_msg;
        int len;
        MPI_Error_string(error_code, error_msg, &len);
        log_msg(std::format("Call to MPI function returned with error:\n\t{}\n", error_msg),
                LogLevel::error, location);
    }
    //=============================================================================
    template <typename T>
    MPI_Datatype to_mpi_datatype()
    {
        MPI_Datatype mpi_type;
        if constexpr (std::is_same_v<int, T>)
        {
            mpi_type = MPI_INT;
        }
        else if constexpr (std::is_same_v<float, T>)
        {
            mpi_type = MPI_FLOAT;
        }
        else if constexpr (std::is_same_v<double, T>)
        {
            mpi_type = MPI_DOUBLE;
        }
        return mpi_type;
    }
    // Explicit instantiations
    // int
    template MPI_Datatype to_mpi_datatype<int>();
    // real_t
    template MPI_Datatype to_mpi_datatype<real_t>();
    //=============================================================================
    template <typename T>
    std::tuple<std::vector<T>, std::vector<int>, std::vector<int>>
    send_to_dest(const std::span<const T> data, const std::span<const int> dest)
    {
        SFEM_CHECK_SIZES(data.size(), dest.size());

        // Compute send buffer counts
        std ::vector<int> send_counts(n_procs(), 0);
        for (std::size_t i = 0; i < dest.size(); i++)
        {
            send_counts[dest[i]]++;
        }

        // Compute the send buffer displacements
        std::vector<int> send_displs(n_procs(), 0);
        std::exclusive_scan(send_counts.cbegin(), send_counts.cend(), send_displs.begin(), 0);

        // Fill the send buffer
        std::vector<T> send_buffer(std::accumulate(send_counts.cbegin(), send_counts.cend(), 0));
        std::fill(send_counts.begin(), send_counts.end(), 0);
        for (std::size_t i = 0; i < dest.size(); i++)
        {
            int owner = dest[i];
            send_buffer[send_displs[owner] + send_counts[owner]++] = data[i];
        }

        // Compute the receive buffer size
        std::vector<int> recv_counts(n_procs(), 0);
        int error_code = MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        SFEM_CHECK_MPI_ERROR(error_code);

        // Compute the recv buffer displacements
        std::vector<int> recv_displs(n_procs(), 0);
        std::exclusive_scan(recv_counts.cbegin(), recv_counts.cend(), recv_displs.begin(), 0);

        // Send the data to their owner process
        std::vector<T> recv_buffer(std::accumulate(recv_counts.cbegin(), recv_counts.cend(), 0));
        error_code = MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displs.data(), to_mpi_datatype<T>(),
                                   recv_buffer.data(), recv_counts.data(), recv_displs.data(), to_mpi_datatype<T>(),
                                   MPI_COMM_WORLD);
        SFEM_CHECK_MPI_ERROR(error_code);

        return {recv_buffer, recv_counts, recv_displs};
    }
    // Explicit instantiations
    // int
    template std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
        send_to_dest<int>(std::span<const int>, std::span<const int>);
    // real_t
    template std::tuple<std::vector<real_t>, std::vector<int>, std::vector<int>>
        send_to_dest<real_t>(std::span<const real_t>, std::span<const int>);
    //=============================================================================
    template <typename T>
    std::vector<T> distribute(const std::span<const T> data, const std::span<const int> dest)
    {
        SFEM_CHECK_SIZES(data.size(), dest.size());

        std ::vector<int> send_counts;
        std::vector<int> send_displs;

        if (rank() == root())
        {
            send_counts.resize(n_procs(), 0);
            send_displs.resize(n_procs(), 0);

            // Compute send buffer counts
            for (std::size_t i = 0; i < dest.size(); i++)
            {
                send_counts[dest[i]]++;
            }

            // Compute the send buffer displacements
            std::exclusive_scan(send_counts.cbegin(), send_counts.cend(), send_displs.begin(), 0);
        }

        // First send the size
        int n_recv;
        int error_code = MPI_Scatter(send_counts.data(), 1, MPI_INT, &n_recv, 1, MPI_INT, root(), MPI_COMM_WORLD);
        SFEM_CHECK_MPI_ERROR(error_code);

        // Receive buffer
        std::vector<T> recv_buffer(n_recv);
        error_code = MPI_Scatterv(data.data(), send_counts.data(),
                                  send_displs.data(), to_mpi_datatype<T>(),
                                  recv_buffer.data(), n_recv, to_mpi_datatype<T>(),
                                  root(), MPI_COMM_WORLD);
        SFEM_CHECK_MPI_ERROR(error_code);

        return recv_buffer;
    }
    // Explicit instantiations
    // int
    template std::vector<int> distribute<int>(const std::span<const int> data,
                                              const std::span<const int> dest);
    // real_t
    template std::vector<real_t> distribute<real_t>(const std::span<const real_t> data,
                                                    const std::span<const int> dest);
}