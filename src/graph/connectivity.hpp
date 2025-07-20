#pragma once

#include <vector>
#include <algorithm>
#include <functional>
#include <span>
#include <string>

namespace sfem::graph
{
    /// This class implements the functionality required to describe the connectivity
    /// between a set of primary entities and a set of secondary entities.
    /// For example, consider a (primary) set of two cells, named C1 and C2,
    /// and a (secondary) set of six vertices, named V1 through V6.
    /// If C1 is comprised of V1, V2, V4 and V5 and C2 of V2, V3, V5, V6,
    /// then an instance this class will store this information as follows:
    //
    /// offsets = [0 | 4 | 8]
    /// array = [V1, V2, V4, V5 | V2, V3, V5, V6]
    //
    /// The connectivity relation can also be "inverted", so instead of
    /// storing the vertices of each cell, it is possible to store
    /// the cells of each vertex. For the above example, this would
    /// correspond to the following:
    ///
    /// inverse_offsets = [0, 1, 3, 4, 5, 7, 8]
    /// inverse_array = [C1 | C1, C2 | C2 | C1 | C1, C2 | C2]
    ///
    /// @note The indices of the secondary entities stored in array (or array_) must belong
    /// to the contiguous range [0, n), where n is the number of secondary entities.
    class Connectivity
    {
    public:
        /// @brief Create a Connectivity
        /// @param offsets Offsets
        /// @param array Connectivity array
        Connectivity(std::vector<int> &&offsets = {0}, std::vector<int> &&array = {});

        /// @brief Get the offsets
        std::vector<int> offsets() const;

        /// @brief Get the connectivity array
        std::vector<int> array() const;

        /// @brief Get the number of primary entities
        int n_primary() const;

        /// @brief Get the number of secondary entities
        int n_secondary() const;

        /// @brief Get the number of total connections
        int n_links() const;

        /// @brief Get the number of secondary entities connected to a primary
        int n_links(int primary) const;

        /// @brief Get the secondary entities connected to a primary
        std::span<const int> links(int primary) const;

        /// @brief Get the offset for a primary
        int offset(int primary) const;

        /// @brief Get the relative index of a secondary from
        /// the perspective of a primary. For example, consider
        /// that the links of 5 are 10, 11, 12, 13.
        /// Then the relative index of 12 in 5 is 2
        int relative_index(int primary, int secondary) const;

        /// @brief Compute the inverse connectivity
        Connectivity invert() const;

        /// @brief Compute the connectivity between primaries
        /// @param n_common Number of required common secondaries
        /// for two primaries to be considered linked
        /// @note n_common must be greater than 0, else an empty
        /// Connectivity is returned
        Connectivity primary_to_primary(int n_common = 1, bool include_self = false) const;

        /// @brief Get a string representation of the connectivity
        std::string str(std::string_view name = "Connectivity") const;

    private:
        /// @brief The connectivity offsets
        std::vector<int> offsets_;

        /// @brief The connectivity array
        std::vector<int> array_;

        /// @brief The number of secondary entities
        int n_secondary_;
    };
}