#include "connectivity.hpp"
#include "../base/error.hpp"
#include "../base/timer.hpp"
#include <numeric>
#include <format>

namespace sfem::graph
{
    //=============================================================================
    Connectivity::Connectivity(std::vector<int> &&offsets, std::vector<int> &&array)
        : offsets_(std::move(offsets)),
          array_(std::move(array)),
          n_secondary_(0)
    {
        SFEM_CHECK_SIZES(offsets_.back(), array_.size());

        // Compute the no. secondary entities
        if (array_.size() > 0)
        {
            n_secondary_ = std::ranges::max(array_) + 1;
        }

        // Check that the connectivity array is valid
        std::vector<bool> is_included(n_secondary_, false);
        for (std::size_t i = 0; i < array_.size(); i++)
        {
            if (array_[i] < 0)
            {
                break;
            }
            else
            {
                is_included[array_[i]] = true;
            }
        }
        if (!std::all_of(is_included.cbegin(),
                         is_included.cend(),
                         [](int e)
                         { return e; }))
        {
            SFEM_ERROR(std::format("Invalid connectivity array\n"));
        }
    }
    //=============================================================================
    std::vector<int> Connectivity::offsets() const
    {
        return offsets_;
    }
    //=============================================================================
    std::vector<int> Connectivity::array() const
    {
        return array_;
    }
    //=============================================================================
    int Connectivity::n_primary() const
    {
        return static_cast<int>(offsets_.size()) - 1;
    }
    //=============================================================================
    int Connectivity::n_secondary() const
    {
        return n_secondary_;
    }
    //=============================================================================
    int Connectivity::n_links() const
    {
        return static_cast<int>(array_.size());
    }
    //=============================================================================
    int Connectivity::n_links(int primary) const
    {
        SFEM_CHECK_INDEX(primary, n_primary());
        return offsets_[primary + 1] - offsets_[primary];
    }
    //=============================================================================
    std::span<const int> Connectivity::links(int primary) const
    {
        return std::span<const int>(array_.data() + offsets_[primary], n_links(primary));
    }
    //=============================================================================
    int Connectivity::offset(int primary) const
    {
        SFEM_CHECK_INDEX(primary, n_primary() + 1);
        return offsets_[primary];
    }
    //=============================================================================
    int Connectivity::relative_index(int primary, int secondary) const
    {
        auto links_ = links(primary);
        auto it = std::find(links_.cbegin(),
                            links_.cend(),
                            secondary);
        if (it == links_.end())
        {
            SFEM_ERROR(std::format("{} is not a link of {}\n", secondary, primary));
            return -1;
        }
        return static_cast<int>(std::distance(links_.begin(), it));
    }
    //=============================================================================
    Connectivity Connectivity::invert() const
    {
        // Compute the inverse connectivity offsets
        std::vector<int> inverse_n_links(n_secondary(), 0);
        for (int i : array_)
        {
            inverse_n_links[i]++;
        }
        std::vector<int> inverse_offsets(n_secondary() + 1, 0);
        std::inclusive_scan(inverse_n_links.cbegin(),
                            inverse_n_links.cend(),
                            inverse_offsets.begin() + 1);

        // Fill the inverse connectivity array
        std::vector<int> inverse_array(inverse_offsets.back(), 0);
        std::fill(inverse_n_links.begin(), inverse_n_links.end(), 0);
        for (int i = 0; i < n_primary(); i++)
        {
            for (int j : links(i))
            {
                inverse_array[inverse_offsets[j] + inverse_n_links[j]++] = i;
            }
        }

        return Connectivity(std::move(inverse_offsets),
                            std::move(inverse_array));
    }
    //=============================================================================
    Connectivity Connectivity::primary_to_primary(int n_common,
                                                  bool include_self) const
    {
        Timer timer;

        // Inverse connectivity
        auto inverse = invert();

        // Primary-to-primary connectivity offsets
        std::vector<int> ptp_offsets(n_primary() + 1, 0);

        // Primary-to-primary connectivity array
        std::vector<int> ptp_array;

        // Create the primary-to-primary connectivity
        // First call computes the offsets
        // Second call fills the connectivity array
        auto create_conn = [&](bool first_call)
        {
            // Check if two primaries are linked
            // Returns true if the two primaries share more than n_common secondaries
            auto compare = [n_common](std::span<const int> lhs, std::span<const int> rhs)
            {
                int _n_common = 0;
                for (int i : lhs)
                {
                    for (int j : rhs)
                    {
                        if (i == j)
                        {
                            if (++_n_common == n_common)
                            {
                                return true;
                            }
                        }
                    }
                }

                return false;
            };

            // Map to check if another primary has already been registered as a link
            std::unordered_map<int, bool> is_link;

            // Number of primary-to-primary connections per primary
            std::vector<int> ptp_n_links(n_primary(), 0);

            // The loop works by iterating first over all the primaries (i),
            // then over each secondary (j) of the primary, and then over all
            // primaries linked to each secondary (k).
            // Then, if i and k share at least n_common secondaries,
            // k is added as a link of i
            for (int i = 0; i < n_primary(); i++)
            {
                for (int j : links(i))
                {
                    for (int k : inverse.links(j))
                    {
                        if (i == k and include_self == false)
                        {
                            continue;
                        }

                        if (is_link.contains(k) == false and (n_common == 1 or compare(links(i), links(k))))
                        {
                            if (!first_call)
                            {
                                ptp_array[ptp_offsets[i] + ptp_n_links[i]] = k;
                            }
                            ptp_n_links[i]++;
                            is_link[k] = true;
                        }
                    }
                }

                // Reset
                is_link.clear();
            }

            if (first_call)
            {
                // Compute the primary-primary connectivity offsets
                std::inclusive_scan(ptp_n_links.cbegin(),
                                    ptp_n_links.cend(),
                                    ptp_offsets.begin() + 1);

                // Initialize the primary-to-primary connectivity array
                ptp_array.resize(ptp_offsets.back(), 0);
            }
        };

        // First loop, count ptp connections
        create_conn(true);

        // Second loop, fill ptp connectivity array
        create_conn(false);

        return Connectivity(std::move(ptp_offsets),
                            std::move(ptp_array));
    }
    //=============================================================================
    std::string Connectivity::str(std::string_view name) const
    {
        auto repr = std::format("{}\n", name);
        repr += std::format("No. primary entities: {}\n", n_primary());
        repr += std::format("No. secondary entities: {}\n", n_secondary());
        repr += std::format("No. connections: {}\n", n_links());
        repr += "Connections:\n";
        for (int i = 0; i < n_primary(); i++)
        {
            repr += std::format("\t{} -> ", i);
            for (int j : links(i))
            {
                repr += std::format("{} ", j);
            }
            repr += "\n";
        }
        return repr;
    }
}