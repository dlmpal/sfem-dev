#include "fe_function.hpp"

namespace sfem::fem
{
    //=============================================================================
    FEFunction::FEFunction(std::shared_ptr<const FESpace> fe_space,
                           const std::vector<std::string> &components)
        : Function(fe_space->index_map(), components),
          fe_space_(fe_space)
    {
        if (components_.size() <= 0)
        {
            SFEM_ERROR(std::format("Cannot create FESpace with {} (<=0) components\n",
                                   components_.size()));
        }
    }
    //=============================================================================
    std::shared_ptr<const FESpace> FEFunction::space() const
    {
        return fe_space_;
    }
    //=============================================================================
    la::DenseMatrix FEFunction::cell_values(int cell_idx) const
    {
        const auto cell_dof = fe_space_->cell_dof(cell_idx);
        const int n_dof = static_cast<int>(cell_dof.size());
        la::DenseMatrix cell_values(n_dof * block_size_, 1);
        for (int i = 0; i < n_dof; i++)
        {
            for (int j = 0; j < block_size_; j++)
            {
                cell_values(i * block_size_ + j, 0) = values_[cell_dof[i] * block_size_ + j];
            }
        }
        return cell_values;
    }
}