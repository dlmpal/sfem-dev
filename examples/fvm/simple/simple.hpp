#include "sfem.hpp"

using namespace sfem;
using namespace sfem::fvm;

std::pair<std::vector<FVField>, FVField> create_fields(std::shared_ptr<const FVSpace> V)
{
    const std::array<std::string, 3> U_names = {"u", "v", "w"};
    std::vector<FVField> U;
    for (int i = 0; i < V->mesh()->pdim(); i++)
    {
        U.push_back(FVField(V, {U_names[i]}, GradientMethod::green_gauss));
    }
    FVField P(V, {"P"}, GradientMethod::green_gauss);
    return {U, P};
}