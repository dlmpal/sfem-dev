// Solve the incompressible Navier-Stokes equations using SIMPLE

#include "simple.hpp"

int main(int argc, char *argv[])
{
    initialize(argc, argv, false, "SIMPLESolver");

    // SIMPLE solver instance and options
    SimpleSolver solver = from_file(argv[1]);

    const real_t time_start = 0.0;
    const real_t time_stop = 1.0;
    real_t dt = time_stop;
    real_t time = time_start;

    while (time < time_stop)
    {
        solver.pre_timestep();
        solver.step(time, dt);
        solver.post_timestep();
        time += dt;
    }

    return 0;
}