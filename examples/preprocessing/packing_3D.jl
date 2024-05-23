using TrixiParticles
using OrdinaryDiffEq

particle_spacing = 0.02
density = 1000.0
tspan = (0, 10.0)

dir = joinpath("Data", "stl-files", "examples")
filename = joinpath(expanduser("~/") * dir, "aorta.stl")

# Returns `Shape`
shape = load_shape(filename)

# Returns `InitialCondition`.
shape_sampled = ComplexShape(shape; particle_spacing, density, hierarchical_winding=true)

signed_distance_field = SignedDistanceField(shape, particle_spacing;
                                            max_signed_distance=4particle_spacing,
                                            neighborhood_search=true)

background_pressure = 1e8 * particle_spacing^3

packing_system = ParticlePackingSystem(shape_sampled; tlsph=true,
                                       signed_distance_field, neighborhood_search=true,
                                       boundary=shape, background_pressure)

boundary_system = ParticlePackingSystem(shape_sampled; tlsph=true,
                                        signed_distance_field,
                                        is_boundary=true, neighborhood_search=true,
                                        boundary=shape, background_pressure)

semi = Semidiscretization(packing_system)

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="")

callbacks = CallbackSet(UpdateCallback(), saving_callback, info_callback)

sol = solve(ode, RK4();
            save_everystep=false, maxiters=100, callback=callbacks, dtmax=1e-2)

packed_ic = InitialCondition(sol, packing_system, semi)
#packed_boundary_ic = InitialCondition(sol, boundary_system, semi);
