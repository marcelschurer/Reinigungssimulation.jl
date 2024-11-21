using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution

tspan = (0.0, 0.1)
particle_spacing = 3.0 * 10^-3 # 1 mm
smoothing_length = 3.0 * particle_spacing

# ==========================================================================================
# ==== Experiment Setup for air at 20°C
reynolds_number = 292685 # [-]
kinematic_viscosity = 1.5 * 10^-5 # [m^2/s]

prescribed_velocity = 55.0 # [m/s]
density = 1.204 # [kg/m^3]
pressure = 101325.0 # [Pa]

sound_speed = 10 * prescribed_velocity

# ==========================================================================================
# ==== Inlet

filename = "inlet"
file = joinpath("examples", "preprocessing", "data", "intake_manifold", filename * ".stl")

geometry = load_geometry(file)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry,
                                                    # winding_number_factor=0.4,
                                                    hierarchical_winding=true)

# Returns `InitialCondition`
inlet_sampled = ComplexShape(geometry; particle_spacing, density=density,
                             boundary_thickness=5 * particle_spacing,
                             create_signed_distance_field=true,
                             sample_boundary=false,
                             point_in_geometry_algorithm)

#trixi2vtk(inlet_sampled, filename="inlet")

# ==========================================================================================
# ==== Intake Manifold

filename = "intake_manifold"
file = joinpath("examples", "preprocessing", "data", "intake_manifold", filename * ".stl")

geometry = load_geometry(file)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry,
                                                    #winding_number_factor=0.4,
                                                    hierarchical_winding=true)

# Returns `InitialCondition`
intake_manifold_sampled = ComplexShape(geometry; particle_spacing, density=density,
                                       point_in_geometry_algorithm)

#trixi2vtk(intake_manifold_sampled, filename="intake_manifold")

# ==========================================================================================
# ==== Outlet Left

filename = "outlet_left"
file = joinpath("examples", "preprocessing", "data", "intake_manifold", filename * ".stl")

geometry = load_geometry(file)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry,
                                                    # winding_number_factor=0.4,
                                                    hierarchical_winding=true)

# Returns `InitialCondition`
outlet_left_sampled = ComplexShape(geometry; particle_spacing, density=density,
                                   point_in_geometry_algorithm)

#trixi2vtk(outlet_left_sampled, filename="outlet_left")

# ==========================================================================================
# ==== Outlet Right

filename = "outlet_right"
file = joinpath("examples", "preprocessing", "data", "intake_manifold", filename * ".stl")

geometry = load_geometry(file)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry,
                                                    # winding_number_factor=0.4,
                                                    hierarchical_winding=true)

# Returns `InitialCondition`
outlet_right_sampled = ComplexShape(geometry; particle_spacing, density=density,
                                    point_in_geometry_algorithm)

#trixi2vtk(outlet_right_sampled, filename="outlet_right")

# ==========================================================================================
# ==== Fluid
smoothing_kernel = WendlandC2Kernel{3}()

fluid_density_calculator = ContinuityDensity()

viscosity = ViscosityAdami(nu=kinematic_viscosity)

n_buffer_particles = 1 * size(intake_manifold_sampled.coordinates, 2)

fluid_system = EntropicallyDampedSPHSystem(intake_manifold_sampled, smoothing_kernel,
                                           smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Boundary
boundary_system = ParticlePackingSystem(intake_manifold_sampled; tlsph=tlsph,
                                        is_boundary=true, background_pressure=pressure)

packing_system = ParticlePackingSystem(intake_manifold_sampled; tlsph=tlsph,
                                       background_pressure=pressure)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(packing_system, boundary_system)

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=100)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="")

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-6, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);