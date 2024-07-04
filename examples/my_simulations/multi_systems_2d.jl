using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.02

# ==========================================================================================
# ==== Experiment Setup
gravity = 9.81
tspan = (0.0, 1.0)
fluid_flow = 0.5  # spawning time of each fluid_system

number_of_fluid_systems = Int(tspan[2] / fluid_flow) # TODO: Wenn tspan[2] nicht durch fluid_flow teilbar, dann?

boundary_layers = 3
tank_size = (200, 60)
fluid_size = (2, tank_size[2] - 2 * boundary_layers)

fluid_density = 1000.0
sound_speed = 10.0

# ==========================================================================================
# ==== Boundary geometry and initial fluid particle positions

state_equation = StateEquationCole(; sound_speed, reference_density=fluid_density,
                                   exponent=7, clip_negative_pressure=false)

tank_bottom = RectangularShape(particle_spacing, (tank_size[1], boundary_layers),
                               (0.0, 0.0), density=fluid_density)

tank_top = RectangularShape(particle_spacing, (tank_size[1], boundary_layers),
                            (0.0, (tank_size[2] - boundary_layers) * particle_spacing),
                            density=fluid_density)
tank_left = RectangularShape(particle_spacing,
                             (boundary_layers, tank_size[2] - 2 * boundary_layers),
                             (0.0, boundary_layers * particle_spacing),
                             density=fluid_density)
tank_right = RectangularShape(particle_spacing,
                              (boundary_layers, tank_size[2] - 2 * boundary_layers),
                              ((tank_size[1] - boundary_layers) * particle_spacing,
                               boundary_layers * particle_spacing),
                              density=fluid_density)
tank = union(tank_bottom, tank_top, tank_left, tank_right)

fluid = RectangularShape(particle_spacing, fluid_size,
                         (boundary_layers * particle_spacing,
                          boundary_layers * particle_spacing),
                         density=fluid_density, velocity=(1.0, 0.0))

# ==========================================================================================
# ==== Fluid
smoothing_length = 1.2 * particle_spacing
smoothing_kernel = SchoenbergCubicSplineKernel{2}()

viscosity = ArtificialViscosityMonaghan(alpha=0.02, beta=0.0)

fluid_density_calculator = ContinuityDensity()

fluid_systems = Vector{Any}(undef, number_of_fluid_systems)
for i in 1:number_of_fluid_systems
    fluid_systems[i] = WeaklyCompressibleSPHSystem(fluid,
                                                   fluid_density_calculator,
                                                   state_equation,
                                                   smoothing_kernel,
                                                   smoothing_length,
                                                   viscosity=viscosity,
                                                   acceleration=(0.0, -gravity),  # -gravity deleted
                                                   source_terms=nothing)
end

# ==========================================================================================
# ==== Boundary

# This is to set another boundary density calculation with `trixi_include`
boundary_density_calculator = AdamiPressureExtrapolation()

# This is to set wall viscosity with `trixi_include`
viscosity_wall = nothing
boundary_model = BoundaryModelDummyParticles(tank.density,
                                             tank.mass,
                                             boundary_density_calculator, smoothing_kernel,
                                             smoothing_length; viscosity=nothing,
                                             state_equation=state_equation,
                                             correction=nothing)

boundary_system = BoundarySPHSystem(tank, boundary_model,
                                    movement=nothing)

# ==========================================================================================
# ==== Simulation
systems = (boundary_system, fluid_systems...)

semi = Semidiscretization(systems...)
ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=fluid_flow, prefix="",
                                         output_directory="out_multi_systems_2d")

# This is to easily add a new callback with `trixi_include`
extra_callback = nothing

callbacks = CallbackSet(info_callback, saving_callback, extra_callback)

# Use a Runge-Kutta method with automatic (error based) time step size control
sol = solve(ode, RDPK3SpFSAL35(), save_everystep=false, callback=callbacks);
