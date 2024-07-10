# 2D backward facing step simulation following:
# http://dx.doi.org/10.1017/S0022112083002839
# https://doi.org/10.1016/j.icheatmasstransfer.2011.07.003

using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.1
smoothing_length = 3 * particle_spacing

# Make sure that the kernel support of fluid particles at a boundary is always fully sampled
boundary_layers = 5

# It is recommended to use `open_boundary_layers > boundary_layers`
open_boundary_layers = 6

# ==========================================================================================
# ==== Experiment Setup
tspan = (0.0, 5.0)
inflow_direction = [1.0, 0.0]
outflow_direction = [0.0, -1.0]

# note that you might have to change velocity_function_outlet(pos, t) for another prescribed_velocity
const prescribed_velocity = 1.0
sound_speed = 10 * prescribed_velocity

reynolds_number = 389 # laminar flow
fluid_density = 997.0 # water at 20°C
kinematic_viscosity = 0.001 # water at 20°C

# For this particular example, it is necessary to have a background pressure.
# Otherwise the suction at the outflow is to big and the simulation becomes unstable.
pressure = 101325 - 0.5 * fluid_density * prescribed_velocity^2 # Bernoulli's equation

l = 30 # length of cannula to the outlet
d = 2 # diameter of the cannula
recess_length = 5 # length of the recess at the outlet
slip_wall_length = 2 * open_boundary_layers

# ==========================================================================================
# ==== Initial Conditions

fluid_size = (Int(floor((l + recess_length + 5) / particle_spacing))),
             Int(floor(d / particle_spacing))

fluid_rectangular = RectangularShape(particle_spacing,
                                     (fluid_size),
                                     (0.0, 0.0),
                                     pressure=pressure,
                                     density=fluid_density,
                                     velocity=[prescribed_velocity, 0.0])

fluid_sphere = SphereShape(particle_spacing, (d / 2),
                           ((l + recess_length + 5),
                            (d / 2)), fluid_density,
                           pressure=pressure, sphere_type=RoundSphere())

ic_fluid = union(fluid_rectangular, fluid_sphere)

boundary_size = (Int(floor((l + recess_length + 5) / particle_spacing +
                           open_boundary_layers -
                           slip_wall_length))),
                Int(floor(d / particle_spacing + 2 * boundary_layers))

boundary_rectangular = RectangularShape(particle_spacing,
                                        (boundary_size),
                                        ((slip_wall_length - open_boundary_layers) *
                                         particle_spacing,
                                         -boundary_layers * particle_spacing),
                                        pressure=pressure,
                                        density=fluid_density)

boundary_sphere = SphereShape(particle_spacing,
                              ((d + 2 * boundary_layers * particle_spacing) / 2),
                              ((l + recess_length + 5),
                               (d / 2)), fluid_density,
                              pressure=pressure, sphere_type=RoundSphere())

cut_out_outlet_size = (Int(floor(recess_length / particle_spacing))),
                      Int(floor(boundary_layers))

cut_out_outlet = RectangularShape(particle_spacing,
                                  cut_out_outlet_size,
                                  (l,
                                   -boundary_layers * particle_spacing),
                                  pressure=pressure,
                                  density=fluid_density)

ic_boundary = setdiff(union(boundary_rectangular, boundary_sphere), ic_fluid,
                      cut_out_outlet)

boundary_slip_wall = RectangularShape(particle_spacing,
                                      (slip_wall_length,
                                       boundary_size[2]),
                                      (-open_boundary_layers * particle_spacing,
                                       -boundary_layers * particle_spacing),
                                      pressure=pressure,
                                      density=fluid_density)

cut_out_slip_wall = RectangularShape(particle_spacing,
                                     (slip_wall_length, fluid_size[2]),
                                     (-open_boundary_layers * particle_spacing,
                                      0.0),
                                     pressure=pressure,
                                     density=fluid_density)

ic_slip_wall = setdiff(boundary_slip_wall, cut_out_slip_wall)

# ==========================================================================================
# ==== Fluid
smoothing_kernel = WendlandC2Kernel{2}()

fluid_density_calculator = ContinuityDensity()

n_buffer_particles = open_boundary_layers * l * d

viscosity = ViscosityAdami(nu=kinematic_viscosity)

fluid_system = EntropicallyDampedSPHSystem(ic_fluid, smoothing_kernel,
                                           smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           transport_velocity=TransportVelocityAdami(pressure),
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Open Boundary
function velocity_function_inlet(pos, t)
    return SVector(prescribed_velocity, 0.0)
end

function velocity_function_outlet(pos, t)
    return SVector(0.0, -prescribed_velocity)
    # time dependent velocity profile at the outlet for stability
    # f(t) is designed for prescribed_velocity = 1.0
    # return SVector(0.15(1 + tanh(5(t - 0.9))), 0)
end

inflow = InFlow(;
                plane=([0.0, 0.0],
                       [0.0, fluid_size[2] * particle_spacing]),
                flow_direction=inflow_direction,
                open_boundary_layers, density=fluid_density,
                particle_spacing)

open_boundary_in = OpenBoundarySPHSystem(inflow; sound_speed, fluid_system,
                                         buffer_size=n_buffer_particles,
                                         reference_pressure=pressure,
                                         reference_velocity=velocity_function_inlet,
                                         boundary_model=BoundaryModelTafuni())

outflow = OutFlow(;
                  plane=([l, 0.0],
                         [l + recess_length, 0.0]),
                  flow_direction=outflow_direction,
                  open_boundary_layers=(open_boundary_layers - 2),
                  density=fluid_density,
                  particle_spacing)

open_boundary_out = OpenBoundarySPHSystem(outflow; sound_speed, fluid_system,
                                          buffer_size=n_buffer_particles,
                                          reference_pressure=pressure,
                                          reference_velocity=velocity_function_outlet,
                                          boundary_model=BoundaryModelTafuni())

# ==========================================================================================
# ==== Boundary

boundary_model = BoundaryModelDummyParticles(ic_boundary.density,
                                             ic_boundary.mass,
                                             AdamiPressureExtrapolation(),
                                             viscosity=viscosity,
                                             smoothing_kernel, smoothing_length)

boundary_system = BoundarySPHSystem(ic_boundary, boundary_model)

slip_wall_model = BoundaryModelDummyParticles(ic_slip_wall.density,
                                              ic_slip_wall.mass,
                                              AdamiPressureExtrapolation(),
                                              viscosity=nothing, # slip wall condition
                                              smoothing_kernel, smoothing_length)

slip_wall_system = BoundarySPHSystem(ic_slip_wall, slip_wall_model)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(fluid_system, open_boundary_in, open_boundary_out,
                          boundary_system, slip_wall_system)

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)
saving_callback = SolutionSavingCallback(dt=0.02, prefix="",
                                         output_directory="out")

callbacks = CallbackSet(info_callback, saving_callback, UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-5, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-3, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
