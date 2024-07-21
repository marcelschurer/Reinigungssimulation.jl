# 2D cannula flow simulation with an optional obstacle sphere at the beginning of the flow.

using TrixiParticles
using OrdinaryDiffEq

# ==========================================================================================
# ==== Resolution
particle_spacing = 0.05
smoothing_length = 3 * particle_spacing

# Make sure that the kernel support of fluid particles at a boundary is always fully sampled
boundary_layers = 4

# It is recommended to use `open_boundary_layers > boundary_layers`
open_boundary_layers = 6

# ==========================================================================================
# ==== Experiment Setup
tspan = (0.0, 20.0)
inflow_direction = [1.0, 0.0]
outflow_direction = [0.0, -1.0]

obstacle = true # true: an obstacle sphere is placed at the beginning of the flow

l = 30 # length of cannula to the outlet [mm]
d = 2 # diameter of the cannula [mm]
recess_length = 5 # length of the recess at the outlet [mm]
slip_wall_layers = 4 * open_boundary_layers

# note that open_boundary_layers <= stage_size_outlet_layers
number_of_stages_outlet = 28
stage_size_outlet_layers = 10

# note that you might have to change velocity_function_outlet(pos, t) for another prescribed_velocity
const prescribed_velocity = 1.0 # [m/s]
sound_speed = 10 * prescribed_velocity

# water at 20°C
fluid_density = 998.2 # [kg/m^3]
kinematic_viscosity = 1.004 * 10^-6 # [m^2/s]
reynolds_number = prescribed_velocity * d * 10^-3 / kinematic_viscosity

# For this particular example, it is necessary to have a background pressure.
# Otherwise the suction at the outflow is to big and the simulation becomes unstable.
pressure = 105000 # [Pa]

# ==========================================================================================
# ==== Initial Conditions

# = IC Fluid
fluid_size = (Int(floor((l + recess_length + 5) / particle_spacing)),
              Int(floor(d / particle_spacing)))

fluid_rectangular = RectangularShape(particle_spacing,
                                     (fluid_size),
                                     (0.0, 0.0),
                                     pressure=pressure,
                                     density=fluid_density,
                                     velocity=[0.0, 0.0])

fluid_sphere = SphereShape(particle_spacing, (d / 2),
                           ((l + recess_length + 5 + 0.5 * particle_spacing),
                            (d / 2)), fluid_density,
                           pressure=pressure, velocity=[0.0, 0.0],
                           sphere_type=RoundSphere(; start_angle=1.5π, end_angle=2.5π))

fluid_recess = RectangularShape(particle_spacing,
                                (Int(recess_length / particle_spacing),
                                 boundary_layers),
                                (l, -boundary_layers * particle_spacing),
                                pressure=pressure,
                                density=fluid_density,
                                velocity=[0.0, 0.0])

fluid_outlet = Vector{Any}(undef, number_of_stages_outlet)
for i in 1:number_of_stages_outlet
    if i == number_of_stages_outlet
        fluid_outlet[i] = RectangularShape(particle_spacing,
                                           (Int(recess_length / particle_spacing +
                                                (2 - 2 * i)),
                                            stage_size_outlet_layers - open_boundary_layers),
                                           ((l + (i - 1) *
                                                 particle_spacing),
                                            -(i * stage_size_outlet_layers +
                                              boundary_layers - open_boundary_layers) *
                                            particle_spacing), pressure=pressure,
                                           density=fluid_density)
    else
        fluid_outlet[i] = RectangularShape(particle_spacing,
                                           (Int(recess_length / particle_spacing +
                                                (2 - 2 * i)),
                                            stage_size_outlet_layers),
                                           ((l + (i - 1) *
                                                 particle_spacing),
                                            -(i * stage_size_outlet_layers +
                                              boundary_layers) *
                                            particle_spacing), pressure=pressure,
                                           density=fluid_density)
    end
end

obstacle_sphere = SphereShape(particle_spacing,
                              (3 * particle_spacing),
                              ((1 * slip_wall_layers * particle_spacing +
                                0.5 * particle_spacing),
                               (d / 2)), fluid_density,
                              pressure=pressure,
                              sphere_type=RoundSphere())

if obstacle == true
    ic_fluid = setdiff(union(fluid_rectangular, fluid_sphere, fluid_recess,
                             fluid_outlet...),
                       obstacle_sphere)
else
    ic_fluid = union(fluid_rectangular, fluid_sphere, fluid_recess, fluid_outlet...)
end

# = IC Boundary
boundary_size = (Int(floor((l + recess_length + 5) / particle_spacing +
                           open_boundary_layers - slip_wall_layers)),
                 Int(floor(d / particle_spacing + 2 * boundary_layers)))

boundary_rectangular = RectangularShape(particle_spacing,
                                        (boundary_size),
                                        ((slip_wall_layers - open_boundary_layers) *
                                         particle_spacing,
                                         -boundary_layers * particle_spacing),
                                        pressure=pressure,
                                        density=fluid_density)

boundary_sphere = SphereShape(particle_spacing,
                              ((d + 2 * boundary_layers * particle_spacing) / 2),
                              ((l + recess_length + 5 + 0.5 * particle_spacing),
                               (d / 2)), fluid_density,
                              pressure=pressure,
                              sphere_type=RoundSphere(; start_angle=1.5π, end_angle=2.5π))

ic_boundary = setdiff(union(boundary_rectangular, boundary_sphere), ic_fluid)

# = IC Slip Wall
boundary_slip_wall_inlet = RectangularShape(particle_spacing,
                                            (slip_wall_layers,
                                             boundary_size[2]),
                                            (-open_boundary_layers * particle_spacing,
                                             -boundary_layers * particle_spacing),
                                            pressure=pressure,
                                            density=fluid_density)

recess_slip_wall_inlet = RectangularShape(particle_spacing,
                                          (slip_wall_layers, fluid_size[2]),
                                          (-open_boundary_layers * particle_spacing,
                                           0.0),
                                          pressure=pressure,
                                          density=fluid_density)

slip_wall_outlet_left = Vector{Any}(undef, number_of_stages_outlet)
slip_wall_outlet_right = Vector{Any}(undef, number_of_stages_outlet)
for i in 1:number_of_stages_outlet
    slip_wall_outlet_left[i] = RectangularShape(particle_spacing,
                                                (boundary_layers,
                                                 stage_size_outlet_layers),
                                                ((l -
                                                  (boundary_layers - i + 1) *
                                                  particle_spacing),
                                                 -(i * stage_size_outlet_layers +
                                                   boundary_layers) *
                                                 particle_spacing), pressure=pressure,
                                                density=fluid_density)
    slip_wall_outlet_right[i] = RectangularShape(particle_spacing,
                                                 (boundary_layers,
                                                  stage_size_outlet_layers),
                                                 ((l + recess_length +
                                                   (1 - i) *
                                                   particle_spacing),
                                                  -(i * stage_size_outlet_layers +
                                                    boundary_layers) *
                                                  particle_spacing), pressure=pressure,
                                                 density=fluid_density)
end

ic_slip_wall = union(setdiff(boundary_slip_wall_inlet, recess_slip_wall_inlet),
                     slip_wall_outlet_left...,
                     slip_wall_outlet_right...)

# ==========================================================================================
# ==== Setup Fluid
smoothing_kernel = WendlandC2Kernel{2}()

fluid_density_calculator = ContinuityDensity()

n_buffer_particles = (boundary_layers + d)^2 * l^2

viscosity = ViscosityAdami(nu=kinematic_viscosity)

fluid_system = EntropicallyDampedSPHSystem(ic_fluid, smoothing_kernel,
                                           smoothing_length,
                                           sound_speed, viscosity=viscosity,
                                           transport_velocity=TransportVelocityAdami(pressure),
                                           density_calculator=fluid_density_calculator,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Setup Boundary

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
# ==== Setup Open Boundary
function velocity_function_inlet(pos, t)
    return SVector(prescribed_velocity, 0.0)
end

function velocity_function_outlet(pos, t)
    return SVector(0.0, -0.8 * prescribed_velocity)
    #return SVector(0.0, -0.7(1 + tanh(0.5(t - 5))))
end

inflow = InFlow(;
                plane=([0.0, 0.0],
                       [0.0, fluid_size[2] * particle_spacing]),
                flow_direction=inflow_direction,
                open_boundary_layers, density=fluid_density,
                particle_spacing)

inflow_system = OpenBoundarySPHSystem(inflow; sound_speed, fluid_system,
                                      buffer_size=n_buffer_particles,
                                      #   reference_pressure=pressure,
                                      reference_velocity=velocity_function_inlet,
                                      boundary_model=BoundaryModelTafuni())

outflow_plane_A = ((l + (number_of_stages_outlet - 1) * particle_spacing),
                   -(number_of_stages_outlet * stage_size_outlet_layers +
                     boundary_layers -
                     open_boundary_layers) * particle_spacing)
outflow_plane_B = (l + recess_length +
                   (1 - number_of_stages_outlet) * particle_spacing),
                  -(number_of_stages_outlet * stage_size_outlet_layers +
                    boundary_layers -
                    open_boundary_layers) * particle_spacing
outflow = OutFlow(;
                  plane=([outflow_plane_A[1], outflow_plane_A[2]],
                         [outflow_plane_B[1], outflow_plane_B[2]]),
                  flow_direction=outflow_direction,
                  open_boundary_layers,
                  density=fluid_density,
                  particle_spacing)

outflow_system = OpenBoundarySPHSystem(outflow; sound_speed, fluid_system,
                                       buffer_size=n_buffer_particles,
                                       reference_pressure=pressure,
                                       #reference_velocity=velocity_function_outlet,
                                       boundary_model=BoundaryModelTafuni())

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(fluid_system, boundary_system, slip_wall_system, inflow_system,
                          outflow_system)

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
