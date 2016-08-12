##### Interfaces #####

# We must convert ydot into an array and assign element-by-element
# This is because we are modifying by-reference an NVector in Sundials...
function interface_ode(t::Sundials.realtype, y::Sundials.N_Vector, ydot::Sundials.N_Vector, user_data::Vector{Any})
  y = Sundials.asarray(y)
  ydot = Sundials.asarray(ydot)
  user_data[1](t, y, user_data[2], user_data[3], ydot)
  return int32(0)
end

##### Simulation functions #####

# Function to simulate the model
@inbounds function simulate(times::Vector{Float64},
                            states::Vector{Float64},
                            parameters::Vector{Float64},
                            inputs::Vector{Function},
                            settings::Dict{ASCIIString,Any},
                            model::Function,
                            observer::Function)
  neq = length(states)

  if settings["method"] == "bdf"
    cvode_mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
  elseif settings["method"] == "adams"
    cvode_mem = Sundials.CVodeCreate(Sundials.CV_ADAMS, Sundials.CV_NEWTON)
  else
    error("Please choose bdf or adams as method")
  end
# Shut up Sundials
  flag = Sundials.CVodeSetErrFile(cvode_mem, C_NULL )
  if flag < zero(flag)
    error("Error in the CVodeSetErrFile function")
  end

  # Initialize the Sundials solver by passing the wrapper function
  flag = Sundials.CVodeInit(cvode_mem, cfunction(interface_ode, Int32,
    (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector,  Vector{Any})),
             times[1], Sundials.nvector(states))
  if flag < zero(flag)
    error("Error in the CVodeInit function")
  end

  # Tell Sundials the tolerance settings for error control
  flag = Sundials.CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"])
  if flag < zero(flag)
    error("Error in the CVodeSStolerances function")
  end

  # Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
  flag = Sundials.CVDense(cvode_mem, neq)
  if flag < zero(flag)
    error("Error in the CVDense function")
  end

  # Pass information for the callback
  flag = Sundials.CVodeSetUserData(cvode_mem,{model, parameters, inputs, observed})
  if flag < zero(flag)
    error("Error in the CVodeSetUserData function")
  end

  # Set maximum number of steps
  flag = Sundials.CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxNumSteps function")
  end

  # Set maximum order of the integration
  flag = Sundials.CVodeSetMaxOrd(cvode_mem, settings["maxord"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxOrd function")
  end

  # Set the initial step size
  flag = Sundials.CVodeSetInitStep(cvode_mem, settings["hini"])
  if flag < zero(flag)
    error("Error in the CVodeSetInitStep function")
  end

  # Set the minimum step size
  flag = Sundials.CVodeSetMinStep(cvode_mem, settings["hmin"])
  if flag < zero(flag)
    error("Error in the CVodeSetMinStep function")
  end

  # Set the maximum step size
  flag = Sundials.CVodeSetMaxStep(cvode_mem, settings["hmax"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxStep function")
  end

  # Set the maximum number of error test fails
  flag = Sundials.CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxErrTestFails function")
  end

  # Set the maximum number of nonlinear iterations per step
  flag = Sundials.CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxNonlinIters function")
  end

  # Set the maximum number of convergence failures
  flag = Sundials.CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"])
  if flag < zero(flag)
    error("Error in the CVodeSetMaxConvFails function")
  end

  # Make a first call to the model to check that everything is ok and retrieve the number of observed variables
  observed = zeros(Float64, settings["nobs"])
  observer(times[1], states, parameters, inputs, observed)

  # Use time as column -> Because of column-ordered
  output = zeros(Float64, (neq + settings["nobs"] + 1, length(times)))
  output[:,1] = [times[1], states, observed]

  # Main time loop. Each timestep call cvode. Handle exceptions and fill up output
  t = times[1]
  y = copy(states)
  for i in 2:length(times)
    # Note that the second t must be an array to build the pointer...
    flag = Sundials.CVode(cvode_mem, times[i], y, [t], Sundials.CV_NORMAL)
    if flag < zero(flag)
       if flag == oftype(flag, -1)
          error("The solver took mxstep internal steps but could not reach tout.")
       elseif flag == oftype(flag, -2)
          error("The solver could not satisfy the accuracy demanded by the user for some internal step")
       elseif flag == oftype(flag, -3)
          error("Error test failures occured too many times during one internal time step or minimum step size was reached")
       elseif flag == oftype(flag, -4)
          error("Convergence test failures occurred too many times during one internal time step or minimum step size was reached.")
       elseif flag == oftype(flag, -5)
          error("The linear solver’s initialization function failed.")
       elseif flag == oftype(flag, -6)
          error("The linear solver’s setup function failed in an unrecoverable manner")
       elseif flag == oftype(flag, -7)
          error("The linear solver’s solve function failed in an unrecoverable manner")
       elseif flag == oftype(flag, -8)
          error("The right hand side function failed in an unrecoverable manner")
       elseif flag == oftype(flag, -9)
          error("The right-hand side function failed at the first call.")
       elseif flag == oftype(flag, -10)
          error("The right-hand side function had repeated recoverable errors.")
       elseif flag == oftype(flag, -11)
          error("The right-hand side function had a recoverable errors but no recovery is possible.")
       elseif flag == oftype(flag, -25)
          error("The time t is outside the last step taken.")
       elseif flag == oftype(flag, -26)
          error("The output derivative vector is NULL.")
       elseif flag == oftype(flag, -27)
          error("The output and initial times are too close to each other.")
       end
     end
    output[1,i] = times[i]
    for h in 1:neq
      output[h+1,i] = y[h]
    end
  end
  # If we have observed variables we call the model function again
  if settings["nobs"] > 0
    for i in 1:length(times)
      observer(times[i], vec(output[2:(neq + 1), i]), parameters, inputs, observed)
      for h in 1:settings["nobs"]
        output[neq + 1 + h, i] = observed[h]
      end
    end
  end

  return transpose(output)
end

## Return the solution of the model at the specific timepoints
function simulate(model::OdeModel, time::Vector{Float64}, settings::Dict{ASCIIString, Any})
    simulation = simulate(time, convert(Array, model.States), convert(Array, model.Parameters),
                            convert(Array, model.Inputs), settings, model.Model, model.Observer)
end



# Generate Dict with default settings. 

function CvodeSettings()
  settings = Dict{ASCIIString, Any}()
  settings["method"] = "bdf"
  settings["rtol"] = 1e-6
  settings["atol"] = 1e-6
  settings["maxsteps"] = 1e5
  settings["maxord"] = 5
  settings["hini"] = 0
  settings["hmin"] = 0
  settings["hmax"] = 0
  settings["maxerr"] = 12
  settings["maxnonlin"] = 12
  settings["maxconvfail"] = 12
  settings["jacobian"] = false
  return settings
end