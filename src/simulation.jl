using Sundials
using StrPack

# Add a struct for calculations with analytical Jacobian

 @struct type J_DlsMat
   typ::Int32
   M::Int
   N::Int
   ldim::Int
   mu::Int
   ml::Int
   s_mu::Int
   data::Ptr{Float64}
   ldata::Int
   cols::Ptr{Ptr{Float64}}
end

@inbounds function linear_interpolate(x::Array{Float64,1}, y::Array{Float64,1}, xout::Float64)
   if xout >= x[end]
      return y[end]
   elseif xout <= x[1]
      return y[1]
   else
      index = searchsortedfirst(x, xout) - 1
      inty = (y[index + 1] - y[index])/(x[index + 1] - x[index])*(xout -x[index]) + y[index]
      return inty
   end
end

# We must convert ydot into an array and assign element-by-element
# This is because we are modifying by-reference an NVector in Sundials...
@inbounds function interface_ode(t::Sundials.realtype, y::Sundials.N_Vector, ydot::Sundials.N_Vector, user_data::Array{Any, 1})
    y = Sundials.asarray(y)
    ydot = Sundials.asarray(ydot)
    forcings = [linear_interpolate(user_data[3][i][:,1], user_data[3][i][:,2], t) for i in 1:length(user_data[3])]
    derivatives = user_data[1](t, y, user_data[2], forcings)[1]
    for i in 1:length(derivatives)
        ydot[i] = derivatives[i]
    end
    return int32(0)
end

# Function to fill in the Jacobian matrix (this increases memory allocation)
# The Jacobian function modifies in-place and it is stored as user_data[4]
function interface_Jac(N, t, y, fy, Jptr, user_data_p::Ptr{Void},tmp1, tmp2, tmp3)
y = Sundials.asarray(y)
neq = length(y)
# Don't really understand this...
dlsmat = unpack(IOBuffer(pointer_to_array(convert(Ptr{Uint8}, Jptr),
               (sum(map(sizeof, J_DlsMat.types))+10,))),J_DlsMat)
J = pointer_to_array(unsafe_load(dlsmat.cols, 1), (int(neq), int(neq)), false)
# Need to convert user_data into Array{Any, 1}...
user_data = unsafe_pointer_to_objref(user_data_p)
# Do the annoying interpolation again...
forcings = [linear_interpolate(user_data[3][i][:,1], user_data[3][i][:,2], t) for i in 1:length(user_data[3])]
# Call Jacobian function and modify in-place
user_data[4](t, y, user_data[2], forcings, J)
return int32(0)
end

# Function to simulate the model
@inbounds function simulate(times, states, parameters, forcings_data, settings, model, jacobian)
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
  if flag < 0.0
    error("Error in the CVodeSetErrFile function")
  end

  # Initialize the Sundials solver by passing the wrapper function
  flag = Sundials.CVodeInit(cvode_mem, cfunction(interface_ode, Int32,
    (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector,  Array{Any, 1})),
             times[1], Sundials.nvector(states))
  if flag < 0.0
    error("Error in the CVodeInit function")
  end

  # Tell Sundials the tolerance settings for error control
  flag = Sundials.CVodeSStolerances(cvode_mem, settings["rtol"], settings["atol"])
  if flag < 0.0
    error("Error in the CVodeSStolerances function")
  end

  # Tell Sundials the number of state variables, so that I can allocate memory for the Jacobian
  flag = Sundials.CVDense(cvode_mem, neq)
  if flag < 0.0
    error("Error in the CVDense function")
  end

  if settings["jacobian"]
    flag = Sundials.CVDlsSetDenseJacFn(cvode_mem, interface_Jac)
    if flag < 0.0
      error("Error in the CVDlsSetDenseJacFn function")
    end
  end

  # Pass information for the callback
  flag = Sundials.CVodeSetUserData(cvode_mem,{model, parameters, forcings_data, jacobian})
  if flag < 0.0
    error("Error in the CVodeSetUserData function")
  end

  # Set maximum number of steps
  flag = Sundials.CVodeSetMaxNumSteps(cvode_mem, settings["maxsteps"])
  if flag < 0.0
    error("Error in the CVodeSetMaxNumSteps function")
  end

  # Set maximum order of the integration
  flag = Sundials.CVodeSetMaxOrd(cvode_mem, settings["maxord"])
  if flag < 0.0
    error("Error in the CVodeSetMaxOrd function")
  end

  # Set the initial step size
  flag = Sundials.CVodeSetInitStep(cvode_mem, settings["hini"])
  if flag < 0.0
    error("Error in the CVodeSetInitStep function")
  end

  # Set the minimum step size
  flag = Sundials.CVodeSetMinStep(cvode_mem, settings["hmin"])
  if flag < 0.0
    error("Error in the CVodeSetMinStep function")
  end

  # Set the maximum step size
  flag = Sundials.CVodeSetMaxStep(cvode_mem, settings["hmax"])
  if flag < 0.0
    error("Error in the CVodeSetMaxStep function")
  end

  # Set the maximum number of error test fails
  flag = Sundials.CVodeSetMaxErrTestFails(cvode_mem, settings["maxerr"])
  if flag < 0.0
    error("Error in the CVodeSetMaxErrTestFails function")
  end

  # Set the maximum number of nonlinear iterations per step
  flag = Sundials.CVodeSetMaxNonlinIters(cvode_mem, settings["maxnonlin"])
  if flag < 0.0
    error("Error in the CVodeSetMaxNonlinIters function")
  end

  # Set the maximum number of convergence failures
  flag = Sundials.CVodeSetMaxConvFails(cvode_mem, settings["maxconvfail"])
  if flag < 0.0
    error("Error in the CVodeSetMaxConvFails function")
  end

  # Make a first call to the model to check that everything is ok and retrieve the number of observed variables
  forcings = zeros(Float64, length(forcings_data))
  for i in 1:length(forcings)
    forcings[i] = linear_interpolate(forcings_data[i][:,1], forcings_data[i][:,2], times[1])
  end
  first_call  = model(times[1], states, parameters, forcings)

  # Fill up the output matrix with the values for the initial time

  nder = length(first_call[2])
  output = zeros(Float64, (length(times), neq + nder + 1))
  output[1,:] = [times[1], states, first_call[2]]

  # Main time loop. Each timestep call cvode. Handle exceptions and fill up output
  t = times[1]
  y = copy(states)
  for i in 2:length(times)
    # Note that the second t must be an array to build the pointer...
    flag = Sundials.CVode(cvode_mem, times[i], y, [t], Sundials.CV_NORMAL)
       if flag == -1
          error("The solver took mxstep internal steps but could not reach tout.")
       elseif flag == -2
          error("The solver could not satisfy the accuracy demanded by the user for some internal step")
       elseif flag == -3
          error("Error test failures occured too many times during one internal time step or minimum step size was reached")
       elseif flag == -4
          error("Convergence test failures occurred too many times during one internal time step or minimum step size was reached.")
       elseif flag == -5
          error("The linear solver’s initialization function failed.")
       elseif flag == -6
          error("The linear solver’s setup function failed in an unrecoverable manner")
       elseif flag == -7
          error("The linear solver’s solve function failed in an unrecoverable manner")
       elseif flag == -8
          error("The right hand side function failed in an unrecoverable manner")
       elseif flag == -9
          error("The right-hand side function failed at the first call.")
       elseif flag == -10
          error("The right-hand side function had repeated recoverable errors.")
       elseif flag == -11
          error("The right-hand side function had a recoverable errors but no recovery is possible.")
       elseif flag == -25
          error("The time t is outside the last step taken.")
       elseif flag == -26
          error("The output derivative vector is NULL.")
       elseif flag == -27
          error("The output and initial times are too close to each other.")
       end
    output[i,1] = times[i]
    for h in 1:neq
      output[i,h+1] = y[h]
    end
  end
  # If we have observed variables we call the model function again
  if nder > 0
    for i in 1:length(times)
      forcings = [linear_interpolate(forcings_data[j][:,1], forcings_data[j][:,2], times[i]) for j in 1:length(forcings)]
      model_call  = model(times[i], vec(output[i,2:(neq + 1)]), parameters, forcings)
      for h in 1:nder
      @inbounds  output[i,neq + 1 + h] = model_call[2][h]
      end
    end
  end

  return output
end

function simulate(model::OdeModel, times, settings)
    # Check settings (if anything missing assign default value)
    !haskey(settings, "method") && (settings["method"] = "bdf")
    !haskey(settings, "rtol") && (settings["rtol"] = 1e-6)
    !haskey(settings, "atol") && (settings["atol"] = 1e-6)
    !haskey(settings, "maxsteps") && (settings["maxsteps"] = 1e5)
    !haskey(settings, "maxord") && (settings["maxord"] = 5)
    !haskey(settings, "hini") && (settings["hini"] = 0)
    !haskey(settings, "hmin") && (settings["hmin"] = 0)
    !haskey(settings, "hmax") && (settings["hmax"] = 0)
    !haskey(settings, "maxerr") && (settings["maxerr"] = 12)
    !haskey(settings, "maxnonlin") && (settings["maxnonlin"] = 12)
    !haskey(settings, "maxconvfail") && (settings["maxconvfail"] = 12)
    !haskey(settings, "jacobian") && (settings["jacobian"] = false)
    forcs = [[i[2][1] i[2][2]] for i in model.Forcings]
    states = Float64[i[2] for i in model.States]
    parameters = Float64[i[2] for i in model.Parameters]
    simulation = convert(DataFrame, simulate(times, states, parameters,
                            forcs, settings, model.Model, model.Jacobian))
    # Use names for indexing...
    names!(simulation, convert(Array{Symbol, 1}, [symbol("$i") for i in
                           ["time", collect(keys(model.States)), model.Observed]]))
    return simulation
end
