# List of macros that generate functions to interpolate data that contain hard-coded data
macro generateInputSine(mean, amplitude, freq, offset)
  quote
    function dummy(t)
      $mean + $amplitude*sin($freq*t + $offset)
    end
  end
end

macro generateInputStep(tstep,i0,i1)
  quote
    function dummy(t)
      if (t <= $tstep)
        $i0
      else 
        $i1 
      end
    end
  end
end

macro generateInputSlowStep(t0,t1,i0,i1)
  quote
    function dummy(t)
      if (t <= $t0)
        $i0
      else if ((t > $t0) & (t <= $t1))
        $i0 + (t - $t0)/($t1 - $t0)$i1 
      else 
        $i1
      end
    end
  end
end