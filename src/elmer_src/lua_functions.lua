function netcdf_cond(t, inter)
  -- floor the float
  year = t - t % 1;

  if year % inter == 0 then
    write = 1;
  else
    write = -1;
  end

  return write;
end

function IfThenElse(condition,t,f)
  if condition then return t else return f end
end


function periodic_surge(omega, time)
  -- @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  --
  -- periodic surge with hard coded periods and beta values
  --
  -- @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  S_p = 5          -- (s)urge      (p)eriod
  Q_p = 45         -- (q)quiescent (p)eriod
  C_p = S_p + Q_p  -- (c)ycle      (p)eriod

  print(time)
  
  -- check if within surge period and that the bed is temperate
  if (time % C_p < S_p) and (omega >= 0.001) then
    beta = 0.001
  else
    beta = 1.000
  end

  return beta
end
