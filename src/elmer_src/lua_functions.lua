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
