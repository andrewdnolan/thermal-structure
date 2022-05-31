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

for t = 0,50,0.1
do
  write = netcdf_cond(t, 10)

  if write == 1 then
    print(t)
  end
end
