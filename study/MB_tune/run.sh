#Execute ElmerSolver
ElmerSolver flux_divergence.sif

# Convert result files into NetCDFs
../../src/elmer2nc/elmer2nc.sh  -r "./mesh/div_flux_lk.result" \
                                -m "./mesh/" \
                                -t 2 \
                                -o "./nc/" 
