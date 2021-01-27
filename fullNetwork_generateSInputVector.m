function X_population_SInputVector=fullNetwork_generateSInputVector(N, rate_x, delta_t)
X_population_SInputVector=1/delta_t*(rand(N,1)<rate_x*delta_t);