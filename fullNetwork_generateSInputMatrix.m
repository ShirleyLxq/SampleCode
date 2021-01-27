function X_population_SInputMatrix=fullNetwork_generateSInputMatrix(N, numTimeSteps, rate_x, delta_t)
X_population_SInputMatrix=1/delta_t*(rand(N,numTimeSteps)<rate_x*delta_t);
