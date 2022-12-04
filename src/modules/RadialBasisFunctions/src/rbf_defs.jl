
"Evaluate the Multiquadric (MQ) RBF kernel."
rbf_mq(r) = sqrt(1+r^2)
rbf_mq_diff1(r) = r / sqrt(1+r^2)

"Evaluate the Inverse Multiquadric (IMQ) RBF kernel."
rbf_imq(r) = 1/sqrt(1 + r^2)
rbf_imq_diff1(r) = -r / sqrt(1 + r^2)^3

"Evaluate the Inverse Quadratic (IQ) RBF kernel."
rbf_iq(r) = 1/(1 + r^2)

"Evaluate the Gaussian (GA) RBF kernel."
rbf_ga(r) = exp(-r^2)
rbf_ga_diff1(r) = -2r*exp(-r^2)
rbf_ga_diff2(r) = (4r^2-2)*exp(-r^2)

"Evaluate the polyharmonic spline (PHS) RBF with degree `p`."
rbf_phs(r, p) = r^(2p+1)

"Evaluate the thin plate spline (TPS) RBF with degree `p`."
rbf_tps(r, p) = r^(2p) * log(r)
