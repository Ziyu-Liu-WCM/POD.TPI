POD_TPI_decision <- function(W, MaxDose,
                             dose_vec, event_time_vec, event_vec,
                             n_d, m_d, r_d, p_T, epsilon = c(0.05, 0.05),
                             niter = 1000,
                             type_dose_decision = "mTPI",
                             type_p_prior = "uniform",
                             type_t_model = "pwuniform",
                             type_suspension = 2,
                             type_safety = 2,
                             q1 = 1, q2 = 0.15){

  # stopifnot(all(sapply(list(W, MaxDose, n_d, m_d, r_d), is.numeric)),
  #           W > 0, MaxDose > 0,
  #           length(dose_vec) == length(event_time_vec) & length(dose_vec) == length(event_vec),
  #           all(sapply(dose_vec, function(x) x > 0 & x <= MaxDose)),
  #           all(sapply(event_time_vec, function(x) x > 0 & x <= W)),
  #           all(event_vec %in% c(0,1,NA)),
  #           n_d >= 0, m_d >= 0, r_d >= 0, p_T <= 1 & p_T >= 0,
  #           (p_T - epsilon[1]) >= 0 & (p_T + epsilon[2]) <= 1,
  #           niter >= 10, type_dose_decision %in% c("mTPI", "i3+3"),
  #           type_p_prior %in% c("uniform", "semiparametric"),
  #           type_t_model %in% c("uniform", "pwuniform", "dhazard", "pwhazard"),
  #           type_suspension %in% 0:3,
  #           q1 <= 1 & q1 >= 0.33, q2 <= 0.5 & q2>= 0)

  # Check if W, MaxDose, n_d, m_d, r_d are numeric
  if (!all(sapply(list(W, MaxDose, n_d, m_d, r_d), is.numeric))) {
    stop("Error: W, MaxDose must all be numeric values.")
  }

  # Check if W is greater than 0
  if (W <= 0) {
    stop("Error: W (DLT assessment time window) must be greater than 0.")
  }

  # Check if MaxDose is greater than 0
  if (MaxDose <= 0) {
    stop("Error: MaxDose must be greater than 0.")
  }

  # Check if dose_vec and event_time_vec have the same length and match event_vec length
  if (length(dose_vec) != length(event_time_vec) || length(dose_vec) != length(event_vec)) {
    stop("Error: dose_vec, event_time_vec, and event_vec must have the same length.")
  }

  # Check if all dose_vec values are within the valid range
  if (!all(sapply(dose_vec, function(x) x > 0 & x <= MaxDose))) {
    stop("Error: All values in dose_vec must be greater than 0 and less than or equal to MaxDose.")
  }

  # Check if all event_time_vec values are within the valid range
  if (!all(sapply(event_time_vec, function(x) x > 0 & x <= W))) {
    stop("Error: All values in event_time_vec must be greater than 0 and less than or equal to W.")
  }

  # Check if event_vec contains only 0, 1, or NA
  if (!all(event_vec %in% c(0, 1, NA))) {
    stop("Error: event_vec must contain only 0, 1, or NA.")
  }

  # Check if n_d, m_d, and r_d are non-negative
  if (n_d < 0 || m_d < 0 || r_d < 0) {
    stop("Error: n_d, m_d, and r_d must all be non-negative.")
  }

  # Check if p_T is within the range [0, 1]
  if (p_T < 0 || p_T > 1) {
    stop("Error: p_T must be between 0 and 1.")
  }

  # Check if (p_T - epsilon[1]) >= 0 and (p_T + epsilon[2]) <= 1
  if ((p_T - epsilon[1]) < 0 || (p_T + epsilon[2]) > 1) {
    stop("Error: The interval [p_T - epsilon[1], p_T + epsilon[2]] must be within [0, 1].")
  }

  # Check if niter is at least 10
  if (niter < 10) {
    stop("Error: niter must be at least 10.")
  }

  # Check if type_dose_decision is valid
  if (!type_dose_decision %in% c("mTPI", "i3+3")) {
    stop("Error: type_dose_decision must be either 'mTPI' or 'i3+3'.")
  }

  # Check if type_p_prior is valid
  if (!type_p_prior %in% c("uniform", "semiparametric")) {
    stop("Error: type_p_prior must be either 'uniform' or 'semiparametric'.")
  }

  # Check if type_t_model is valid
  if (!type_t_model %in% c("uniform", "pwuniform", "dhazard", "pwhazard")) {
    stop("Error: type_t_model must be one of 'uniform', 'pwuniform', 'dhazard', or 'pwhazard'.")
  }

  # Check if type_suspension is within the valid range
  if (!type_suspension %in% 0:3) {
    stop("Error: type_suspension must be an integer between 0 and 3.")
  }

  # Check if q1 is within the range [0.33, 1]
  if (q1 < 0.33 || q1 > 1) {
    stop("Error: q1 must be between 0.33 and 1.")
  }

  # Check if q2 is within the range [0, 0.5]
  if (q2 < 0 || q2 > 0.5) {
    stop("Error: q2 must be between 0 and 0.5.")
  }

  current_dose = dose_vec[length(dose_vec)]
  z = dose_vec
  v = event_time_vec
  y = event_vec
  D = MaxDose
  h = c(W/3, 2*W/3, W)
  decision_vec <- c("De-escalation", "Stay", "Escalation")

  hp = list()
  hp$eta = rep(1, length(h))
  if(type_p_prior == "semiparametric"){
    hp$kappa = rep(1/D, D)
    hp$theta = sapply(1:D, getprior, halfwidth = mean(epsilon), target = p_T, nlevel = D)
    hp$c = 1
  }


  n_vec = unname(table(factor(z[y == 1], levels = 1:MaxDose)))
  m_vec = unname(table(factor(z[y == 0], levels = 1:MaxDose)))
  r_vec = unname(table(factor(z[is.na(y)], levels = 1:MaxDose)))

  n_m_vec = n_vec + m_vec
  safety_bound = pbeta(p_T, shape1 = n_vec+1, shape2 = m_vec+1)
  safety_index = which(as.numeric(n_m_vec >= 3 & safety_bound < 0.05) == 1)


  EarlyStop_trigger = (n_m_vec[1] >= 3 & safety_bound[1] < 0.05)

  if((type_safety %in% 1:2) & EarlyStop_trigger){

    if(pbeta(p_T, shape1 = n_vec[1] + 1, shape2 = m_vec[1] + r_vec[1] +1) < 0.05){

      return(sprintf("Safety rule 1 (early termination) triggered. The trial should be terminated."))

    } else {

      return(sprintf("Patient SUSPEND(based on safety rule 1 of early suspend). First dose, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d. No probability calculated",
                     n_vec[1], m_vec[1], r_vec[1]))

    }
  }

  if(type_safety == 2 & length(safety_index) != 0){

    return(sprintf("Safety rule 2 (dose exclusion) triggered. Dose level %d and higher should be suspended from the trial.", min(safety_index)))

  }


  if(r_d == 0){

    # If there is no pending patient, just use dose decision based on complete data

    dose_decision = TITE_dose_decision(n_d, m_d, p_T, epsilon, type_dose_decision)

    if(dose_decision == 1 & current_dose == MaxDose){
      dose_decision = 0
    }

    if(dose_decision == -1 & current_dose == 1){
      dose_decision = 0
    }

    return(list(Decision = sprintf("There is no pending patient. Dose decision is made based on complete data and no probability will be calculated. Dose decision: %s. Dose assignment for this patient: %d. Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d.",
                                   decision_vec[dose_decision+2], current_dose + dose_decision, current_dose, n_d, m_d, r_d)))

  } else {

    ## r_d >= 1. Assign according to POD

    ################################################################################
    ### SUSPENSION RULE

    if(type_suspension == 1 & r_d > (n_d + m_d + r_d) / 2){


      return(list(Decision = sprintf("Patient SUSPEND(based on r_d > (n_d + m_d + r_d) / 2). Current Dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d. No probability calculated",
                                     current_dose, n_d, m_d, r_d)))



    } else if (type_suspension %in% c(0, 2) & n_d + m_d == 0) {

      return(list(Decision = sprintf("Patient SUSPEND(based on n_d + m_d = 0, no patient with observed outcome). Current Dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d. No probability calculated",
                                     current_dose, n_d, m_d, r_d)))

    } else if(type_suspension == 3) {

      dose_decisions_s_d = apply(rbind(n_d + 0:r_d, m_d + r_d:0), 2, function(x){
        TITE_dose_decision(n = x[1], m = x[2], p_T = p_T, epsilon = epsilon, type = type_dose_decision)
      })

      if(current_dose == MaxDose){
        # not possible to escalate
        dose_decisions_s_d[dose_decisions_s_d == 1] = 0
      }

      if(current_dose == 1){
        # not possible to de-escalate
        dose_decisions_s_d[dose_decisions_s_d == -1] = 0
      }

      # safety_all = pbeta(p_T, shape1 = n_d + 0:r_d + 1, shape2 = m_d + r_d:0 +1) < 0.05

      # if pending outcomes do have an impact on decision
      # if(length(unique(dose_decisions_s_d)) > 1 | length(unique(safety_all)) > 1) {
      if(length(unique(dose_decisions_s_d)) > 1) {

        return(list(Decision = sprintf("Patient SUSPEND(based on multiple possible decisions). Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d. No probability calculated",
                                       current_dose, n_d, m_d, r_d)))
      }
    }


    if(type_suspension == 3) {

      # r_d does not affect decision anyway
      dose_decision = TITE_dose_decision(n_d + r_d, m_d + 0, p_T, epsilon, type_dose_decision)

      if(dose_decision == 1 & current_dose == MaxDose){
        dose_decision = 0
      }

      if(dose_decision == -1 & current_dose == 1){
        dose_decision = 0
      }

      prob_dose_decisions <- c(0, 0, 0)
      prob_dose_decisions[dose_decision+2] <- 1

      return(list(Decision = sprintf("Dose decision: %s. Dose assignment for this patient: %d. Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d.",
                                     decision_vec[dose_decision+2], current_dose + dose_decision, current_dose, n_d, m_d, r_d),
                  'Probability of De-escalation' = prob_dose_decisions[1],
                  'Probability of Stay' = prob_dose_decisions[2],
                  'Probability of Escalation' = prob_dose_decisions[3]))

    } else {

      # If type_suspension != 3
      MCMC_spls = TITE_MCMC(z, v, y, D, W, h, p_T, epsilon, hp, niter, p_prior = type_p_prior, t_model = type_t_model)

      v_unobserve_d = v[z == current_dose & is.na(y)]

      p_spls = MCMC_spls$p_spls
      rho_unobserve_spls = MCMC_spls$rho_unobserve_spls

      if(is.null(dim(rho_unobserve_spls))){
        rho_unobserve_spls = matrix(rho_unobserve_spls, 1, length(rho_unobserve_spls))
      }


      # p_d_spls & rho_unobserve_d_spls: r_d * niter matrix
      p_d_spls = p_spls[rep(current_dose, length(v_unobserve_d)), ]

      rho_unobserve_d_spls = rho_unobserve_spls[z[is.na(y)] == current_dose, ]

      prob_DLT_unobserve_d_spls = (p_d_spls - rho_unobserve_d_spls * p_d_spls) / (1 - rho_unobserve_d_spls * p_d_spls)

      # prob_DLT_unobserve_d: r_d length vector
      if(r_d == 1){
        prob_DLT_unobserve_d = mean(prob_DLT_unobserve_d_spls)
      } else{
        prob_DLT_unobserve_d = apply(prob_DLT_unobserve_d_spls, 1, mean)
      }



      # sometimes, dpoisbinom returns negative values if prob is too small
      prob_s_d = dpoisbinom(x = 0:r_d, pp = prob_DLT_unobserve_d)

      dose_decisions_s_d = apply(rbind(n_d + 0:r_d, m_d + r_d:0), 2, function(x){
        TITE_dose_decision(n = x[1], m = x[2], p_T = p_T, epsilon = epsilon, type = type_dose_decision)
      })

      dose_decisions = c(-1, 0, 1)
      prob_dose_decisions = rep(0, 3)
      prob_dose_decisions[1] = sum(prob_s_d[dose_decisions_s_d == -1])
      prob_dose_decisions[2] = sum(prob_s_d[dose_decisions_s_d == 0])
      prob_dose_decisions[3] = sum(prob_s_d[dose_decisions_s_d == 1])

      if(current_dose == MaxDose){
        prob_dose_decisions[2] = prob_dose_decisions[2] + prob_dose_decisions[3]
        prob_dose_decisions[3] = 0
      }

      if(current_dose == 1){
        prob_dose_decisions[2] = prob_dose_decisions[2] + prob_dose_decisions[1]
        prob_dose_decisions[1] = 0
      }

      max_index = which.max(prob_dose_decisions)
      max_prob = max(prob_dose_decisions)


      if(dose_decisions[max_index] == 1){

        if( type_suspension == 2 & max_prob < q1 ){

          # SUSPEND
          return(list(Decision = sprintf("Patient SUSPEND(based on probability thresholds, probability of escalation < q1). Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d.",
                                         current_dose, n_d, m_d, r_d),
                      'Probability of De-escalation' = prob_dose_decisions[1],
                      'Probability of Stay' = prob_dose_decisions[2],
                      'Probability of Escalation' = prob_dose_decisions[3]))


        } else {

          dose_decision = 1

        }

      } else if(dose_decisions[max_index] == 0){

        if(type_suspension == 2 & prob_dose_decisions[1] > q2 & current_dose != 1){

          # SUSPEND
          return(list(Decision = sprintf("Patient SUSPEND(based on probability thresholds, probability of de-escalation > q2). Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d.",
                                         current_dose, n_d, m_d, r_d),
                      'Probability of De-escalation' = prob_dose_decisions[1],
                      'Probability of Stay' = prob_dose_decisions[2],
                      'Probability of Escalation' = prob_dose_decisions[3]))


        } else {

          dose_decision = 0

        }

      } else {

        # if(dose_decisions[max_index] == -1)

        dose_decision = -1
      }

      if(dose_decision == 1 & current_dose == MaxDose){
        dose_decision = 0
      }

      if(dose_decision == -1 & current_dose == 1){
        dose_decision = 0
      }

      return(list(Decision = sprintf("Dose decision: %s. Dose assignment for this patient: %d. Current dose: %d, DLT Patients: %d, Non-DLT Patients: %d, Pending Patients: %d.",
                                     decision_vec[dose_decision+2], current_dose + dose_decision, current_dose, n_d, m_d, r_d),
                  'Probability of De-escalation' = prob_dose_decisions[1],
                  'Probability of Stay' = prob_dose_decisions[2],
                  'Probability of Escalation' = prob_dose_decisions[3]))

    }
  }
}


TITE_dose_decision = function(n, m, p_T, epsilon, type = "i3+3"){

  # n: number of DLTs
  # m: number of non-DLTs

  if(type == "mTPI"){

    interval_length = sum(epsilon)
    K_1 = ceiling((p_T - epsilon[1]) / interval_length)
    K_2 = ceiling((1 - p_T - epsilon[2]) / interval_length)
    intervals = seq(p_T - epsilon[1] - K_1 * interval_length,
                    p_T + epsilon[2] + K_2 * interval_length, length.out = K_1 + K_2 + 2)
    intervals[1] = 0
    intervals[K_1 + K_2 + 2] = 1
    interval_lower = intervals[1:(K_1 + K_2 + 1)]
    interval_upper = intervals[2:(K_1 + K_2 + 2)]

    weights = pbeta(interval_upper, shape1 = n+1, shape2 = m+1) -
      pbeta(interval_lower, shape1 = n+1, shape2 = m+1)

    chosen_model = which.max(weights)
    if(chosen_model <= K_1){
      return(1)
    } else if(chosen_model ==(K_1 + 1)){
      return(0)
    } else{
      return(-1)
    }

  } else if(type == "i3+3"){

    if( n / (n+m) <  p_T - epsilon[1] ){
      return(1)
    } else if(n / (n+m) <= p_T + epsilon[2]){
      return(0)
    } else{
      if( (n-1) / (n+m) < p_T - epsilon[1] ){
        return(0)
      } else{
        return(-1)
      }
    }

  }

}



TITE_MCMC = function(z, v, y, D, W, h, p_T, epsilon, hp, niter,
                     burnin = NULL, p_prior = "uniform", t_model = "pwuniform"){

  # z: n dimensional vector (n is the number of patients), dose levels
  # v: n dimensional vector, min of DLT time and censoring time (that we can observe)
  # y: n dimensional vector, binary DLT indicator. NA if DLT not assessed (follow-up time < W and DLT not occured yet)
  # D: number of dose levels
  # W: DLT assessment window
  # h: K dimensional vector (h_1, ..., h_K) that partitions the (0, W] interval. h_K = W, but h_0 = 0 not included
  # p_T: the target DLT rate
  # epsilon: 2 dimensional vector, the bounds of the equivalence interval [p_T - epsilon_1, p_T + epsilon_2]
  # hp: a list of hyperparameters
  # hp$theta: a D * D matrix, each row corresponds to a dose, and each column corresponds to a (unknown) true MTD
  # niter: number of MCMC iterations
  # burnin: number of MCMC iterations treated as burnin

  # p_prior: prior for p. Choices: uniform or semiparametric
  # t_model: model for time-to-toxicity. Choices: uniform, pwuniform (piecewise-uniform),
  #          dhazard (discrete hazard), pwhazard (piecewise constant hazard)

  # RETURN: (niter - burnin) samples

  if(is.null(burnin)){
    burnin = floor(niter/3)
  }

  # data counts
  z_observe = z[!is.na(y)]
  v_observe = v[!is.na(y)]
  y_observe = y[!is.na(y)]
  n_vec = unname(table(factor(z_observe[y_observe == 1], levels = 1:D)))
  m_vec = unname(table(factor(z_observe[y_observe == 0], levels = 1:D)))

  # n_dot_vec: count number of observed DLTs in each of the K sub-interval
  n_dot_vec = unname(table(cut(v_observe[y_observe == 1], breaks = c(0, h))))

  z_unobserve = z[is.na(y)]
  v_unobserve = v[is.na(y)]

  r = length(z_unobserve)

  #cat(sprintf("%d\n", length(v_unobserve)))
  eta = hp$eta

  if(p_prior == "semiparametric"){
    kappa = hp$kappa
    theta = hp$theta  # D * D matrix
    c = hp$c
  }


  K = length(h)

  if(p_prior == "uniform"){

    p_pro = matrix(rbeta(n = D*niter, shape1 = rep(n_vec + 1, niter),
                         shape2 = rep(m_vec + 1, niter)), D, niter)

    if(t_model == "uniform"){

      # no omega
      # rho_unobserve_pro: r * niter

      rho_unobserve_pro = matrix(rep(v_unobserve / W, niter), r, niter)


    } else if(t_model == "pwuniform"){

      omega_pro = t(rdirichlet(niter, n_dot_vec + eta))

      # sapply(v_unobserve, beta_fun, h = h) is a K * r matrix
      # K: length of h; r: length of v_unobserve

      rho_unobserve_pro = crossprod(sapply(v_unobserve, beta_fun, h = h), omega_pro)


    } else if(t_model == "dhazard"){

      # v_toxic: n length vector. time-to-toxicity for the patients with observed DLTs
      v_toxic = v_observe[y_observe == 1]

      if(length(v_toxic) == 0){

        # then pending patients receive no weight, and just use observed data for inference

        rho_unobserve_pro = matrix(0, r, niter)

      } else {

        h_new = sort(unique(c(v_toxic, W)))

        # r length vector, for each i = 1, ..., r, finding the largest k such that h[k] <= v[i]
        k_lgst_unobserve = sapply(v_unobserve, function(y) { Position(function(x) x <= y, h_new, right = TRUE) } )

        if(all(is.na(k_lgst_unobserve))){
          # this means, if the follow-up times for all pending patients are shorter than
          # the observed DLTs, then all pending patients won't receive weights

          rho_unobserve_pro = matrix(0, r, niter)


        } else {

          K_new = length(h_new) - 1
          # h_freq: K_new + 1 length, calculates number of v[i] = h[k]
          h_freq = unname(table(factor(v_toxic, levels = h_new)))

          omega_post_shape1 = 1/2 + h_freq[1:K_new]
          omega_post_shape2 = 1/2 + rev(cumsum(rev(h_freq)) - rev(h_freq))[1:K_new]

          # omega_pro: K_new * niter
          omega_pro = matrix(rbeta(n = K_new * niter, shape1 = rep(omega_post_shape1, niter),
                                   shape2 = rep(omega_post_shape2, niter)), K_new, niter)

          rho_pro = 1 - apply(1 - omega_pro, 2, cumprod)
          rho_pro = matrix(rho_pro, K_new, niter)

          # rho_unobserve_pro: r * niter
          rho_unobserve_pro = rho_pro[k_lgst_unobserve, ]
          rho_unobserve_pro[is.na(rho_unobserve_pro)] = 0

        }


      }



    } else if(t_model == "pwhazard"){

      # v_toxic: n length vector. time-to-toxicity for the patients with observed DLTs
      v_toxic = v_observe[y_observe == 1]

      # h_length: K length vector, h[k] - h[k-1]
      h_length = h - c(0, h[-K])

      if(length(v_toxic) == 0){

        # omega_post_shape: K length vector
        omega_post_shape = K / (2 * (W * (K - 1:K + 0.5)))
        omega_post_rate = rep(1/2, K)

      } else {

        # omega_post_shape: K length vector
        omega_post_shape = K / (2 * (W * (K - 1:K + 0.5))) + n_dot_vec
        # beta_fun_toxic: K * n matrix
        beta_fun_toxic = sapply(v_toxic, beta_fun, h = h)

        omega_post_rate = 1/2 + apply(beta_fun_toxic, 1, sum) * h_length

      }

      # omega_pro: K * niter matrix
      omega_pro = matrix(rgamma(n = K * niter, shape = rep(omega_post_shape, niter),
                                rate = rep(omega_post_rate, niter)), K, niter)

      # beta_fun_unobserve: K * r matrix
      beta_fun_unobserve = sapply(v_unobserve, beta_fun, h = h)

      rho_unobserve_pro = (1 - exp( -t(beta_fun_unobserve) %*% (omega_pro * h_length)))



    }

    rho_unobserve_pro = matrix(rho_unobserve_pro, r, niter)
    log_accept_probs = apply(log(1 - p_pro[z_unobserve, ] * rho_unobserve_pro), 2, sum)

    log_accept_probs = log_accept_probs - c(0, log_accept_probs[-niter])
    log_accept_probs[1] = 0

    accept_spls = as.numeric(log(runif(niter)) < log_accept_probs)

    # cat(sprintf("Acceptance rate for the independent M-H: %.2f%%.\n", mean(accept_spls)*100))
    index_spls = 1:niter
    index_spls[accept_spls == 0] = NA
    index_spls = na.locf(index_spls)

    index_spls = index_spls[(burnin+1):niter]

    # return samples of p and samples of rho_unobserve (a function of omega)
    list_spls = list(p_spls = p_pro[ , index_spls], rho_unobserve_spls = rho_unobserve_pro[ , index_spls])

    return(list_spls)

  } else if(p_prior == "semiparametric"){

    # D * D matrix
    p_lower_bound = matrix(0, D, D)
    p_lower_bound[lower.tri(p_lower_bound)] = p_T + epsilon[2]
    diag(p_lower_bound) = p_T - epsilon[1]

    p_upper_bound = matrix(1, D, D)
    p_upper_bound[upper.tri(p_upper_bound)] = p_T - epsilon[1]
    diag(p_upper_bound) = p_T + epsilon[2]

    p_prior_shape1 = c * theta + 1
    p_prior_shape2 = c * (1 - theta) + 1

    p_post_shape1 = sweep(p_prior_shape1, 1, n_vec, "+") # add n_vec by row
    p_post_shape2 = sweep(p_prior_shape2, 1, m_vec, "+") # add m_vec by row



    gamma_pro = sample(x = 1:D, size = niter, replace = TRUE, prob = kappa)

    # D * niter matrix
    p_lower_bound = p_lower_bound[ , gamma_pro]
    p_upper_bound = p_upper_bound[ , gamma_pro]
    p_prior_shape1 = p_prior_shape1[ , gamma_pro]
    p_prior_shape2 = p_prior_shape2[ , gamma_pro]
    p_post_shape1 = p_post_shape1[ , gamma_pro]
    p_post_shape2 = p_post_shape2[ , gamma_pro]

    # (D * niter) * 4 matrix
    p_post_bound_shape = matrix(c(p_lower_bound, p_upper_bound, p_post_shape1, p_post_shape2), D * niter, 4)

    p_pro = apply(p_post_bound_shape, 1, function(x){
      rtrunc(n = 1, spec = "beta", a = x[1], b = x[2], shape1 = x[3], shape2 = x[4])
    } )

    p_pro = matrix(p_pro, D, niter)

    omega_pro = t(rdirichlet(niter, n_dot_vec + eta))

    rho_unobserve_pro = crossprod(sapply(v_unobserve, beta_fun, h = h), omega_pro)

    # r * niter matrix, r = r_1 + ... + r_D, i.e. number of all pending patients
    log_accept_probs1 = log(1 - p_pro[z_unobserve, ] * rho_unobserve_pro)

    log_accept_probs1 = apply(log_accept_probs1, 2, sum)

    # D * niter matrix
    log_accept_probs2 = - lbeta(p_prior_shape1, p_prior_shape2) -
      log(pbeta(p_upper_bound, shape1 = p_prior_shape1, shape2 = p_prior_shape2) -
            pbeta(p_lower_bound, shape1 = p_prior_shape1, shape2 = p_prior_shape2)) +
      lbeta(p_post_shape1, p_post_shape2) +
      log(pbeta(p_upper_bound, shape1 = p_post_shape1, shape2 = p_post_shape2) -
            pbeta(p_lower_bound, shape1 = p_post_shape1, shape2 = p_post_shape2))

    log_accept_probs2 = apply(log_accept_probs2, 2, sum)

    log_accept_probs = log_accept_probs1 + log_accept_probs2
    log_accept_probs = log_accept_probs - c(0, log_accept_probs[-niter])
    log_accept_probs[1] = 0

    accept_spls = as.numeric(log(runif(niter)) < log_accept_probs)

    index_spls = 1:niter
    index_spls[accept_spls == 0] = NA
    index_spls = na.locf(index_spls)

    index_spls = index_spls[(burnin+1):niter]

    return(list(gamma_spls = gamma_pro[index_spls], p_spls = p_pro[ , index_spls],
                rho_unobserve_spls = rho_unobserve_pro[ , index_spls]))

  }


}
