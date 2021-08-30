# import the necessarily libraries
library(ggplot2)

#################################################################################################
#### Gillespie algorithm ####
#################################################################################################

# define the function to output the results from gillespie algorithm
gillespie <- function(S, I, R, t, beta, r, endsim, stoptime) {
  # set names for the transition states
  states <- c("infected", "recovered")
  # create vectors to store all states
  Susceptible <- S
  Infected <- I
  Recovered <- R
  transition_times <- t
  counter <- 1
  
  if (endsim == "yes") {
    while ((beta*S*I + r*I) > 0 && (counter < stoptime)) {
      # increase the counter
      counter <- counter + 1
      # calculate the time spent in this transition state
      time_step <- rexp(1, rate=(beta*S*I + r*I))
      # add it to the current time
      transition_times <- c(transition_times, transition_times[length(transition_times)] + time_step)
      # calculates probabilities
      p <- beta*S*I/(beta*S*I + r*I)
      q <- 1-p
      # sample a transition
      transition <- sample(states, 1, prob = c(p,q))
      
      if (transition == "infected") {
        S <- S - 1
        I <- I + 1
        Susceptible <- c(Susceptible, S)
        Infected <- c(Infected, I)
        Recovered <- c(Recovered, R)
      } else {
        I <- I - 1
        if (I == 0) {
          # restore back the original values and start again
          S <- Susceptible[1]
          Susceptible <- S
          I <- Infected[1]
          Infected <- I
          R <- Recovered[1]
          Recovered <- R
          transition_times <- t
          counter <- 1
          next
        }
        R <- R + 1
        Susceptible <- c(Susceptible, S)
        Infected <- c(Infected, I)
        Recovered <- c(Recovered, R)
      }
    }
  } else {
    while ((beta*S*I + r*I) > 0) {
      # calculate the time spent in this transition state
      time_step <- rexp(1, rate=(beta*S*I + r*I))
      # add it to the current time
      transition_times <- c(transition_times, transition_times[length(transition_times)] + time_step)
      # calculates probabilities
      p <- beta*S*I/(beta*S*I + r*I)
      q <- 1-p
      # sample a transition
      transition <- sample(states, 1, prob = c(p,q))
      
      if (transition == "infected") {
        S <- S - 1
        I <- I + 1
        Susceptible <- c(Susceptible, S)
        Infected <- c(Infected, I)
        Recovered <- c(Recovered, R)
      } else {
        I <- I - 1
        R <- R + 1
        Susceptible <- c(Susceptible, S)
        Infected <- c(Infected, I)
        Recovered <- c(Recovered, R)
      }
    }
  }
  
  # return a list of values
  return(list("Susceptible" = Susceptible, "Infected" = Infected, "Recovered" = Recovered, "time" = transition_times))
}

# run the algorithm
# set the random seed for reproducible results
set.seed(1234)
gil_output_fixed <- gillespie(S = 999, I = 1, R = 0, t = 0, beta = 0.05, r = 10, endsim = "no", stoptime = Inf)

# store the results in a data frame for plotting them with ggplot
gil_results <- data.frame(T = gil_output_fixed$time, S = gil_output_fixed$Susceptible, I = gil_output_fixed$Infected, R = gil_output_fixed$Recovered)

# plot the results
ggplot(data = gil_results, aes(x=T)) +
  geom_line(aes(y=S, colour = "Susceptible")) +
  geom_line(aes(y=I, colour = "Infectious")) +
  geom_line(aes(y=R, colour = "Recovered")) +
  scale_colour_manual("",
                      breaks = c("Susceptible", "Infectious", "Recovered"),
                      values = c("blue", "orange", "green")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Number of people")

#################################################################################################
#### Noisy observations ####
#################################################################################################

# run the gillespie algorithm
set.seed(1234)
gil_output_fixed <- gillespie(S = 999, I = 1, R = 0, t = 0, beta = 0.05, r = 10, endsim = "no", stoptime = Inf)

# create a function to generate noise
noise_gen <- function(gen_data, steps) {
  # extract the final time
  t_final <- gen_data$time[length(gen_data$time)]
  # calculate the time step
  time_step <- t_final/steps
  # extract the closest time indices to time observations
  time_line <- time_step
  for (i in 1:(steps-20)) {
    time_line <- c(time_line, time_line[i]+time_step)
  }
  indices <- NULL
  for (i in 1:length(time_line)) {
    print(i)
    for (j in 1:length(gen_data$time)) {
      if (gen_data$time[j] < time_line[i] && time_line[i] < gen_data$time[j+1]) {
        indices <- c(indices, j)
        break
      }
    }
  }
  # extract the number of true infected individuals at observations times
  infected <- gen_data$Infected[indices]
  
  # extract the exact times of the observations
  inf_times <- gen_data$time[indices]
  
  # generate noisy observations
  noisy_obs <- NULL
  for (i in 1:length(infected)) {
    noisy_obs <- c(noisy_obs, rpois(1, lambda = infected[i]))
  }
  
  # return the noisy observations and the time axis
  return(list("noisy" = noisy_obs, "time" = inf_times, "steps" = steps, "step" = time_step, "last_index" = indices[length(indices)]))
}

# generate noisy data
noise <- noise_gen(gil_output_fixed, 50)

# store the noisy data into a data frame
noise_frame <- data.frame(t = noise$time, n = noise$noisy)

# remove the times after the last noisy observation
noise_results <- data.frame(T = gil_output_fixed$time[1:noise$last_index], I = gil_output_fixed$Infected[1:noise$last_index])

# plot the true data with the generated noise points in red
ggplot(data = noise_results, aes(x=T)) +
  geom_line(aes(y=I), col = "orange") +
  geom_point(data=noise_frame, aes(x=t, y=n), col = "red", shape = 4, size = 2) +
  labs(x = "Time", y = "Number of people")