# Cellular Automata Simulation of SIR Model

This project is about implementing cellular automata simulations on a SIR Model, which considers three major population groups in the presence of a disease: 
+ susceptibles (S) that have no immunity from the disease
+ infecteds (I) that have the disease and can spread it to others, and 
+ recovereds (R) that have recovered from the disease and are immune to further infection. 

# Model 1: 
Suppose an individual is at each grid point. An individual can be well and susceptible (value = 0) to a disease, sick with the disease that has two phases (values 1 to 2), or immune (values 3 to 7). The infection lasts exactly 2 days, and immunity lasts exactly 5 days before the individual becomes susceptible again.

+	susceptible: value = 0
+	infectious: value = 1 or 2, indicating the day of infection
+	immune: value = 3 to 7, where the day of immunity is the cell value minus 2. For example, on day 1 of immunity, the cell value is 3

probSusceptible: the probability the individual is initially susceptible
probInfectious: the probability that an individual that is not susceptible is infectious initially
probImmune: the probability the individual is initially immune

Uniformly distribute the infected individuals between day 1 and day 2 of the infection. In the initialization, uniformly distribute immune individuals with values 3 through 7. 
The following rules apply: 
+	If an individual is susceptible and a neighbor is infected, the individual becomes infected;
+	The infection lasts for 2 days exactly;
+	Immunity lasts for 5 days exactly, after which time the individual again becomes susceptible 100%

Color the graphic as follows:
+	Susceptible: full green;
+	Infectious: red; full red on the first day, paler red on the second;
+	Immune: blue; full red on the first day, paler blue thereafter. 

The simulation does not occur around the edges of the grid – the edges are designed to be walls that do not interact with populations in any of the stages . 

# Model 2: 
Model 2 is a modified model based on Exercise 1. There’re additional assumptions:
+	Immune period longer than 5 days is allowable. At the end of the immune period, there is the probBeSusceptible, the probability that an individual who has been immune for 5 days will become susceptible. Thus, someone who is exposed to the virus might not become sick, and a person might have longer immunity than 5 days.
+	The probability that a susceptible individual will get sick is the percentage of sick neighbors. It’s reasonable to assume that presence of more infected people would result in a higher chance of being infected. 

# Model 3: 
Model 3 is a modified model based on Model 2. There’re additional assumptions:
+	Periodic variation of a parameter in a model is presented. Implement periodic seasonal infection probabilities in the exercise. In particular, the period of the flu season is 60 days (present for 60 days and absent for 60 days). For off-season span in the cycle, the infection rate is 0; for the 60 days flu season, the infection rate starts at 0, rises to its maximal value at 30 days, and then fall back to 0 at the end of the 60-day season. 
+	This infection rate will result in a longer infection period than usual (at the end of the infection period, the infection rate decides whether you would become immune or you stay being infected for another day. 

# Model 4: 
Model 4 is a modified model based on Model 1. There’re additional assumptions:
+	Implementation of a vaccine with certain threshold for total infection rate of the population. Assume the government would put a vaccine with certain success rate into practice when the infection rate (total infection/total population). 
+	The vaccine would be taken away or ignored when the total infection rate is low enought
 
## Reference consulted: 
Angela B. Shiflet and George W. Shiflet. 2014. Introduction to Computational Science: Modeling and Simulation for the Sciences (2nd. ed.). Princeton University Press, USA.
