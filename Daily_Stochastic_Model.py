    #!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)

import random as r
import numpy as np
import scipy
from scipy.stats import chi2_contingency
import scipy.stats as stats
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import time
import pandas as pd
from scipy.optimize import curve_fit
plt.switch_backend('agg')
print('Script started', flush=True)
########################## model global parameters ############################

# To test the effect of a parameter, comment it and modify the parameter later : 
# Modify the parameter you want to test in the if (line 418) and else (line 428), as well as in the Output (line 510).

beta = 3                      # encounter rate between active males and active virgin females
d = 0.1                       # death rate applied to whole population
dlt_e = 0.1                   # energy linked death rate of active individuals
dlt_c = 0.1                   # density dependant competition death rate based on male-male competition
dlt_l = 0.1                   # egg laying cost for females
offspring_number = 5          # mean offspring number per female per day 
mu = 1/25                     # sexual activity window position standard deviation
mating_threshold = 1          # maximum number of mating events a female can go through


lower = 0                     # used to bound drawing in truncated normal distributions
upper = 1                     # used to bound drawing in truncated normal distribution

K = 1000                    # carrying capacity regarding competition between females for access to host plants to lay eggs

e1 = 0.5                      # peak emergence time 

ve1 = 0.05                    # emergence time standard deviation

mu_g = 0.1                    # mutation rate for genes
genome_size = 10              # number of genes of each individual
# incompatibility_lim = 1      # maximum genetic distance allowing breeding


######################### individuals attributes ##############################

class individuals:                           # Individuals in our population belongs to this class
    """
    Class to define the attributes of individals in the model

    """
    def __init__(self, sex, ha, ha_male, wa, wa_male, f, sp, gene, gene_male, nb_offsprings,emerg,birthday,generation): # Each individual is defined by the following attributes :
        self.sex = sex                       # Male or Female
        self.ha = ha                         # center of sexual activity timing
        self.ha_male = []                    # sexual activity timing storage for mating events
        self.wa = 0.15                       # sexual activity window  (fixed)
        self.wa_male = []                    # sexual activity window width storage for mating events
        self.f = f                           # fertilisation state 0 (virgin), 1 (mated)
        self.sp = sp                         # deprecated
        self.gene = gene                     # np.array([a,b,c,d]) where a,b,c,d are between 0 and 1
        self.gene_male = []                  # for females to save the male genes when mated        
        self.nb_offsprings = nb_offsprings   # number of offsprings birthed for a female
        self.emerg = emerg                   # emergence time of the individual
        self.birthday = birthday
        self.generation = generation
        
####################### mutations related functions ###########################

def g_transfer(femelle, male): 

    gene = [0]* len(femelle)        # we intialize a dummy gene to fill
    for i in range(len(femelle)):   
        if r.random() > 0.5:        # with a 50/50 probability on each gene
            gene[i] = femelle[i]    # either an allele from the mother
        else:
            gene[i] = male[i]       # or the father
    
    for i in range(len(gene)):      # we then give it a random change to mutate
        if r.random() < mu_g:
            gene[i] = gene[i]^1     # by flipping the value
            
    return gene
        
########### functions to obtain the population categories at time t ###########


def list_virgin_patrolling_females(pop,ti): # get a list of active virgin females in the population
    fem_v = []
    for i in range(len(pop)): 
        if pop[i].sex == "F" and pop[i].f == 0 and ((pop[i].ha - pop[i].wa) <= ti <= (pop[i].ha + pop[i].wa)):
            fem_v.append(pop[i])
    return(fem_v)

def list_virgin_nonpatrolling_females(pop,ti): # get a list of non-active virgin females in the population
    fem_v = []
    for i in range(len(pop)): 
        if pop[i].sex == "F" and pop[i].f == 0 and ((pop[i].ha - pop[i].wa) > ti 
                                    or ti > (pop[i].ha + pop[i].wa)):
            fem_v.append(pop[i])
    return(fem_v)

def list_mated_females_1(pop): # get a list of mated females within species 1 in the population
    fem_f1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "F" and pop[i].f == 1 and pop[i].sp == 1:
            fem_f1.append(pop[i])
    return(fem_f1)


def list_patrolling_males_1(pop,ti): # get a list of active males in sp 1 at time "ti" in the population
    mal_a1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "M" and ((pop[i].ha - pop[i].wa) <= ti <= (pop[i].ha + pop[i].wa)):
            mal_a1.append(pop[i])
    return(mal_a1)


def list_nonpatrolling_males(pop,ti): # get a list of non-active males at time "ti" in the population
    mal_i = []
    for i in range(len(pop)): 
        if pop[i].sex == "M" and ((pop[i].ha - pop[i].wa) > ti 
                                    or ti > (pop[i].ha + pop[i].wa)):
            mal_i.append(pop[i])
    return(mal_i)

def list_eggs_weight(pop):  # Compute the egg weight when selecting a female to die
    egg_weight= []
    for i in range(len(pop)): 
            egg_weight.append(pop[i].nb_offsprings)
    return(egg_weight)

def logistic(diff,incompatibility_lim,k=50):    #Compute the number of offsprings based on genetic difference
    return 1 - 1 / (1 + np.exp(-k * (diff - incompatibility_lim)))


########################### time period cutting functions #####################

def patrolling_schedule(pop):          # get a dictionary of hours of transition from non-active to active and vice versa
    ht = {}                                  # dictionnaries allows us to link events with their associated time
    for i in range(len(pop)): 
        
        if pop[i].emerg > pop[i].ha + pop[i].wa:   #sexual activity is before emergence, skip individual
            pop[i].emerg=0       
            continue
        if pop[i].ha - pop[i].wa < pop[i].emerg < pop[i].ha + pop[i].wa:   #emergence is during sexual activity
            if pop[i].sex == "M":
                ht[pop[i].emerg] = "entering male " + str(i)
                pop[i].emerg=0   
            if pop[i].ha + pop[i].wa > 1 :   # avoids to exit [0,1]
                ht[1] = "exiting male " + str(i)
                pop[i].emerg=0   
            else:
                ht[pop[i].ha + pop[i].wa] = "exiting male " + str(i)    
                pop[i].emerg=0  
                
            if pop[i].sex == "F":
                ht[pop[i].emerg] = "entering female " + str(i)
                pop[i].emerg=0   
            if pop[i].ha + pop[i].wa > 1 :   # avoids to exit [0,1]
                ht[1] = "exiting female " + str(i)
                pop[i].emerg=0   
            else:
                ht[pop[i].ha + pop[i].wa] = "exiting female " + str(i)  
                pop[i].emerg=0   
                
                
        else : 
            if pop[i].sex == "M":
                if pop[i].ha - pop[i].wa < 0 :   # avoids to exit [0,1]
                    ht[0] = "entering male " + str(i)
                else:
                    ht[pop[i].ha - pop[i].wa] = "entering male " + str(i) 
                    
                if pop[i].ha + pop[i].wa > 1 :   # avoids to exit [0,1]
                    ht[1] = "exiting male " + str(i)
                else:
                    ht[pop[i].ha + pop[i].wa] = "exiting male " + str(i)  
                    
            if pop[i].sex == "F":
                if pop[i].ha - pop[i].wa < 0 :   # avoids to exit [0,1]
                    ht[0] = "entering female " + str(i)
                else:
                    ht[pop[i].ha - pop[i].wa] = "entering male " + str(i) 
                    
                if pop[i].ha + pop[i].wa > 1 :   # avoids to exit [0,1]
                    ht[1] = "exiting female " + str(i)
                else:
                    ht[pop[i].ha + pop[i].wa] = "exiting female " + str(i)  
    return(ht,pop)


    
def emergences_and_schedule(pop,incompatibility_lim,ve1,e1,K,mu,current_day):   
    Mated_females = [individu for individu in pop if individu.sex == "F" and len(individu.gene_male)>0]         
    b = offspring_number * (1-len(Mated_females)/K) # number of offspring per female per day is density dependant depending on mated females 
    generation_time = []
    if b < 0 :
        b = 0
    emergence_list = []                                         # list of all emergences of the day
    emergence_schedule = {}                                     # dict of emergence schedule of offsprings
    for i in range(len(Mated_females)): 
        r_ind = r.randint(0, len(Mated_females[i].gene_male)-1)  #Choosing father
        dist = scipy.spatial.distance.hamming(Mated_females[i].gene, Mated_females[i].gene_male[r_ind], w=None)
        
        if dist <= incompatibility_lim:

            offsp = np.random.poisson(b)   #Number of offsprings
            genes = g_transfer(Mated_females[i].gene, Mated_females[i].gene_male[r_ind])  #Generating offspring genotype
        
            
            for j in range(offsp):  
                ha = stats.truncnorm.rvs((lower - (Mated_females[i].ha+Mated_females[i].ha_male[r_ind])/2)/mu,(upper-(Mated_females[i].ha+Mated_females[i].ha_male[r_ind])/2)/mu, (Mated_females[i].ha+Mated_females[i].ha_male[r_ind])/2, mu)
                emergence_list.append(individuals(r.sample(["M","F"],1)[0],         # creation of new emerging individuals
                                              ha, 
                                              None,
                                              0.15,
                                              None,
                                              0,                              # emerging individual are virgins
                                              Mated_females[i].sp,            # deprecated
                                              genes,
                                              [],
                                              0,
                                              stats.truncnorm.rvs((lower - e1)/ve1,(upper-e1)/ve1, e1, ve1),
                                              current_day,
                                              None))
                
                if Mated_females[i].nb_offsprings == 0:
                    generation_time.append((current_day - Mated_females[i].birthday) + (emergence_list[-1].emerg - Mated_females[i].emerg))
                              
                Mated_females[i].nb_offsprings += 1
        
            
    for k in range(len(emergence_list)):                        # drawing of emergence time for each new individual
        emergence_schedule[emergence_list[k].emerg] = "emergence individu " + str(k) + emergence_list[k].sex
    
    return(emergence_list,emergence_schedule,generation_time)                   # returns emergence list AND emergence with associated time

########################### event drawing #####################################


def event(pop,ti, T, dlt_c, mating_threshold,current_day):                          # determine if an event happens between ti and T, and chooses which one
    nb_patrolling_males_1 = len(list_patrolling_males_1(pop,ti))        # calculating each population category abundance
    nb_nonpatrolling_males = len(list_nonpatrolling_males(pop,ti))
    nb_virgin_patrolling_females = len(list_virgin_patrolling_females(pop,ti))
    nb_virgin_nonpatrolling_females = len(list_virgin_nonpatrolling_females(pop,ti))
    nb_mated_females_1 = len(list_mated_females_1(pop))
    Etot = np.sum(list_eggs_weight(pop))

    Ei = []
    lifespan = None
    
    lr = beta * (nb_patrolling_males_1) * nb_virgin_patrolling_females         # calculating lambda_r and lambda_d used in exponential law
    ld = d * (nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1 + nb_nonpatrolling_males) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot

    if lr + ld == 0:                                            # specific case if population go extinct
        print("population extinct")                             # return t = 100 to indicate error
        return(lifespan,100)
    draw = np.random.exponential(1/(lr + ld))                 # drawing in exponential law of parameter lambda = lambda_r + lambda_d
    if draw < T - ti:
        uni = np.random.uniform(0,1)                            # drawing in unif(0,1) to set which event is happening
        
        if uni < lr/(lr+ld):                                  # event : meeting between active male and active virgin female
            father = r.sample(list_patrolling_males_1(pop,ti), 1)[0]  # father is randomly sampled among active males
            mother = r.sample(list_virgin_patrolling_females(pop,ti),1)[0]    # mother is randomly sampled among active virgin females
            
            mother.ha_male.append(father.ha)                        # mother stores father's patrolling window position
            mother.wa_male.append(father.wa)                        # mother stores father's patrolling window width
            mother.gene_male.append(father.gene)                    # mother stores father's genome
            
            if len(mother.gene_male) == mating_threshold:     # mother becomes fertilized if she reaches the max number of mating events
                mother.f = 1 
            if mother.nb_offsprings == 0:
                mother.generation = [current_day,ti]
                
        elif uni < (lr/(lr+ld)) + (d * nb_nonpatrolling_males)/(lr + ld):       # event : death of non-active male
            dead = r.sample(list_nonpatrolling_males(pop,ti),1)[0]              # the dead is randomly chosen among non-active males at time ti
            lifespan = current_day - dead.birthday
            pop.remove(dead) 
            
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_patrolling_females) + dlt_e * nb_virgin_patrolling_females)/(lr+ld):    # event : death of virgin active female  
            dead = r.sample(list_virgin_patrolling_females(pop,ti),1)[0]    
            lifespan = current_day - dead.birthday
            pop.remove(dead) 
        
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females) + dlt_e * nb_virgin_patrolling_females)/(lr+ld): # event : death of virgin non-active female    
            dead = r.sample(list_virgin_nonpatrolling_females(pop,ti),1)[0]                                      
            lifespan = current_day - dead.birthday
            pop.remove(dead) 
        
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1) + dlt_e * nb_virgin_patrolling_females + dlt_l * Etot)/(lr+ld):
            Esum = 0
            for l in range(nb_mated_females_1):
                Ei = list_mated_females_1(pop)[l]
                Esum += Ei.nb_offsprings 
        
                if uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females) + d*(l+1) + dlt_e * nb_virgin_patrolling_females  + dlt_l*Esum)/(lr+ld): # the dead is chosen among mated females
                    lifespan = current_day - Ei.birthday
                    pop.remove(Ei)
                    break

        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot)/(lr+ld):
                dead = r.sample(list_patrolling_males_1(pop,ti),1)[0]
                lifespan = current_day - dead.birthday
                pop.remove(dead)  
                    
        return(lifespan,ti + draw)
    else:
        return(lifespan,T) # no event happens when the draw in exponential law falls after T - ti interval


######################## simulation study #####################################

def list_ha(pop): # get a list of the center of sexual activity of the population
    mal1 = []
    for i in range(len(pop)): 
            mal1.append(abs(pop[i].ha))
    return(mal1)


def list_ha_and_g_of_males(pop): # get a list of the center of sexual activity of males
    ha_l = []
    g_l = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            ha_l.append(abs(pop[i].ha))
            g_l.append(abs(pop[i].gene))
    return(ha_l,g_l)

def list_ha_of_males_sp1(pop): # get a list of the center of sexual activity of males
    mal1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            mal1.append(abs(pop[i].ha))
    return(mal1)

def list_ha_of_females_sp1(pop): # get a list of the center of sexual activity of females
    fem1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "F":
            fem1.append(abs(pop[i].ha))
    return(fem1)


def list_wa_of_males_sp1(pop): # get a list of the window of sexual activity of males (deprecated)
    mal3 = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            mal3.append(abs(pop[i].wa)) 
    return(mal3)


def list_gene_of_males(pop): # get a list of the genome of males
    gene_data = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            gene_data.append(pop[i].gene) 
    return(gene_data)

def list_gene_of_females(pop): # get a list of the genome of females
    gene_data = []
    for i in range(len(pop)): 
        if pop[i].sex == "F":
            gene_data.append(pop[i].gene) 
    return(gene_data)

def sex_ratio_sp1(pop): # get the sex ratio of the population
    m = 0
    f = 0
    for i in range(len(pop)):
        if pop[i].sex == "M":
            m += 1
        elif pop[i].sex == "F":
            f += 1
        else:
            pass
    try :
        result = f/(m+f)
    except ZeroDivisionError:
        result = "Population extinct"
    return(result)


def size_sp1_sp2(pop): # get the size of the population 
    size_sp1 = 0
    for i in range(len(pop)):
        if pop[i].sp == 1:
            size_sp1 += 1
    return(size_sp1,len(pop)-size_sp1)

def list_gene(pop): # get a list of the genome of males
    gene_data = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            gene_data.append(pop[i].gene) 
    return(gene_data)

def list_e(pop): # get a list of the window of sexual activity of males (deprecated)
    mal3 = []
    for i in range(len(pop)): 
            mal3.append(abs(pop[i].emerg)) 
    return(mal3)
#################### simulation over several days with replicates  #####################


# Modify the parameter you want to test in the if (line 364) and else (line 371), as well as in the Output (line 444).

def ReplicatesSimulation(queue,qn,simulation_days, Output):     # Main function to loop through replications/parameters
    repli = queue.get()    
    print("Started replication number " + str(repli), flush=True)          
    if repli !=0 and repli%(replications-1)==0:     # When we're not done with the replications of a parameter
        n = qn.get()
        incompatibility_lim = start_step + n*step          
        qn.put(n+1)      
        queue.put(0)

    else :                                # When we have to move on to the next parameter
        n = qn.get()
        qn.put(n)      
        incompatibility_lim = start_step + n*step         
    
        queue.put(repli+1)        
    
    mean_ha_evolution_sp1 = []           # used to track the evolution of mean sexual activity timing
    mean_wa_evolution_sp1 = [] 
    ha_list = []                          # used to track the evolution of mean sexual activity window width (deprecated)
    male_ha_list = []
    female_ha_list = []
    pop_size_sp1 = []                    # used to track the evolution of population size
    total_lifespan = []
    total_generation = []
    generation_time = []
    survival_time = []
    mean_e_evolution_sp1 = []
    
    # initialisation of the population
    population = []                                               # list of all individuals
    initial_abund_1 = 100                                         # initial effective
                                           
    
    for j in range(initial_abund_1):                              # creation of the population
        initial_ha1 = 0.5
        ha = stats.truncnorm.rvs((lower - initial_ha1)/mu,(upper-initial_ha1)/mu, initial_ha1, mu)
        population.append(individuals(r.sample(["M","F"],1)[0],
                                ha,
                                None,
                                0.15,
                                None,
                                0,
                                1,
                                genome_size * [0],
                                [],
                                0,
                                0,
                                0,
                                None))
    
    
    
    
    
    for k in range(simulation_days): # for desired duration of replicates
        
        t = 0 # reset initial time to 0, corresponds to the start of the day
        
        emergence_list,emergence_schedule,generation_times = emergences_and_schedule(population,incompatibility_lim,ve1,e1,K,mu,k)                      # retrieves emergence list and emergence schedule of the day
        total_generation = generation_times
        ht,pop = patrolling_schedule(population)
        known_events = {**ht , **emergence_schedule }       # creates a dictionary of events (emergences and sexual activity shifts) supposedly happening on this day 
        known_events[1] = "end day"                          # adding the end of the day event to known_events
        chronological_events = sorted(known_events.keys())   # sort events of the day chronologically
        
        for i in range(len(chronological_events)):                          
            if known_events.get(chronological_events[i])[2] == "e":     # determines if "known" event is an emergence           
                population.append(emergence_list[int(known_events.get(chronological_events[i])[19:-1])]) # emerging individual is added to the population
    
            while t < chronological_events[i]:                          # until the next "known" event happens
                lifespan,t = event(population,t,chronological_events[i],dlt_c,mating_threshold,k)   # determine if an event happens and which event, or if no event happens between ti and T
                if lifespan !=None:
                    total_lifespan.append(lifespan)
            
        if t == 100:
            break
        
        generation_time.append(total_generation)
        survival_time.append(total_lifespan)
        mean_ha_evolution_sp1.append(np.mean(list_ha(population)))       # registers the mean sexual activity timing of males at the end of the day
        mean_e_evolution_sp1.append(np.mean(list_e(population)))  # deprecated
        pop_size_sp1.append(size_sp1_sp2(population)[0])                 # registers the population size at the end of the day
        ha_list.append(list_ha(population))            # registers male sexual activity timing
        male_ha_list.append(list_ha_of_males_sp1(population))        # registers female sexual activity timing

    
    
    Output.put([repli, 
                [pop_size_sp1],
                [mean_ha_evolution_sp1],
                [mean_e_evolution_sp1],
                [[list_ha_of_males_sp1(population)]],
                [list_e(population)],
                [ha_list],[male_ha_list],
                [sex_ratio_sp1(population)],
                [survival_time],
                list_gene(population),[generation_time], 
                population,
                incompatibility_lim,
                n])



simulation_days = 500    # Number of days for each replicate
param = 9                # Number of parameters tested
    

replications = 100       # Number of replications for each parameter
start_step = 0.1         # First parameter value 
step = 0.1               # Step size in parameter value

Test_value = "incompatibility"         # Put "e" if you want to test the effect of emergence values
print(Test_value) 

import multiprocessing
from time import perf_counter

if __name__ == "__main__":
    
    print("Number of CPUs available : " + str(multiprocessing.cpu_count()), flush=True)
    pool = multiprocessing.Pool(multiprocessing.cpu_count())       # Number of cores used
    
    m = multiprocessing.Manager()
    
    queue = m.Queue()                   # queue for the current replication number
    qn = m.Queue()                      # queue for the current parameter number
    Output = m.Queue()                  # Output queue for the data
    
    queue.put(0)                        # Initialize the queues at 0
    qn.put(0)
    
    start_time = perf_counter()         # Check run time for test purposes
    
    # Main loop :
    processes = [pool.apply_async(ReplicatesSimulation, args = (queue,qn,simulation_days, Output)) for x in range(param*replications)]
    
    result = [p.get() for p in processes]   
    finish_time = perf_counter()

    print("Program finished in " + str(finish_time-start_time) +"  seconds")
   
 
bufferpop = []                                        # Temporary variable to store all of the output data
for i in range(Output.qsize()):
    bufferpop.append(Output.get())
from operator import itemgetter
print("append to bufferpop")
popsort = sorted(bufferpop, key=itemgetter(13))        # Sort by replicate number 

print("All results obtained and sorted")

# Initialize output variable 
pop1 = []
ha_list1  = []
wa_list1  = []
ha_final_list1 = []
e_final_list1 = []
ha_list = []
male_ha_list = []
sex1 = []
generation_time = []
survival_time = []
population = []
beta = []
gene = []
gene_f = []


for i in range(replications*param):         # Sort all the output data into their variables
    pop1.append(popsort[i][1][0])
    ha_list1.append(popsort[i][2][0])
    wa_list1.append(popsort[i][3][0])
    ha_final_list1.append(popsort[i][4])
    e_final_list1.append(popsort[i][5])
    ha_list.append(popsort[i][6])
    male_ha_list.append(popsort[i][7])
    sex1.append(popsort[i][8][0])
    generation_time.append(popsort[i][9])
    gene.append(popsort[i][10])
    survival_time.append(popsort[i][11]) 
    population.append(popsort[i][12]) 
    beta.append(popsort[i][13])
    n = popsort[i][14]



def gaus(x,a,x0,sigma):
    """
    For gaussian curve fitting
    
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def std_binom(total,success):
    th_prob_success = success/total
    std = np.sqrt(th_prob_success * (1 - th_prob_success) / total)
    return std



print('Input data sorted')
val_n = qn.get()            # Get the number of parameters
for i in range(val_n):      # Generate data for each parameter
    
    mean_ha1 = 0
    mean_ha2 = 0
    mean_wa1 = 0
    mean_wa2 = 0
    mean_var2 = 0
    mean_var1 = 0
    mean_w_var1 = 0
    mean_w_var2 = 0
    count1 = 0
    count2 = 0
    cocount = 0
    cowa1 = 0
    cowa2 = 0
    coha1 = 0
    coha2 = 0
    # Get mean and variance of emergence peak and window duration
    for rep in range(i*replications,(i+1)*replications):
        if pop1[rep][-1] != 0:
            if len(ha_final_list1[rep][0][0]) != 0:
                mean_ha1 += np.mean(ha_final_list1[rep][0][0])
                # mean_wa1 += np.mean(wa_final_list1[rep][0])
                mean_var1 += np.var(ha_final_list1[rep][0][0])
                # mean_w_var1 += np.var(wa_final_list1[rep][0])
                count1 += 1

                    
        if pop1[rep][-1] != 0 :
            if len(ha_final_list1[rep][0][0]) != 0:
                # cowa1 += np.mean(wa_final_list1[rep][0])
                coha1 += np.mean(ha_final_list1[rep][0][0])   
                cocount += 1 
                
    if count1 == 0:
        mean_ha1 = "Extinct 100% of the time"
        # mean_wa1 = "Extinct 100% of the time"
        mean_var1 = "Extinct 100% of the time" 
        # mean_w_var1 = "Extinct 100% of the time"
    else: 
        mean_ha1 = mean_ha1/count1
        # mean_wa1 = mean_wa1/count1
        mean_var1 = mean_var1/count1
        # mean_w_var1 = mean_w_var1/count1

        
    if cocount == 0:
        # cowa1 = "Never coexisted"
        coha1 = "Never coexisted"
    else : 
        # cowa1 = cowa1/cocount
        coha1 = coha1/cocount
        
    print('incompatibility_lim: '+str(beta[i*replications]))
    
    print('ha1 : '+str(mean_ha1),
          # '\nwa1 : '+str(mean_wa1)
          )
    
    # print('coha1 : '+str(coha1),
    #       '\ncowa1 : '+str(cowa1))
    
    print('Variance final ha1 : ' + str(mean_var1))
    # print('Variance final wa1 : ' + str(mean_w_var1))

        

    # =====
    
    # Get extinction rate, mean population size and sex ratio
    taux_extinction_sp1 = 0
    taux_extinction_sp2 = 0
    pop_sp1 = 0
    pop_sp2 = 0
    sex_ratio_sp1 = 0
    sex_ratio_sp2 = 0
    
    for j in range(i*replications,(i+1)*replications):
        if len(ha_final_list1[j][0][0]) == 0:
            taux_extinction_sp1 +=1
        if len(ha_final_list1[j][0][0]) > 0:
            pop_sp1 += pop1[j][-1]
            sex_ratio_sp1 += sex1[j]
            
            
    print('Extinction rate : '+ str(taux_extinction_sp1/replications))
    
    if replications-taux_extinction_sp1 ==0:
        print('Average population size : Extinct')
    else :  
        print('Average population size : ' +str(pop_sp1/(replications-taux_extinction_sp1)))
        print('Sex ratio : ' +str(sex_ratio_sp1/(replications-taux_extinction_sp1)))

    
    moy_peak = []
    double_peaks = 0
    mean_fst = []
    peaks_distance = []
    fwmh_tot_1_peak = []
    fwmh_tot_2_peak = []
    generation_mean = []
    lifespan_mean = []
    linkage_des_percent = 0
    linkage_des_size1 = []
    linkage_des_percent2 = 0
    linkage_des_size2 = []
    E=0
    L=0
    I=0
    
    end_fv_1_peak = []
    end_fv_2_peak = []
    for rep in range(i*replications,(i+1)*replications):
        
        # Uncomment to save figures of simulation outcomes
        
        gen = []
        for g in range(len(ha_list[rep][0])):
            gen.extend([g]*len(ha_list[rep][0][g]))
            
        yCM = np.array([j for i in ha_list[rep][0] for j in i])
        # x = np.array(gen)
        
        plt.figure()
        plt.ylabel("Daily time")
        plt.xlabel("Day")
        img = plt.hist2d(gen, yCM, bins = (simulation_days,60))
        plt.savefig("HIST2D_"+str(i)+str(rep)+".png")
        
              
        #Detecting the number of peaks in the distribution of the timing of sexual activity
        #Detection the position of the peaks to classify them as early or late
        peaks_nb = []
        fv = []
        # plt.close('all')

        
        generation_mean.append(np.mean(generation_time[rep]))
        # lifespan_mean.append(np.mean(survival_time[rep]))



        peaks_x = []
        peaks_y = []
        cluster_length = []
        peak1_x = []
        peak1_y = []
        currE = 0
        currL = 0
        Elen = []
        Llen = []
        lastpeakE = 0
        lastpeakL = 0
        peaks = []
        if len(ha_list[rep][0])==simulation_days:
            
            plt.figure()
            colors = cm.viridis(np.linspace(0,1,simulation_days))
            
            for j in range(len(ha_list[rep][0])):

                MF_ha_list = ha_list[rep][0][j]
                fv = []
    

                #Finding peaks with kdeplot
                for k in range(len(MF_ha_list)):
                    fv.append(MF_ha_list[k])
                   
                sns.kdeplot(np.array(fv),bw_method = 0.15, clip = [0,1],color = "lightgray", legend = 'day' + str(j))
                plt.xlim(0,1)
                
                x = plt.gca().lines[-1].get_xdata() # Get the x data of the distribution
                y = plt.gca().lines[-1].get_ydata() # Get the y data of the distribution


                peaks = find_peaks(y, prominence = 0.5, distance = 30)[0]
                
                
                #Separating unimodal and bimodal cases
                
                    
                    
                if len(peaks)==2 and abs(x[peaks[0]]-x[peaks[1]]) > 0.15:

                    peaks_x.append(x[peaks])
                    peaks_y.append(j)
                    if j == len(ha_list[rep][0])-1:
                        double_peaks += 1
                        
                        print('Early sub-pop sexual activity timing : ' + str(x[peaks][0])) 

                        print('Late sub-pop sexual activity timing : ' + str(x[peaks][1]))  
                        end_fv_2_peak.append(fv)
                        
                        
                elif len(peaks)==2 and abs(x[peaks[0]]-x[peaks[1]]) <= 0.15:
                    peaks = np.array([round(np.mean(peaks))])

                    
                if len(peaks)==1:
                    peak1_x.append(x[peaks])
                    peak1_y.append(j)
                    if j == len(ha_list[rep][0])-1:
                        end_fv_1_peak.append(fv)
                
                    
                    
                peaks_nb.append(len(peaks))
            


            key1 = peak1_y
            value1 = peak1_x
            key2 = peaks_y
            value2 = peaks_x
            
            dicpeak1 = {}
            for key in key1:
                for value in value1:
                    dicpeak1[key] = value.tolist()
                    value1.remove(value)
                    break
        
            dicpeak2 = {}
            for key in key2:
                for value in value2:
                    dicpeak2[key] = value.tolist()
                    value2.remove(value)
                    break
                        
                
            #Finding sub-population coexistence length
            dictfused = dicpeak2 | dicpeak1
            iterdic2 = iter(dicpeak2)
            iterdic1 = iter(dicpeak1)
            lastpeak = 0
            lastpeakE = 0
            lastpeakL = 0
            late = 0
            early = 0
            CLate = []
            CEarly = []
            for d in np.sort(list(dictfused.keys())):
                    
                if len(dictfused[d]) == 2:
                    c = next(iterdic2)

                    if abs(dicpeak2[d][1]-dicpeak2[c][1])<0.3 and d-lastpeakL<15:
                        late += 1
                        lastpeakL = d

                    if abs(dicpeak2[d][0]-dicpeak2[c][0])<0.3 and d-lastpeakE<15:
                        early += 1
                        lastpeakE = d 
                        
                if len(dictfused[d]) ==1:   
                    if dictfused[d][0] < 0.3 :
                        
                        early += 1
                        lastpeakE = d
                        lastpeakL +=1

                        if late>40:
                            CLate.append(late)
                            late = 0
                            lastpeakL = 0

                    if dictfused[d][0] > 0.3 :
                        late += 1
                        
                        lastpeakE +=1
                        lastpeakL = d

                        if early>40:
                            CEarly.append(early)
                            early = 0
                            lastpeakE = 0

            if late > 40:
                CLate.append(late)
                
            if early > 40:
                CEarly.append(early)
                
            
                
                
            print("Early sub-pop lenght, Late sub-pop lenght : ")  
            print(CEarly, CLate)    
            
            
            
            if len(peaks) == 1:
                
                if Test_value == "e":
                    
                    if beta[i*replications]-0.05 < x[peaks][0] < beta[i*replications]+0.05 :
                            
                        print("Immediate population only, average sexual activity time :" + str(x[peaks][0]))
                        I+=1
                        
                    elif x[peaks][0] > beta[i*replications]+0.05 :
                        print('Delayed-dusk population only, average sexual activity time : ' + str(x[peaks][0]))  
                        L+= 1
                        
                    elif x[peaks][0] < beta[i*replications]-0.05 :
                        print('Delayed-dawn population only, average sexual activity time : ' + str(x[peaks][0]))  
                        E+= 1
                    
                # Errors are normal here if test value = e
                else :
                     
                      if e1-0.05 < x[peaks][0] <e1+0.05 :
                              
                          print("Immediate population only, average sexual activity time :" + str(x[peaks][0]))
                          I+=1
                          
                      elif x[peaks][0] > e1+0.05 :
                          print('Delayed-dusk population only, average sexual activity time : ' + str(x[peaks][0]))  
                          L+= 1
                          
                      elif x[peaks][0] < e1-0.05 :
                          print('Delayed-dawn population only, average sexual activity time : ' + str(x[peaks][0]))  
                          E+= 1
                
                # Fitting gaussian curve
                n = len(x)                          #the number of data points
                mean = sum(x*y)/n                   
                sigma = sum(y*(x-mean)**2)/n        
                
                try:
                    popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
                    fwmh = abs(2*np.sqrt(2*np.log(2))*popt[-1])
                    
                    half_max = gaus(popt[1], *popt) /2
                    l_half_max = popt[1] - fwmh/2
                    r_half_max = popt[1] + fwmh/2
                    
                    plt.plot(x,y,'b+:',label='data')
                    plt.plot(x,gaus(x,*popt),'ro:',label='fit')
                    plt.vlines(x=l_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.vlines(r_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.hlines(half_max, l_half_max, r_half_max, linestyles='--')
                    # plt.show()
                    plt.savefig("Gaussian_Width_"+str(i)+str(rep)+".png")
                    print('Full width at half maximum of the peak : ' + str(fwmh))
                    
                    fwmh_tot_1_peak.append(fwmh)
                
                except RuntimeError:
                    print("Could not fit curve")
                    
                
                
            elif len(peaks) == 2:
                        
                
                x1 = np.array([])
                y1 = np.array([])
                x2 = np.array([])
                y2 = np.array([])
                

                mid_peaks = x[min(range(len(y[min(peaks):max(peaks)])), key=y[min(peaks):max(peaks)].__getitem__)+min(peaks)]
                
                for p in range(len(x)):
                    
                    if x[p] < mid_peaks : 
                        x1 = np.append(x1,x[p])
                        y1 = np.append(y1,y[p])
                    
                    else:
                        x2 = np.append(x2,x[p])
                        y2 = np.append(y2,y[p])
                    
                n = len(x1)                          
                mean = sum(x1*y1)/n                   
                sigma = sum(y1*(x1-mean)**2)/n    
                
                try:
                    popt,pcov = curve_fit(gaus,x1,y1,p0=[1,mean,sigma])
                    
                    
                    fwmh = abs(2*np.sqrt(2*np.log(2))*popt[-1])
                    
                    half_max = gaus(popt[1], *popt) /2
                    l_half_max = popt[1] - fwmh/2
                    r_half_max = popt[1] + fwmh/2
                    
                    plt.plot(x1,y1,'b+:',label='data')
                    plt.plot(x1,gaus(x1,*popt),'ro:',label='fit')
                    plt.vlines(x=l_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.vlines(r_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.hlines(half_max, l_half_max, r_half_max, linestyles='--')
                    
                    fwmh_tot_2_peak.append(fwmh)
                    
                    print('Full width at half maximum of the early peak : ' + str(fwmh))
                    
                    n = len(x2)                          
                    mean = sum(x2*y2)/n                   
                    sigma = sum(y2*(x2-mean)**2)/n  
                    popt,pcov = curve_fit(gaus,x2,y2,p0=[1,mean,sigma])
                    
                    fwmh = abs(2*np.sqrt(2*np.log(2))*popt[-1])
                    
                    half_max = gaus(popt[1], *popt) /2
                    l_half_max = popt[1] - fwmh/2
                    r_half_max = popt[1] + fwmh/2
                    
                    
                    plt.plot(x2,y2,'b+:',label='data')
                    plt.plot(x2,gaus(x2,*popt),'ro:',label='fit')
                    plt.vlines(x=l_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.vlines(r_half_max,ymin=0,ymax=half_max, linestyles='--')
                    plt.hlines(half_max, l_half_max, r_half_max, linestyles='--')
                    # plt.show()
                    plt.savefig("Gaussian_Width_"+str(i)+str(rep)+".png")
                    
                    fwmh_tot_2_peak.append(fwmh)
                    
                    print('Full width at half maximum of the late peak : ' + str(fwmh))
                    
                    peaks_distance.append(x[max(peaks)]-x[min(peaks)])
                    
                    print('Distance between peaks : '+str(x[max(peaks)]-x[min(peaks)]))
                
                except RuntimeError:
                    print("Could not fit curve")
                    
                
    
            if len(peaks_x)>39:
                temp = 1
                for i in range(len(peaks_x)-1):
                    if abs(peaks_x[i][0] - peaks_x[i+1][0])<0.15 and peaks_y[i+1]-peaks_y[i]<15:
                        temp+=1
                        
                    else :
                        cluster_length.append(temp)
                        temp = 1
                cluster_length.append(temp)
            for i in range(len(cluster_length)-1,-1,-1):
                if cluster_length[i]<40:
                    del cluster_length[i]
                 

            moy_peak.append(peaks_nb)  
            peaks_track = moy_peak
    
            # FST computing
            
            
            # get gene values and ha values
            dic = {}
            

            for j,k in enumerate(male_ha_list[rep][0][-1]):
                dic[k] = gene[rep][j]

            
            # get median ha
            ha_mid = 0
    
            if len(peaks)==2 and abs(x[peaks[0]]-x[peaks[1]]) > 0.15:
                
                
                # ha_mid = (x[peaks[0]] + x[peaks[1]])/2
                ha_mid = x[min(range(len(y[min(peaks):max(peaks)])), key=y[min(peaks):max(peaks)].__getitem__)+min(peaks)]
                print("HA MID : " + str(ha_mid))
                
            ha_sup = []
            ha_inf = []
            linkage_des_size = 0
                
            ha_sup = np.array([x for key, x in dic.items() if key > ha_mid])
            ha_inf = np.array([x for key, x in dic.items() if key <= ha_mid])
            
            
            
            
            
            

            # FST calculation
                
    
                        
            size_sup = len(ha_sup)
            size_inf = len(ha_inf)
                
            if size_sup*size_inf !=0:
                found_linkage = False
            
                # Linkage desequilibrium test
                group_0 = []
                group_1 = []
                
                for loci in range(genome_size):
                    for ind in range(len(population[rep])):
                        #if population[rep][ind].sex == "M": uncomment for faster sim with similar results
                            if population[rep][ind].gene[loci] == 0:
                                group_0.append(population[rep][ind].ha)
                            else:
                                group_1.append(population[rep][ind].ha)

                    t_stat, p_value_ttest = stats.ttest_ind(group_0, group_1)
                    
                    if p_value_ttest <0.05:
                        linkage_des_size += 1
                        
                        if found_linkage == False:
                            linkage_des_percent += 1
                            found_linkage=True
                  
                if linkage_des_size !=0:
                    linkage_des_size1.append(linkage_des_size)
                

                
                
                df = 0
                for sup in range(size_sup):
                    for inf in range(size_inf):
                        df += scipy.spatial.distance.hamming(ha_sup[sup],ha_inf[inf], w=None)
                
                df = df/(size_sup*size_inf)
                
                
                
                dfsup = 0
                for sup in range(size_sup):
                    for sup2 in range(size_sup):
                        dfsup += scipy.spatial.distance.hamming(ha_sup[sup],ha_sup[sup2])
                
                dfsup = dfsup/(size_sup**2)
                
                dfinf = 0
                for inf in range(size_inf):
                    for inf2 in range(size_inf):
                        dfinf += scipy.spatial.distance.hamming(ha_inf[inf],ha_inf[inf2])
                
                dfinf = dfinf/(size_inf**2)
                
                
                fstinf = max((df - dfsup)/df,0)
                fstsup = max((df - dfinf)/df,0)
                
        
                mean_fst.append((fstsup+fstinf)/2)
                print("FST"+str(rep)+" : " + str((fstsup+fstinf)/2))
                
                
                
    plt.close('all')
    plt.figure()
    
    for f in range(len(end_fv_1_peak)):
        sns.kdeplot(
                np.array(end_fv_1_peak[f]),
                bw_method=0.15,
                clip=[0, 1],
                color="lightgray",  # Overlay all distributions in light gray
                alpha=0.7,
                legend=False        # Avoid multiple legends
            )
    
    plt.xlabel("Activity Timing", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("Overlapping KDE Plots of Simulation Endings", fontsize=14)
    plt.savefig("combined_1peak_kdeplot.png", dpi=300)  
    
    plt.close('all')
    plt.figure()
    
    for f in range(len(end_fv_2_peak)):
        sns.kdeplot(
                np.array(end_fv_2_peak[f]),
                bw_method=0.15,
                clip=[0, 1],
                color="lightgray",  # Overlay all distributions in light gray
                alpha=0.7,
                legend=False        # Avoid multiple legends
            )
    
    plt.xlabel("Activity Timing", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title("Overlapping KDE Plots of Simulation Endings", fontsize=14)
    plt.savefig("combined_2peak_kdeplot.png", dpi=300)  
           
    print("Mean FST for this parameter value : "+ str(np.mean(mean_fst)))
    print("Standart deviation of the Fst values : " + str(np.var(mean_fst)))
    total = L + E + double_peaks
    print("==========================================================================")
    if double_peaks!=0:
        print("Percentage of simulation with linkage desequilibrium : " + str(round((linkage_des_percent/double_peaks) * 100,2))+"% of the time, var : " + str(std_binom(replications,linkage_des_percent)))  
        print("Average number of loci with linkage desequilibrium : " + str(np.mean(linkage_des_size1)))
        print("Standard deviation of the number of loci with linkage desequilibrium : " + str(np.std(linkage_des_size1)))
        print(linkage_des_size1)
    print("==========================================================================")
    print("Only delayed-dusk " + str((L/replications)*100)+"% of the time, var : " + str(std_binom(replications,L)))  
    print("Only delayed-dawn : " + str((E/replications)*100)+"% of the time, var : " + str(std_binom(replications,E)))  
    print("Only immediate : " + str((I/replications)*100)+"% of the time, var : " + str(std_binom(replications,I)))  
    print("Bimodal case " + str(((double_peaks/replications))*100)+"% of the time, var : " + str(std_binom(replications,double_peaks)))  
    print("==========================================================================")
    print("Average distance between peaks : " + str(np.mean(peaks_distance)))
    print("Standart deviation of the distance between peaks : " + str(np.std(peaks_distance)))
    print("==========================================================================")
    print("Average width of double peaks : " + str(np.mean(fwmh_tot_2_peak)))
    print("Standard deviation of the double width : " + str(np.std(fwmh_tot_2_peak)))
    print("==========================================================================")
    print("Average width of single peaks : " + str(np.mean(fwmh_tot_1_peak)))
    print("Standard deviation of the single width : " + str(np.std(fwmh_tot_1_peak)))
    print("==========================================================================")
    # print('Average life time : '+str(np.mean(lifespan_mean)))
    # print('Std life time : '+str(np.std(lifespan_mean)))
    print("==========================================================================")
    print('Average generation time : '+str(np.mean(generation_mean)))
    print('Std generation time : '+str(np.std(generation_mean)))
    # print("Survived for " +str(np.mean(len(ha_list[rep]))) + " days")
    print('\n')
    
