# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 18:43:50 2024

@author: adm.violaine
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 11:02:00 2023

@author: morpho2019
"""

import random as r
import numpy as np
import scipy
import scipy.stats as stats
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from time import perf_counter
# import time
# plt.switch_backend('agg')
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 09:43:27 2022

@author: bouinier
"""


########################## model global parameters ############################


beta = 3                     # encounter rate between patrolling males and virgin females independant of species
d = 0.1                       # death rate applied to whole population
dlt_e = 0.1                  # energy linked death rate of patrolling individuals
dlt_c = 0.1                # density dependant competition death rate patrolling male
dlt_l = 0.1
offspring_number = 5          # mean offspring number per female per day 
mu = 1/25                     # mutation size of patrolling window position AND width (unimplemented for width)

lower = 0                     # used to bound drawing in truncated normal distributions
upper = 1                     # used to bound drawing in truncated normal distribution

K = 5000                      # carrying capacity regarding competition between females for access to host plants to lay eggs

e1 = 0.5                     # peak emergence time species 1

ve1 = 0.12                    # emergence time variance species 1

# initial_wa1 = 0.05            # initial patrol time width sp 1

mu_g = 0.1                    # mutation rate for genes
genome_size = 10              # number of genes of each individual
incompatibility_lim = 0.4     # maximum genetic distance allowing breeding
#
######################### individuals attributes ##############################


class individuals:                           # Individuals in our population belongs to this class
    def __init__(self, sex, ha, ha_male, wa, wa_male, f, sp, gene, gene_male, nb_offsprings,emerg): # Each individual is defined by the following attributes :
        self.sex = sex                       # Male or Female
        self.ha = ha                         # patrolling window position
        self.ha_male = ha_male               # patrolling window position storage for mating events
        self.wa = 0.15                       # patrolling window width
        self.wa_male = wa_male               # patrolling window width storage for mating events
        self.f = f                           # fertilisation state 0 (virgin), 1 (mated)
        self.sp = sp                         # deprecated
        self.gene = gene                     # np.array([a,b,c,d]) where a,b,c,d are between 0 and 1
        self.gene_male = gene_male           # np.array([a,b,c,d]) where a,b,c,d are between 0 and 1 (for females to save the male genes when mated)                  
        self.nb_offsprings = nb_offsprings   # number of offsprings birthed for a female
        self.emerg = emerg                   # emergence time of the individual
        
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

def list_virgin_patrolling_females(pop,ti): # get a list of patrolling virgin females in the population
    fem_v = []
    for i in range(len(pop)): 
        # print(pop[i].ha - pop[i].wa)
        if pop[i].sex == "F" and pop[i].f == 0 and ((pop[i].ha - pop[i].wa) <= ti < (pop[i].ha + pop[i].wa)):
            fem_v.append(pop[i])

    return(fem_v)

def list_virgin_nonpatrolling_females(pop,ti): # get a list of patrolling virgin females in the population
    fem_v = []
    for i in range(len(pop)): 
        if pop[i].sex == "F" and pop[i].f == 0 and ((pop[i].ha - pop[i].wa) > ti 
                                    or ti >= (pop[i].ha + pop[i].wa)):
            fem_v.append(pop[i])
    return(fem_v)

def list_mated_females_1(pop): # get a list of mated females within species 1 in the population
    fem_f1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "F" and pop[i].f == 1 and pop[i].sp == 1:
            fem_f1.append(pop[i])
    return(fem_f1)


def list_patrolling_males_1(pop,ti): # get a list of patrolling males in sp 1 at time "ti" in the population
    mal_a1 = []
    
    for i in range(len(pop)): 
        # print(pop[i].ha - pop[i].wa)
        if pop[i].sex == "M" and ((pop[i].ha - pop[i].wa) <= ti < (pop[i].ha + pop[i].wa)):
            mal_a1.append(pop[i])
    return(mal_a1)


def list_nonpatrolling_males(pop,ti): # get a list of non-patrolling males at time "ti" in the population, unregarding the species
    mal_i = []
    for i in range(len(pop)): 
        if pop[i].sex == "M" and ((pop[i].ha - pop[i].wa) > ti 
                                    or ti >= (pop[i].ha + pop[i].wa)):
            mal_i.append(pop[i])
    return(mal_i)

def list_eggs_weight(pop): 
    egg_weight= []
    for i in range(len(pop)): 
        # if pop[i].nb_offsprings !=0:
            egg_weight.append(pop[i].nb_offsprings)
    return(egg_weight)

########################### time period cutting functions #####################

def patrolling_schedule(pop):          # get a dictionary of hours of transition from non-patrolling to patrolling and vice versa for males of both sp
    ht = {}                                  # dictionnaries allows us to link events with their associated time
    for i in range(len(pop)): 
        
        if pop[i].emerg > pop[i].ha + pop[i].wa:   #patrol is before emergence, skip individual
            pop[i].emerg=0       
            continue
        if pop[i].ha - pop[i].wa < pop[i].emerg < pop[i].ha + pop[i].wa:   #emergence is during patrol
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



def emergences_and_schedule(pop):   
    Mated_females = list_mated_females_1(pop)            
    b = offspring_number * (1-len(Mated_females)/K) # number of offspring per female per day is density dependant depending on mated females 
    if b < 0 :
        b = 0
    emergence_list = []                                         # list of all emergences of the day
    emergence_schedule = {}                                     # dict of emergence schedule of offsprings
    for i in range(len(Mated_females)): 
        # if pop[i].sex == "F" and pop[i].f == 1:
            if scipy.spatial.distance.hamming(Mated_females[i].gene, Mated_females[i].gene_male, w=None) <= incompatibility_lim:
                offsp = np.random.poisson(b) 
                genes = g_transfer(Mated_females[i].gene, Mated_females[i].gene_male)
                for j in range(offsp):   
                    emergence_list.append(individuals(r.sample(["M","F"],1)[0],         # creation of new emerging individuals
                                                  stats.truncnorm.rvs((lower - (Mated_females[i].ha+Mated_females[i].ha_male)/2)/mu,(upper-(Mated_females[i].ha+Mated_females[i].ha_male)/2)/mu, (Mated_females[i].ha+Mated_females[i].ha_male)/2, mu), # emerging individual inherit the patrolling window position of their father
                                                  None,
                                                  0.15,
                                                  None,
                                                  0,                                                                               # emerging individual are virgins
                                                  Mated_females[i].sp,            # emerging individual belong to the species of the parents
                                                  genes,
                                                  genome_size * [0],
                                                  0,
                                                  stats.truncnorm.rvs((lower - e1)/ve1,(upper-e1)/ve1, e1, ve1)
                                                  ))
                                                  
                    Mated_females[i].nb_offsprings += 1
                    
    # for k in range(len(emergence_list)):                        # drawing of emergence time for each new individual
    #     emergence_schedule[stats.truncnorm.rvs((lower - e1)/ve1,(upper-e1)/ve1, e1, ve1)] = "emergence individu " + str(k) + emergence_list[k].sex
    for k in range(len(emergence_list)):                        # drawing of emergence time for each new individual
        emergence_schedule[emergence_list[k].emerg] = "emergence individu " + str(k) + emergence_list[k].sex
    
    return(emergence_list,emergence_schedule)                   # returns emergence list AND emergence with associated time

########################### event drawing #####################################


def event(pop,ti, T):                                                   # determine if an event happens between ti and T, and chooses which one
    nb_patrolling_males_1 = len(list_patrolling_males_1(pop,ti))        # calculating each population category abundance
    nb_nonpatrolling_males = len(list_nonpatrolling_males(pop,ti))
    nb_virgin_patrolling_females = len(list_virgin_patrolling_females(pop,ti))
    nb_virgin_nonpatrolling_females = len(list_virgin_nonpatrolling_females(pop,ti))
    nb_mated_females_1 = len(list_mated_females_1(pop))
    Etot = np.sum(list_eggs_weight(pop))

    Ei = []
    # print(nb_mated_females_1)
    lr = beta * (nb_patrolling_males_1) * nb_virgin_patrolling_females         # calculating lambda_r and lambda_d used in exponential law
    ld = d * (nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1 + nb_nonpatrolling_males) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot
    # print((lr/(lr+ld)) + (d * (nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1 + nb_nonpatrolling_males) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot)/(lr+ld))
    # print((lr/(lr+ld)) + (d * (nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1 + nb_nonpatrolling_males) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot)/(lr+ld))
    
    if lr + ld == 0:                                            # specific case if population go extinct
        print("population extinct")                             # return t = 100 to indicate error
        return(100)
    draw = np.random.exponential(1/(lr + ld))                 # drawing in exponential law of parameter lambda = lambda_r + lambda_d
    if draw < T - ti:
        uni = np.random.uniform(0,1)                            # drawing in unif(0,1) to set which event is happening
        
        if uni < lr/(lr+ld):                                  # event : meeting between patrolling male and virgin female
            father = r.sample(list_patrolling_males_1(pop,ti), 1)[0]  # father is randomly sampled among patrolling males
            mother = r.sample(list_virgin_patrolling_females(pop,ti),1)[0]    # mother is randomly sampled among virgin females
            if father.sp == mother.sp:
                mother.f = 1                                      
                mother.ha_male = father.ha                        # mother stores father's patrolling window position
                mother.wa_male = father.wa                        # mother stores father's patrolling window width
                mother.gene_male = father.gene                    # mother stores father's genome
            
        elif uni < (lr/(lr+ld)) + (d * nb_nonpatrolling_males)/(lr + ld):       # event : death of non-patrolling male
            dead = r.sample(list_nonpatrolling_males(pop,ti),1)[0]              # the dead is randomly chosen among non-patrolling males at time ti
            pop.remove(dead)
            
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_patrolling_females) + dlt_e * nb_virgin_patrolling_females)/(lr+ld):    # event : death of virgin patrolling female  
            dead = r.sample(list_virgin_patrolling_females(pop,ti),1)[0]    
            pop.remove(dead)
        
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females) + dlt_e * nb_virgin_patrolling_females)/(lr+ld): # event : death of virgin nonpatrolling female    
            dead = r.sample(list_virgin_nonpatrolling_females(pop,ti),1)[0]                                      
            pop.remove(dead)
        
        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1) + dlt_e * nb_virgin_patrolling_females + dlt_l * Etot)/(lr+ld):
            Esum = 0
            E = list_mated_females_1(pop)
            for l in range(nb_mated_females_1):
                Ei = E[l]
                Esum += Ei.nb_offsprings 
        
                if uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females) + d*(l+1) + dlt_e * nb_virgin_patrolling_females  + dlt_l*Esum)/(lr+ld): # the dead is chosen among mated females of species 1
                    pop.remove(Ei)
                    break

        elif uni < (lr/(lr+ld)) + (d * (nb_nonpatrolling_males + nb_virgin_nonpatrolling_females + nb_virgin_patrolling_females + nb_mated_females_1 + nb_patrolling_males_1) + dlt_e * (nb_patrolling_males_1 + nb_virgin_patrolling_females) + dlt_c * nb_patrolling_males_1 * (nb_patrolling_males_1 -1) + dlt_l * Etot)/(lr+ld):
                dead = r.sample(list_patrolling_males_1(pop,ti),1)[0]
                pop.remove(dead) 
                    
        return(ti + draw)
    else:
        return(T)                                               # no event happens when the draw in exponential law falls after T - ti interval


######################## simulation study #####################################


def list_ha_of_males_sp1(pop): # get a list of patrolling window position of males of sp 1
    mal1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            mal1.append(abs(pop[i].ha))
    return(mal1)

def list_ha_of_females_sp1(pop): # get a list of patrolling window position of females of sp 1
    fem1 = []
    for i in range(len(pop)): 
        if pop[i].sex == "F":
            fem1.append(abs(pop[i].ha))
    return(fem1)


def list_wa_of_males_sp1(pop): # get a list of patrolling window width of males of sp 1
    mal3 = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            mal3.append(abs(pop[i].wa)) 
    return(mal3)


def list_gene_of_males(pop): # get a list of patrolling window width of males of sp 2
    gene_data = []
    for i in range(len(pop)): 
        if pop[i].sex == "M":
            gene_data.append(pop[i].gene) 
    return(gene_data)

def list_gene_of_females(pop): # get a list of patrolling window width of males of sp 2
    gene_data = []
    for i in range(len(pop)): 
        if pop[i].sex == "F":
            gene_data.append(pop[i].gene) 
    return(gene_data)

def sex_ratio_sp1(pop): # get the sex ratio of sp1
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
        result = "pop 1 extinct"
    return(result)


def size_sp1_sp2(pop):
    size_sp1 = 0
    for i in range(len(pop)):
        if pop[i].sp == 1:
            size_sp1 += 1
    return(size_sp1,len(pop)-size_sp1)



#################### simulation over several days with replicates  #####################

                   

m_mean_ha_evolution_sp1 = []           # used to track the evolution of mean patrolling window position of sp1
f_mean_ha_evolution_sp1 = []           # used to track the evolution of mean patrolling window position of sp1

mean_wa_evolution_sp1 = []           # used to track the evolution of mean patrolling window width of sp 1
male_ha_list = []
female_ha_list = []
pop_size_sp1 = []                    # used to track the evolution of population size of sp 1
sex_evo = []

# initialisation of the population
population = []                                               # list of all individuals
initial_abund_1 = 100                                         # initial effective of sp 1
                                       

for j in range(initial_abund_1):                              # creation of sp 1 
    # initial_ha1 = np.random.uniform(0,1)                      # initial patrolling window position is random for each individual
    initial_ha1 = 0.5 
    population.append(individuals(r.sample(["M","F"],1)[0],
                            stats.truncnorm.rvs((lower - initial_ha1)/mu,(upper-initial_ha1)/mu, initial_ha1, mu),
                            None,
                            0.15,
                            None,
                            0,
                            1,
                            genome_size * [0],
                            genome_size * [0],
                            0,
                            0))




start_time = perf_counter()
for k in range(500): # for desired duration of replicates
    print(k)
    t = 0 # reset initial time to 0, corresponds to the start of the day
    
    emergences_today = emergences_and_schedule(population)                        # retrieves emergence list and emergence schedule of the day
    ht,pop = patrolling_schedule(population)
    known_events = {**ht , **emergences_today[1] }       # creates a dictionary of events (emergences and patrol shifts) supposedly happening on this day 
    known_events[1] = "end day"                                                   # adding the end of the day event to known_events
    chronological_events = sorted(known_events.keys())                            # sort events of the day chronologically

    for i in range(len(chronological_events)):                          
        if known_events.get(chronological_events[i])[2] == "e":                   # determines if "known" event is an emergence           
            population.append(emergences_today[0][int(known_events.get(chronological_events[i])[19:-1])]) # emerging individual is added to the population

        while t < chronological_events[i]:                                        # until the next "known" event happens
            # print(t)
            t = event(population,t,chronological_events[i])                       # determine if an event happens and which event, or if no event happens between ti and T
        
        
    if t == 100:
        break
    m_mean_ha_evolution_sp1.append(np.mean(list_ha_of_males_sp1(population)))       # registers the mean patrolling window position of males of sp 1 at the end of the day
    f_mean_ha_evolution_sp1.append(np.mean(list_ha_of_females_sp1(population)))       # registers the mean patrolling window position of males of sp 1 at the end of the day
    
    mean_wa_evolution_sp1.append(np.mean(list_wa_of_males_sp1(population)))       # registers the mean patrolling window width of males of sp 1 at the end of the day
    pop_size_sp1.append(size_sp1_sp2(population)[0])                                          # registers the population size at the end of the day
    male_ha_list.append(list_ha_of_males_sp1(population))
    female_ha_list.append(list_ha_of_females_sp1(population))
    sex_evo.append(sex_ratio_sp1(population))
        
gen = []
for g in range(len(female_ha_list)):
    gen.extend([g]*len(female_ha_list[g]))
yCM = np.array([j for i in female_ha_list for j in i])
plt.hist2d(gen, yCM, bins = (100,60))


print("done")
plt.figure()
plt.plot(pop_size_sp1)
finish_time = perf_counter()
print("Program finished in " + str(finish_time-start_time) +"  seconds")
