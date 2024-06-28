#About the project 

This python code simulates a stochastic, individual-based birth-and-death process (e.g. M ́el ́eard, 2016; Otto and Day, 2011), to study the evolution of daily periods of sexual activity timings within species, and its potential effect on population divergence and speciation.

There is also a variant modelling the evolution of seasonal non-perennial species, using the similar equations for mating and death events.

#How to use
If you wish to run it using the multiprocessing feature, a linux build is recommanded, or a cluster as it is ressource intensive.

To test the effet of one parameter, you must comment it at the start of the script, and modifiy the parameter at lines 371 and 378 (parameter_to_change = start_step + n*step), and line 451 in the Output element.

This code will output text based results, made to be outputed by a cluster. The results will only be displayed on the console if run localy ! 
