#About the project 

The python code "DailyStochModel.py"  simulates a stochastic, individual-based birth-and-death process (e.g. M ́el ́eard, 2016; Otto and Day, 2011), to study the evolution of daily periods of sexual activity timings within species, and its potential effect on population divergence and speciation.

There is also a variant modelling the evolution of seasonal non-perennial species, using the similar equations for mating and death events. ("SeasonalStochModel.py")

#How to use

If you wish to run it using the multiprocessing feature, a linux build is recommanded, or a cluster as it is ressource intensive.

To test the effet of one parameter, you must comment it at the start of the script, and modifiy the parameter later as explained at the start of the script.

This code will return text based results, made to be outputed by a cluster. The results will only be displayed on the console if run localy ! 

#Data

All of the resulting data is available in the "Data and graph" folder, along with R code to generate the figures of the main text and supplementary graphs.
