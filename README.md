#About the project 

The python code "Daily_Stochastic_Model.py"  simulates a stochastic, individual-based birth-and-death process (e.g. Méléard, 2016; Otto and Day, 2011), to study the evolution of daily periods of sexual activity timings within species, and its potential effect on population divergence and speciation.


Also included are multiple variation of the base daily model :
-DailyStochModel_sex_specific.py : Displaying the effects of sex-specific expression of the reproductive activity timing (ha)

-DailyStochModel_coevolution.py : Displaying the effects of allowing coevolution of the reproductive activity timing (ha) and the emergence timing (e)

-DailyStochModel_singlelociG : Displaying the effects of using a single continuous loci to represent the neutral genetic information, instead of a list of discrete loci.

-DailyStochModel_logisticG.py :  Displaying the effects of using a logistic function to model the survival chance of offspring rather than a threshold, to model a more continuous effect of incompatibilites.

There is also a variant modelling the evolution of seasonal non-perennial species, using the similar equations for mating and death events. ("SeasonalStochModel.py")

#How to use

If you wish to run it using the multiprocessing feature, a linux build is recommended, or a cluster as it is ressource intensive.

To test the effet of one parameter, you must comment it at the start of the script, and modify the parameter later as explained at the start of the script.

This code will return text based results, made to be outputted by a cluster. The results will only be displayed on the console if run locally ! 

#Data

All of the resulting data is available in the "Data and graph" folder, along with R code to generate the figures of the main text and supplementary graphs.
