README file for using the code accompanying "Estimating the ecological drivers of insect abundance when detection is imperfect"

##--------------------------------------------------------------------------
## HOW TO RUN AN ANLALYSIS ON REAL DATA

We collected mark-recapture data for 8 species of wild bees from 10 field sites in summer 2022 and 2023. See the metadata "METADATA.txt" for more info about the real data that we collected "./data/abundance_data.csv". The mark-recapture data could be treated as mark-recapture or devolved into simple counts - how many individuals were detected regardless of whether they were marked.

A multinomial N-mixture model can be fit to the mark-recapture data - "./multimix/models/multimix_negbin_NAs.stan". The model considers abundance to be drawn from a negative binomial distribution. The model simultaneously estimates abundance while estimating detection rate from the encounter histories. The model can also take a vector of NA data, which will tell the model to treat outcomes as NAs when we didn't conduct a survey. Use "./multimix/run_model/run_model.R" to fit this model to the data.

A binomial N-mixture model can be fit to the devolved count data - "./binmix/models/binmix_negbin_NAs.stan". The model considers abundance to be drawn from a negative binomial distribution. The model simultaneously estimates abundance while estimating detection rate from variation in the count sizes, which are assumed to be drawn from a static true population size. The model can also take a vector of NA data, which will tell the model to treat outcomes as NAs when we didn't conduct a survey. Use "./binmix/run_model/run_model.R" to fit this model to the data.

A GLMM can be fit to the mark-recapture devolved count data - "./GLM/models/GLM_negbin_NAs.stan". The model considers abundance to be drawn from a negative binomial distribution. The model does not estimate detection rate. The model can also take a vector of NA data, which will tell the model to treat outcomes as NAs when we didn't conduct a survey. Use "./GLM/run_model/run_model.R" to fit this model to the data.

##--------------------------------------------------------------------------
## HOW TO RUN AN ANLALYSIS ON SIMULATED DATA

We also prepared simulation functions to create communities of species where the underlying properties governing their population sizes are known. We can then fit models to these data to see how accurately and precisely they estimate the known quantities of interest (e.g., the effect size of a hypothetical driver of abundance).

Open "./simulation_full/run_simulations.R" to start simulating some communities. This file will call a simulation function held in "simulate_function.R". Use the source() function to make sure R knows where to find this helper file.

Choose a model type that you will use to analyze the data (GLMM, binmix or multimix). Then select some simulation scenarios. We varied the effect of habitat on detection while simultaneously varying either baseline detection rate or the number of study sites and then looked at the accuracy and precision of each of the models in response.

Open "./simulation_full/run_simulations_closure_violation.R" to start simulating some communities. This file will call a simulation function held in "simulate_function_closure_violations.R". Use the source() function to make sure R knows where to find this helper file.

Choose a model type that you will use to analyze the data (binmix or multimix). Then select some simulation scenarios. We varied the effect of tendency of individuals to remain in a plot (theta0) while simultaneously varying how a hypothetical driver of abundance also impacts the tendency of individuals to remain in a study plot (theta1) then looked at the accuracy and precision of each of the models in response.

##--------------------------------------------------------------------------
## VISUALIZING THE RESULTS

open the "./make_figures/" folder to start plotting results from the real data analysses or simulation analyses.