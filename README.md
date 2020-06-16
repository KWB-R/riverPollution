# riverPollution

Functions for processing data on pollutant decay in rivers obtained through Lagrangian sampling in reaches with high potential for interference from ungauged lateral discharges (e.g., combined sewer overflows in an urban river).

Each function is contained in a separate script:

- probInterf.R: Uses (externally) simulated discharges at several points along flow path and, based on a user-given error, computes kernel-density-based probabilities of interference for all samples at all lateral inputs. 

- lagrComp.R: Performs Lagrangian (travel-time based) comparisons of samples at fixed stations along river, computes 1st order rate (k, k>0 indicates concentration increase and vice versa) and writes out results table containing the resulting k values.
