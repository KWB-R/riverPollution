# riverPollution

Functions for processing data on pollutant decay in rivers obtained through Lagrangian sampling in reaches with high potential for interference from ungauged lateral discharges (e.g., combined sewer overflows in an urban river). The functions are based on a sampling design where several composite samples are taken at >=2 fixed stations during passage of a pollution plume.

These functions were developed in the context of research project FLUSSHYGIENE: https://www.kompetenz-wasser.de/de/project/flusshygiene/.

Each function is contained in a separate script:

- probInterf.R: Uses (externally) simulated discharges at several points along flow path and, based on a user-given error, computes kernel-density-based probabilities of interference for all samples at all lateral inputs. 

- lagrComp.R: Performs Lagrangian (travel-time based) comparisons of samples at fixed stations along river, computes 1st order rate (k, k>0 indicates concentration increase and vice versa) and writes out results table containing the resulting k values.

example.R shows an example based on the following text files (also in the repository):
- latQ.txt: Simulated discharges at lateral input points
- latQx.txt: List of lateral input points and their downstream coordinate (column dx0) and reach. A reach is a river segment between two consecutive sampling stations
- reachInfo.txt: Describes the reaches, their lengths (reachLength), flow velocity (v), travel time (dt). Reaches are numbered in the downstream direction (1, 2, 3....) 
- sampleTimes.txt: Contains the starting (tBeg) and end times (tEnd) of the composite samples

station;tBeg;tEnd;sample;dx0
