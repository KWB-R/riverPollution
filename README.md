# riverPollution

Functions for quantifying pollutant decay in rivers based on data obtained through Lagrangian sampling in reaches with high potential for interference from ungauged lateral inputs (e.g., discharges from combined sewer overflows in an urban river). The functions are based on a sampling design where several composite samples are taken at >=2 fixed stations during passage of a pollution plume.

Each function is contained in a separate script:

- probInterf.R: Uses (externally) simulated discharges at several points along flow path and, based on a user-given error, computes kernel-density-based probabilities of interference for all samples at all lateral inputs. 

- lagrComp.R: Performs Lagrangian (travel-time based) comparisons of samples at fixed stations along river, computes 1st order rate (k, k>0 indicates concentration increase and vice versa) and writes a table containing the resulting k values and the interference probability of the underlying composite samples.

example.R shows an example based on the following text files (also in the repository):
- latQ.txt: Simulated discharges at lateral input points
- latQx.txt: List of lateral input points and their downstream coordinate (column dx0) and reach. A reach is a river segment between two consecutive sampling stations
- reachInfo.txt: Describes the reaches, their lengths (reachLength), flow velocity (v), travel time (dt). Reaches are numbered in the downstream direction (1, 2, 3....) 
- sampleTimes.txt: Contains the starting (tBeg) and end times (tEnd) of the composite samples at all sampling stations (which are numbered in the downstream direction, such as in latQx for lateral inputs), sample IDs and the downstream coordinate of the sampling stations (dx0)
- sampleConcentrations.txt: Contains each sample's ID, sampling station number, the name of the water quality parameter (variable) and its concentration value

