# zombie_forests
Trailing-edge zombie forests can increase population persistence in the face of climate change.

The code solves a spatial integrodifference equation model of forests shifting in response to climate change, and examines the role of individuals at the trailing range edge in promoting population persistance. The model draws on two example tree species from the Plant Matrix Database (Compadre), but the model is general enough to work for any species in the database, or more generally for any transition matrix describing plant demography.

The pre-print manuscript associated with this project is 'zombie_forests_biorxiv.pdf'.

The analyses of the model are run and plotted in 'zombieforests.R', which sources the R files 'zombieNumerics.R' and 'zombiefunctions.R'. The data 'COMPADRE_v.4.0.1.Rdata' was downloaded from the Plant Matrix Database at https://compadre-db.org/.
