#####################
Article: Exposure to short dry phases prevents hatching of mayfly egg 
Authors: Emmanuel Jaulin Gautronneau, Abdelkader Azougui, Antoine Callu, Hervé Capra, Maxime Logez, Hervé Pella, Maria Alp

Correspondence author: Emmanuel Jaulin Gautronneau
Email addresses: jaulinemmanuel@orange.fr


#####################
Content of DATA

1. exposure_dry_phase.txt
	Text file with the characteristics of each Baetis egg mass exposed to one dry phase duration. 
	1 row = 1 egg mass.
	
	  Column    |                             Meaning                           |
	
	hatch       |   Number of hatched eggs per mass                             |
	unhatch     |   Number of fertilised unhatched eggs per mass                |
	totegg      |   Number of fertilised eggs per mass                          |
	prop        |	Hatching proportion per mass                                |
	flume       |	ID of the flume in which the egg mass was located           |
	sub         |	ID of the substrate on which the egg mass was located       |
	dry_time    |   Dry phase duration (h) for which the egg mass was exposed   |

2. R_script.R
	R file with the script to:
				  - Calculate the hatching proportion of egg masses for dry phase durations and/or slates detailed in 				    the article
				  - Run the generalized linear mixed model

				  - Display Figures 6 and S5
