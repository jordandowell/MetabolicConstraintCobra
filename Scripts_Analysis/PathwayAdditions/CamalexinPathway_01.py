## adding the camalexin pathway 
#makesure you have came from running script 1 

from cobra import Model, Reaction, Metabolite

#Reaction 1
reaction01 = Reaction('Camalexin001_c0')
reaction01.name = 'Tryptophan N-monooxygenase'
reaction01.subsystem = 'c0'
reaction01.lower_bound = 0.  # This is the default
reaction01.upper_bound = 1000.  # This is the default


AT_model.get_by_id("cpd00005")
