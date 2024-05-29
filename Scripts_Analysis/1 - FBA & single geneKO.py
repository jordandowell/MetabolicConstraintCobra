#### Python script

#imports
import cobra
import micom
import os
#specific imports?

from cobra.io import load_model



#import model from desktop 
#specific imports
from pathlib import Path
from cobra.io import read_sbml_model
#define data path
data_dir = Path("Data")
data_dir.resolve()
#import and save model
AT_model_path = data_dir / "ATH_Masscorrected_01.xml"
AT_model = read_sbml_model(str(AT_model_path.resolve()))
AT_model_Universal = read_sbml_model(str(AT_model_path.resolve()))
#assess basic descriptions

#check the compartment names if they are missing add them here
#{'c0': 'Cytosol', 'd0': 'Plastid', 'g0': 'Golgi', 'v0': 'Vacuole', 'w0': 'Wall', 'x0': 'Peroxisome', 'm0': 'Mitochondria', 'n0': 'Nucleus', 'r0': 'EndoplasmicReticulum', 'e0': 'External', 'j0': 'Mitochondrialintermembrane'}



#number of reactiuons
print(len(AT_model.reactions))
#number of Metbaolites
print(len(AT_model.metabolites))
#number of genes 

print(len(AT_model.genes))

print(AT_model_Universal.reactions.get_by_id("bio2_biomass"))

print(AT_model.objective)
##### cut an paste pathways here. 
AT_model.compartments |= {'c0': 'Cytosol', 'd0': 'Plastid', 'g0': 'Golgi', 'v0': 'Vacuole', 'w0': 'Wall', 'x0': 'Peroxisome', 'm0': 'Mitochondria', 'n0': 'Nucleus', 'r0': 'EndoplasmicReticulum', 'e0': 'External', 'j0': 'Mitochondrialintermembrane'}
#AT_model.compartments

## adding the camalexin pathway exchange reationcs
#you need exchange reactions because the metabolite needs to be 'consumed'

from cobra import Model, Reaction, Metabolite


print((AT_model.genes))

AT_model_Universal.reactions.get_by_id("Camalexin001_c0")


##### cut an paste pathways here. 



## adding the camalexin pathway 
#makesure you have came from running script 1 

from cobra import Model, Reaction, Metabolite

#Reaction 1
reaction01 = Reaction('Camalexin001_c0')
reaction01.name = 'Tryptophan N-monooxygenase'
reaction01.subsystem = 'c0'
reaction01.lower_bound = 0.  # This is the default
reaction01.upper_bound = 1000.  # This is the default

print(AT_model_Universal.reactions.get_by_id("Camalexin008_c0"))

#Reaction 1 Metabolites * Check if metabolite exists first!

print("exchanges", AT_model.exchanges)
print("demands", AT_model.demands)
print("sinks", AT_model.sinks)

#replace compound with any that you have added to check if you need to add a type of reaction

print("exchanges", AT_model.compartments)
print("demands", AT_model.demands)
print("sinks", AT_model.sinks)


from cobra.flux_analysis import gapfill
universal = cobra.Model("universal_reactions")
with AT_model:
    AT_model.objective = AT_model.add_boundary(AT_model.metabolites.cpd28220_c0, type='demand')
    solution = gapfill(AT_model, universal)
    for reaction in solution[0]:
        print(reaction.id)


#trying to add reactions 

# create demand reactions for HCN and camalexin so can be  'excreted'
AT_model.add_boundary(AT_model.metabolites.get_by_id("HCN_c0"), type="demand")
AT_model.add_boundary(AT_model.metabolites.get_by_id("cpd28220_c0"), type="demand")

#change biomass?
print(AT_model.reactions.bio1_biomass)

NGpermm = -6.2
CamalexinValue =  -NGpermm*10**-6

#add camalexin
AT_model.reactions.bio1_biomass.add_metabolites({"cpd28220_c0": CamalexinValue})
#Remove current Camalexin Value 
CurrentCamalexin = AT_model.reactions.bio1_biomass.metabolites.popitem()
CurrentCamalexin[-1]
AT_model.reactions.bio1_biomass.subtract_metabolites({"cpd28220_c0": CurrentCamalexin[-1]})



AT_model.reactions.bio1_biomass.metabolites = Biomget()AT_model.reactions.bio1_biomass.metabolites = BiomassMetabolites
AT_model.reactions.bio1_biomass.reaction

AT_model.reactions.bio1_biomass.subtract_metabolites({"cpd28220_c0": CamalexinValue})

BiomassMetabolites["cpd28220_c0"] = 0.0



# create exchange reaction
AT_model.add_boundary(AT_model.metabolites.get_by_id("HCN_c0"), type="sink")
AT_model.add_boundary(AT_model.metabolites.get_by_id("cpd28220_c0"), type="sink")

#export genelist
import pandas as pd

df = pd.DataFrame(AT_model.genes)
#we will trim the expression data using this list
df.to_csv('Data/Model_genelist.csv', index=False,header=False) 

#outside check if any need to be translated to the correct name and subset 

#import expression list filtered 


#set objective function
from cobra.util.solver import linear_reaction_coefficients
linear_reaction_coefficients(AT_model)

print(AT_model.summary())

#reducing avliable carbon and oxygen media
#name media
Media = AT_model.medium
#remove oxygen
Media["EX_cpd00007_e0"] = 0.0
#remove sucrose
Media["EX_cpd00076_e0"] = 0.0
#remove Urea
Media["EX_cpd00073_e0"] = 0.0
#add new media amounts back!
AT_model.medium = Media


# change the objective to reaction of choice
AT_model.objective = "Camalexin008_c0"
AT_model.objective = "bio1_biomass"
#bio1 plant leaf 
#bio2


print(AT_model.boundary)
#create new biomass function
print(AT_model.objective)

EX_cpd00076_e0



# change the objective to reaction of choice
AT_model.objective = "Camalexin008_c0"

#run basic flux balance analysis to see if model can grow

AT_solution = AT_model.optimize()
print(AT_model.summary())
print(AT_solution)


#check objective value

AT_solution.objective_value


#adding a new pathway










#check summary


#revisit changing objective functions


#run flux variability analysis 

#imports 

from cobra.flux_analysis import flux_variability_analysis

#FVA (or flux variability analysis) finds the ranges of each metabolic flux at the optimum
#fraction_of_optimium=0.90  
#gives the flux ranges for reactions at 90% optimality.
cobra.flux_analysis.flux_variability_analysis(AT_model, AT_model.reactions[:10], fraction_of_optimum=0.9)

#fun the FVA but looples
 cobra.flux_analysis.geometric_fba(AT_model)


print(AT_solution)
print(AT_model_Universalsolution.fluxes)
#FVA in summar7
print(AT_model.summary(fva=0.95))




from pathlib import Path
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import logging


###########Performing single gene knockouts

#imports
import pandas
from time import time
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

#perform all single gene deletions
deletion_results = single_gene_deletion(AT_model)

print(AT_model.gene_reaction_rule)
print(deletion_results)

#pairwise deletion unecesary tho
double_gene_deletion(
    AT_model, AT_model.genes[-5:]).round(4)
    
    #check medium 
    AT_model.compartments
AT_model.sinks
    
cobra.medium.find_external_compartment(AT_model)
