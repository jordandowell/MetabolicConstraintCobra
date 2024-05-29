##### Using RNA for Biosyntheic constraint

#Make sure in python 2.7
#imports
import cobra
import micom
import os
#specific imports

from cobra.io import load_model
import math
#import model from desktop 
#specific imports
from pathlib import Path
from cobra.io import read_sbml_model

#importmodel if necessary 
AT_model = read_sbml_model(str(AT_model_path.resolve()))
#Make sure that it will run with 0 constraints

##### Add compartments. 
AT_model.compartments |= {'c0': 'Cytosol', 'd0': 'Plastid', 'g0': 'Golgi', 'v0': 'Vacuole', 'w0': 'Wall', 'x0': 'Peroxisome', 'm0': 'Mitochondria', 'n0': 'Nucleus', 'r0': 'EndoplasmicReticulum', 'e0': 'External', 'j0': 'Mitochondrialintermembrane'}

#add camalexin
NGpermm = -6.2 #camalexin EC50 for botrytis
CamalexinValue =  NGpermm*10**-6
AT_model.reactions.bio1_biomass.add_metabolites({"cpd28220_c0": CamalexinValue})
AT_model.reactions.bio1_biomass.reaction
# create exchange reactions for HCN and camalexin so they are 'excreted'
AT_model.add_boundary(AT_model.metabolites.get_by_id("HCN_c0"), type="sink")

#reducing avliable carbon sources and oxygen media
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
from cobra.medium import minimal_medium
Minimalsolution=minimal_medium(AT_model, 1)
print(Minimalsolution)
minimal_medium(AT_model, AT_model.slim_optimize())
AT_model.slim_optimize()
AT_model.optimize()

####import Bioconstraint
BiosyntheticConstraints = pd.read_csv("Data/BiosyntheticConstraints.csv")
BiosyntheticConstraints.head()


#import Camalexin Constraints
CamalexinConstraints = pd.read_csv("Data/Camalexin_Constraints.csv")
CamalexinConstraints.head()
#CamalexinValues = CamalexinConstraints.iloc[j,1]

BiosyntheticConstraints.drop(columns=BiosyntheticConstraints.columns[0]).columns.values.tolist()




#test replacement of constraints for one column

  df = pd.DataFrame(AT_solution.fluxes)
#we will trim the expression data using this list
df.to_csv('Data/Model_fluxes.csv', index=True,header=False) 
  
# we need to get model solutions in a single dataframe.
#REMINDER min flux is = to 1
#get the path to the outputs
#define data path
Output_dir = Path("Output")
Output_dir.resolve()
#create a set of minmum flux values to be assessed 
Minimum_FLUX = {0,1,10,100,1000,10000,100000}
os.walk(Output_dir)[1]

Minimum_FLUX_Value = float(1)





FLUX_DATA = []
j = 0
#get iterate of expressions set for each individuals 
for j in range(1,len(BiosyntheticConstraints.columns)):
  # iterate over each individual reactions
  i = 0
  for i in range(0,len(AT_model.reactions)):
  
    #check if reaction has associated gene 
    if bool(AT_model.reactions[i].gene_reaction_rule):
      #get a list of the genes in the reactions 
      Genes = AT_model.reactions[i].gene_name_reaction_rule.split()
      #get the expression for the genes needed
      rows = pd.DataFrame(BiosyntheticConstraints[BiosyntheticConstraints['gene'].isin(Genes)])
  
  
      #thing to repaplace in string
      
      
      #change the string into a useful equation to find the bound to use 
      GENE_String = AT_model.reactions[i].gene_reaction_rule
      # ( = min(
      GENE_String = GENE_String.replace('(','min(1000000000000000, ')
      # and = ,
      GENE_String = GENE_String.replace('and',',')
      # or = +
      GENE_String = GENE_String.replace('or','+')
      GENE_String = GENE_String.replace(')',' )')
      
      # simple add genes
          # ( = min(
     # GENE_String = GENE_String.replace('(','min(1000000000000000, ')
      # and = ,
    #  GENE_String = GENE_String.replace('and',',')
      # or = +
     # GENE_String = GENE_String.replace('or','+')
    #  GENE_String = GENE_String.replace('and','+')
    #  GENE_String = GENE_String.replace(')',' )')
      
      
      
      #replace the avliable genes with expression values
      for index in range(0,rows.shape[0]): 
        GENE_String = GENE_String.replace(rows.iat[index,0],str(rows.iat[index,j]))
        GENE_String_LIST = GENE_String.split()
      #add any 0s for missing genes
      for Emptygene in range(0,len(GENE_String_LIST)):
        if  'A' in GENE_String_LIST[Emptygene]:
         GENE_String = GENE_String.replace(GENE_String_LIST[Emptygene],str(0))
       #extra value is added just incase there is no min
      NewBound =  min(eval('min(1000000000000000,' + GENE_String + ',1000000000000000)'),1000000000000000)
      #ensure that at least some flux can occur +1
      NewBound = float(math.ceil(NewBound+Minimum_FLUX_Value))
      #setting bounds 
     
      AT_model.reactions[i].upper_bound = float(NewBound)
      #if the lower bound is 0 keep it at 0
      if bool(AT_model.reactions[i].reversibility):
        AT_model.reactions[i].lower_bound = float(-(NewBound))
      else:
        AT_model.reactions[i].lower_bound = float(0)
        #if the upper bound is 0 keep it at 0
      #if bool(AT_model.reactions[i].upper_bound == 0):
        #AT_model.reactions[i].upper_bound = 0
      #else: 
        # AT_model.reactions[i].upper_bound = NewBound
        
        
        ##Add Camalexin biomass addition
        #add camalexin
      #find camalexin row to add
     
       
  #find host camalexin row data 
       
  HostCamelxinRow = CamalexinConstraints.loc[CamalexinConstraints['SampleName'] == BiosyntheticConstraints.columns[j]]
  HostCamelxinValue = HostCamelxinRow.iloc[0,2]
  
  #Clear camalexin
  #Remove current Camalexin Value 
  CurrentCamalexin = AT_model.reactions.bio1_biomass.metabolites.popitem()
  AT_model.reactions.bio1_biomass.subtract_metabolites({"cpd28220_c0": CurrentCamalexin[-1]})
  #add new Camalexin value on a ng/g basis 
  AT_model.reactions.bio1_biomass.add_metabolites({"cpd28220_c0": -HostCamelxinValue*10**-6})

    ##SOLVE Model

  AT_solution = AT_model.optimize()
  # store DataFrame in list
  FLUX_DATA.append(pd.DataFrame(AT_solution.fluxes))
  print(AT_solution)



# see pd.concat documentation for more info
FLUX_DATA = pd.concat(FLUX_DATA,axis=1)
#add Column names
FLUX_DATA.columns = BiosyntheticConstraints.drop(columns=BiosyntheticConstraints.columns[0]).columns.values.tolist()

# write DataFrame to an excel sheet 
FLUX_DATA.to_csv('Data/Model_fluxes.csv', index=True,header=True) 






#create Control where all reactions flux is set to min_Flux 
  # iterate over each individual reactions
i = 0
for i in range(0,len(AT_model.reactions)):
  NewBound = float(math.ceil(Minimum_FLUX_Value))
  #setting bounds 
  AT_model.reactions[i].upper_bound = float(NewBound)
  #if the lower bound is 0 keep it at 0
  if bool(AT_model.reactions[i].reversibility):
    AT_model.reactions[i].lower_bound = float(-(NewBound))
  else:
     AT_model.reactions[i].lower_bound = float(0)
      
AT_solution = AT_model.optimize()
  # store DataFrame in list
  FLUX_DATA.append(pd.DataFrame(AT_solution.fluxes))
  print(AT_solution)









# see pd.concat documentation for more info
FLUX_DATA = pd.concat(FLUX_DATA,axis=1)
# write DataFrame to an excel sheet 
FLUX_DATA.to_csv('Data/Model_fluxes.csv', index=True,header=False) 

### need max solution everything set to 1000 * normal model

#### all bounds set to 1
#### all bounds set to 10
#### all bounds set to 100
#### all bounds set to 1000
#### all bounds set to 10000
#### all bounds set to 100000

a = [0, 0, 0, 0, 0]

min(enumerate(a), key=lambda x: x[1] if x[1] > 0 else float(0))
enumerate(a):[1]
min([value for value in a if value > 0 else float(0)])

float(min(a))

#############
 i = 0
   AT_solution = AT_model.optimize()
    print(AT_solution)
     AT_model_Universalsolution = AT_model_Universal.slim_optimize()
    print(AT_model_Universalsolution)

#did it grow?



AT_model_Universal.tolerance
AT_model.problem.Constraint.expression



model.compartments['e']
print("exchanges", AT_model.exchanges)
print("demands", AT_model.demands)
print("sinks", AT_model.sinks)

AT_solution = AT_model.slim_optimize()
print(AT_solution.)
AT_model_Universalsolution = AT_model_Universal.optimize()
print(AT_model_Universalsolution)

#find missing reactions

cobra.flux_analysis.find_blocked_reactions(AT_model)

solution = gapfill(AT_model, AT_model_Universal, demand_reactions=False)
for reaction in solution[0]:
    print(reaction.id)
consistent_model = cobra.flux_analysis.fastcc(AT_model)

#for loop

#for i in 2:length(Bioconstrants)

#replace max constraint by Gene ID

#solve

#export solved model 





#head back to R for GWAS
from cobra.io import load_model
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(AT_model,["EX_cpd11632_e0","EX_cpd00076_e0"])
prod_env.head()



# see pd.concat documentation for more info

print(prod_env)

# write DataFrame to an excel sheet 
prod_env.to_csv('Data/Production_fluxes.csv', index=False,header=True) 
