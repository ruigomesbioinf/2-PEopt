from cobra.io import read_sbml_model
from cobra.io import write_sbml_model
from cobra import Reaction, Metabolite, Model

modelec = read_sbml_model("models/iML1515.xml")
model = Model(modelec)

# CREATE CHECKPOINT TO SEE IF THE MODIFICATIONS WERE SUCESSFULL 
# (initial quantity of reactions, metabolites and genes in our model)
print("INITIAL QUANTITIES")
print(20*"-")
print(f'{len(model.reactions)} reactions initially')
print(f'{len(model.metabolites)} metabolites initially')
print(f'{len(model.genes)} genes initially')
print(20*"-")

# ADD METABOLITE 2-PHENYLETHANOL AS "2phetoh_c"
phenylethanol_c = Metabolite(
    "2phetoh_c",
    formula="C8H10O",
    name="2-phenylethanol",
    compartment="c")

phenylethanol_e = Metabolite(
    "2phetoh_e",
    formula="C8H10O",
    name="2-phenylethanol",
    compartment="e")

model.add_metabolites([phenylethanol_c, phenylethanol_e])


# CREATE EHRLICH PATHWAY REACTIONS INSTANCES (3 REACTIONS)
# phenylalanine transaminase (PHETA1) -> akg_c + phe__L_c ⇌ glu__L_c + phpyr_c (já está no modelo)

# phenylpyruvate decarboxylase (PPYRDC) -> h_c + phpyr_c ⇌ co2_c + pacald_c
ppyrdc_r = Reaction("PPYRDC")
ppyrdc_r.name = "phenylpyruvate decarboxylase (PPYRDC)"

M_h_c = model.metabolites.get_by_id(id="h_c")
M_phpyr_c = model.metabolites.get_by_id(id="phpyr_c")
M_glu__L_c = model.metabolites.get_by_id(id="glu__L_c")
M_pacald_c = model.metabolites.get_by_id(id="pacald_c")

ppyrdc_r.add_metabolites({
    M_h_c: -1.0, 
    M_phpyr_c: -1.0, 
    M_glu__L_c: 1.0, 
    M_pacald_c: 1.0})

# aldehyde dehydrogenase (ALCD25xi) -> h_c + nadh_c + pacald_c ⇌ nad_c + 2phetoh_c
ALCD25xi_r = Reaction("ALCD25xi")
ALCD25xi_r.name = "Aldehyde dehydrogenase 2 phenylethanol NAD"

M_nadh_c = model.metabolites.get_by_id(id="nadh_c")
M_nad_c = model.metabolites.get_by_id(id="nad_c")
M_2phetoh_c = model.metabolites.get_by_id(id="2phetoh_c")

ALCD25xi_r.add_metabolites({
    M_h_c: -1, 
    M_nadh_c: -1, 
    M_pacald_c: -1, 
    M_nad_c: 1, 
    M_2phetoh_c: 1})

# ADD 2PHENYLETHANOL EXCHANGE REACTION (R_EX_2phetoh_c)
R_EX_2phetoh_e = Reaction("EX_2phetoh_e")
R_EX_2phetoh_e.name = "2-phenylethanol exchange"

M_2phetoh_e = model.metabolites.get_by_id(id = "2phetoh_e")
R_EX_2phetoh_e.add_metabolites({
    M_2phetoh_e: -1})

# ADD 2PHENYLETHANOL TRANSPORT REACTION (R_EX_2phetoh_trans)
R_EX_2phetoh_trans = Reaction("EX_2phetoh_trans")
R_EX_2phetoh_trans.name = "2-phenylethanol transport"

R_EX_2phetoh_trans.add_metabolites({
    M_2phetoh_e: -1,
    M_2phetoh_c: 1})

# ADD REACTIONS AND METABOLITES
model.add_reactions([ppyrdc_r, ALCD25xi_r, R_EX_2phetoh_e, R_EX_2phetoh_trans])

# CHECK IF ADITIONS WERE MADE TO THE MODEL
print("FINAL QUANTITIES")
print(20*"-")
print(f'{len(model.reactions)} reactions at the end')
print(f'{len(model.metabolites)} metabolites at the end')
print(f'{len(model.genes)} genes at the end')
print(20*"-")

# SAVE EDITED MODEL
write_sbml_model(model, "edited_ECOLI.xml")

# CHECK IF IT'S A VALIDE SBML FORMAT MODEL
from cobra.io import validate_sbml_model

validate_sbml_model("models/edited_ECOLI.xml")