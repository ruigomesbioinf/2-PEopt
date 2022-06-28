from cobra.io import read_sbml_model
from mewpy.optimization.evaluation import BPCY, TargetFlux, WYIELD
from mewpy.optimization import EA
from mewpy.problems import GOUProblem
from mewpy.util.io import population_to_csv

import warnings
warnings.filterwarnings("ignore")

CPUS = 32
ITERATIONS = 20

def GOU_problem(filename):
    """
    Function that loads a configuration for the optimization of 2-phenylethanol on Escherichia coli model iML1515 edited
    by adding Ehrlich pathway and defines and run a gene over/under expression optimization problem.

    Args:
        filename (str): Name of the file returned by this function.
    """
    
    # load model
    model=read_sbml_model("models/edited_ECOLI.xml")

    # set default solver as cplex
    model.solver="cplex"

    BIOMASS_ID="BIOMASS_Ec_iML1515_core_75p37M"
    PRODUCT_ID="EX_2phetoh_e"

    O2="EX_o2_e"
    GLC = "EX_glc__D_e"

    # environmental conditions
    envcond = envcond = {
        O2: (-20.0, 100000.0), 
        GLC: (-10.0, 100000.0)
        }

    # optimization objectives
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method = "pFBA")
    evaluator_2 = TargetFlux(BIOMASS_ID, method = "ROOM")
    evaluator_3 = WYIELD(BIOMASS_ID, PRODUCT_ID, method = "ROOM")

    # GOUproblem (Gene Over and Under expression problem)
    gouproblem = GOUProblem(model, [evaluator_1, evaluator_2, evaluator_3], envcond = envcond, candidate_max_size = 10)

    # run and save output to csv file
    ea = EA(gouproblem, max_generations = ITERATIONS, algorithm = "NSGAIII")
    final_population = ea.run()
    population_to_csv(gouproblem, final_population, filename, simplify = True)

if __name__ == "__main__":
    from mewpy.util.constants import EAConstants
    EAConstants.NUM_CPUS = CPUS
    for i in range(1):
        filename = "output_ECOLI_GOU.csv"
        GOU_problem(filename)