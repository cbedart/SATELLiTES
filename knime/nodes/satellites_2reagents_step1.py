##############################################################################################################################################
# SATELLiTES - KNIME extension
# Copyright (C) Corentin BEDART, 2024
##############################################################################################################################################

import knime.extension as knext
from satellites_extension import *

import os
import shutil
import pandas as pd

import SATELLiTES.SATELLiTES_enumeration as stenum
import SATELLiTES.SATELLiTES_gui as stgui

##############################################################################################################################################

@knext.node(name="2-reagents - Step #1", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_2reagents)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated Ab/aB compounds", description="Table of enumerated Ab/aB compounds.")

class Enumeration_2reagents_step1:
    """ **SATELLiTES - 2-reagents - Step #1 - **
    This node produces the enumeration of compounds Ab/aB by the combination of reagents A/B with the defined representatives reagents a/b. 
    
    As inputs tables, you need to provide:

    - The **reagents A table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents B table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.


    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds Ab and aB** will be generated.

    The selection of these best compounds is determined by the method(s) of your choice. The next step will include the \"2-reagents - Step #2\" node, as well as the "Chosen compounds from previous step" node.
    """

    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")

    def configure(self, configure_context, input_1, input_2):
        pass

    def execute(self, execute_context, input_1, input_2):
        input_1_pandas = input_1.to_pandas()
        input_2_pandas = input_2.to_pandas()

        #############################~

        workdir = execute_context.get_workflow_data_area_dir()
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

        workdir_temp = workdir + "/tmp"
        os.makedirs(workdir_temp, exist_ok=True)

        #############################~

        representative_checker = [0,0]

        if "Representative" in input_1_pandas.columns:
            representative_1_pandas = input_1_pandas.loc[input_1_pandas["Representative"] == True]
            representative_1_temp = workdir_temp + "/representative_1_temp.smi"

            if representative_1_pandas["ID"].to_list()[0] == "Custom":
                representative_1_pandas["ID"] = "CustomA"
                input_1_pandas = input_1_pandas.drop(input_1_pandas.loc[input_1_pandas["Representative"] == True].index)
                
            representative_1_pandas.to_csv(representative_1_temp, sep = "\t", header = None, index=None)
            representative_checker[0] = 1

        if "Representative" in input_2_pandas.columns:
            representative_2_pandas = input_2_pandas.loc[input_2_pandas["Representative"] == True]
            representative_2_temp = workdir_temp + "/representative_2_temp.smi"

            if representative_2_pandas["ID"].to_list()[0] == "Custom":
                representative_2_pandas["ID"] = "CustomB"
                input_2_pandas = input_2_pandas.drop(input_2_pandas.loc[input_2_pandas["Representative"] == True].index)

            representative_2_pandas.to_csv(representative_2_temp, sep = "\t", header = None, index=None)
            representative_checker[1] = 1

        #############################~

        if representative_checker != [1,1]:
            TypeError("No defined representative.")

        #############################~

        input_1_temp = workdir_temp + "/input_1_temp.smi"
        input_1_pandas.to_csv(input_1_temp, sep = "\t", header = None, index=None)

        input_2_temp = workdir_temp + "/input_2_temp.smi"
        input_2_pandas.to_csv(input_2_temp, sep = "\t", header = None, index=None)

        #############################~

        stenum.run_SATELLiTES_2reagents(input_1_temp, representative_2_temp, self.reaction_param, "temp_Ab1.smi", workdir_temp, workdir_temp, "Ab")
        stenum.run_SATELLiTES_2reagents(representative_1_temp, input_2_temp, self.reaction_param, "temp_aB2.smi", workdir_temp, workdir_temp, "aB")
        
        #############################~

        output_Ab = pd.read_csv(workdir_temp + "/temp_Ab1.smi", sep = "\t", names=["Smiles", "ID"])
        output_aB = pd.read_csv(workdir_temp + "/temp_aB2.smi", sep = "\t", names=["Smiles", "ID"])
        output_cpds = pd.concat([output_Ab, output_aB]).reset_index(drop = True)

        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~
        
        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################