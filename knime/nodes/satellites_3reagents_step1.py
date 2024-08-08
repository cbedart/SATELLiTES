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

@knext.node(name="3-reagents - Step #1", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_3reagents)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents C", description="Input of reagents C with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated Abc/aBc/abC compounds", description="Table of enumerated Abc/aBc/abC compounds.")

class Enumeration_3reagents_step1:
    """ **SATELLiTES - 3-reagents - Step #1 - **
    This node produces the enumeration of compounds Abc/aBc/abC by the combination of reagents A/B/C with the defined representatives reagents a/b/c. 

    As inputs tables, you need to provide:

    - The **reagents A table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents B table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents C table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds Abc, aBc, and abC** will be generated.

    The selection of these best compounds is determined by the method(s) of your choice. The next step will include the \"3-reagents - Step #2\" node, as well as the "Chosen compounds from previous step" node.
    """

    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")

    def configure(self, configure_context, input_1, input_2, input_3):
        pass

    def execute(self, execute_context, input_1, input_2, input_3):
        input_1_pandas = input_1.to_pandas()
        input_2_pandas = input_2.to_pandas()
        input_3_pandas = input_3.to_pandas()

        #############################~

        workdir = execute_context.get_workflow_data_area_dir()
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

        workdir_temp = workdir + "/tmp"
        os.makedirs(workdir_temp, exist_ok=True)

        #############################~

        representative_checker = [0,0,0]

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

        if "Representative" in input_3_pandas.columns:
            representative_3_pandas = input_3_pandas.loc[input_3_pandas["Representative"] == True]
            representative_3_temp = workdir_temp + "/representative_3_temp.smi"

            if representative_3_pandas["ID"].to_list()[0] == "Custom":
                representative_3_pandas["ID"] = "CustomC"
                input_3_pandas = input_3_pandas.drop(input_3_pandas.loc[input_3_pandas["Representative"] == True].index)

            representative_3_pandas.to_csv(representative_3_temp, sep = "\t", header = None, index=None)
            representative_checker[2] = 1

        #############################~

        if representative_checker != [1,1,1]:
            TypeError("No defined representative.")

        #############################~

        input_1_temp = workdir_temp + "/input_1_temp.smi"
        input_1_pandas.to_csv(input_1_temp, sep = "\t", header = None, index=None)

        input_2_temp = workdir_temp + "/input_2_temp.smi"
        input_2_pandas.to_csv(input_2_temp, sep = "\t", header = None, index=None)

        input_3_temp = workdir_temp + "/input_3_temp.smi"
        input_3_pandas.to_csv(input_3_temp, sep = "\t", header = None, index=None)

        #############################~

        stenum.run_SATELLiTES_3reagents(input_1_temp, representative_2_temp, representative_3_temp, self.reaction_param, "temp_Abc.smi", workdir_temp, workdir_temp, "Abc")
        stenum.run_SATELLiTES_3reagents(representative_1_temp, input_2_temp, representative_3_temp, self.reaction_param, "temp_aBc.smi", workdir_temp, workdir_temp, "aBc")
        stenum.run_SATELLiTES_3reagents(representative_1_temp, representative_2_temp, input_3_temp, self.reaction_param, "temp_abC.smi", workdir_temp, workdir_temp, "abC")

        #############################~

        output_Abc = pd.read_csv(workdir_temp + "/temp_Abc.smi", sep = "\t", names=["Smiles", "ID"])
        output_aBc = pd.read_csv(workdir_temp + "/temp_aBc.smi", sep = "\t", names=["Smiles", "ID"])
        output_abC = pd.read_csv(workdir_temp + "/temp_abC.smi", sep = "\t", names=["Smiles", "ID"])
        output_cpds = pd.concat([output_Abc, output_aBc, output_abC]).reset_index(drop = True)

        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################