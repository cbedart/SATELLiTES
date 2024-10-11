##############################################################################################################################################
# SATELLiTES - KNIME extension
# Copyright (C) Corentin BEDART, 2024
##############################################################################################################################################

import knime.extension as knext
from satellites_extension import *

import os
import shutil
import pandas as pd
import pickle

import SATELLiTES.SATELLiTES_enumeration as stenum
import SATELLiTES.SATELLiTES_gui as stgui

##############################################################################################################################################

@knext.node(name="3-reagents - Step #2", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_3reagents)
@knext.input_port(name="Chosen compounds from previous step", description="Input of chosen compounds from the previous step, using the node \"Chosen compounds from previous step\".", port_type=port_type_SATELLiTES)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents C", description="Input of reagents C with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated ABc/AbC/aBC compounds", description="Table of enumerated ABc/AbC/aBC compounds.")

class Enumeration_3reagents_step2:
    """ **SATELLiTES - 3-reagents - Step #2 - **
    This node produces the enumeration of compounds ABc/AbC/aBC by the combination of reagents A/B/C with the defined representatives reagents a/b/c **and** the chosen Abc/aBc/abC compounds from the previous step. 

    As inputs tables, you need to provide:

    - The **chosen Abc/aBc/abC compounds from the previous step** using the squared input port, obtained from the \"Chosen compounds from the previous step\" node.
    - The **reagents A table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents B table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents C table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds ABc, AbC, and aBC** will be generated.

    The selection of these best compounds is determined by the method(s) of your choice. The next step will include the \"3-reagents - Step #3\" node, as well as the "Chosen compounds from previous step" node.
    """
    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")

    def configure(self, configure_context, input_0, input_1, input_2, input_3):
        pass

    def execute(self, execute_context, input_0, input_1, input_2, input_3):
        input_0_pandas = input_0._model
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

        input_0_temp_A = workdir_temp + "/input_0_temp_A.smi"
        input_0_pandas_A = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("Abc")]
        input_0_pandas_A_IDs = input_0_pandas_A.iloc[:,1].str.split("-", expand=True).iloc[:,1].to_list()
        with open(input_0_temp_A, "w") as file_output_A:
            for input_0_pandas_A_ID in input_0_pandas_A_IDs:
                file_output_A.write("{0}\t{1}\n".format(input_1_pandas.loc[input_1_pandas.iloc[:,1] == input_0_pandas_A_ID].iloc[0,0], input_0_pandas_A_ID))

        input_0_temp_B = workdir_temp + "/input_0_temp_B.smi"
        input_0_pandas_B = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("aBc")]
        input_0_pandas_B_IDs = input_0_pandas_B.iloc[:,1].str.split("-", expand=True).iloc[:,2].to_list()
        with open(input_0_temp_B, "w") as file_output_B:
            for input_0_pandas_B_ID in input_0_pandas_B_IDs:
                file_output_B.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_B_ID].iloc[0,0], input_0_pandas_B_ID))

        input_0_temp_C = workdir_temp + "/input_0_temp_C.smi"
        input_0_pandas_C = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("abC")]
        input_0_pandas_C_IDs = input_0_pandas_C.iloc[:,1].str.split("-", expand=True).iloc[:,3].to_list()
        with open(input_0_temp_C, "w") as file_output_C:
            for input_0_pandas_C_ID in input_0_pandas_C_IDs:
                file_output_C.write("{0}\t{1}\n".format(input_3_pandas.loc[input_3_pandas.iloc[:,1] == input_0_pandas_C_ID].iloc[0,0], input_0_pandas_C_ID))


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

        # ABc
        stenum.run_SATELLiTES_3reagents(input_1_temp, input_0_temp_B, representative_3_temp, self.reaction_param, "temp_ABc1.smi", workdir_temp, workdir_temp, "ABc")
        stenum.run_SATELLiTES_3reagents(input_0_temp_A, input_2_temp, representative_3_temp, self.reaction_param, "temp_ABc2.smi", workdir_temp, workdir_temp, "ABc")

        # AbC
        stenum.run_SATELLiTES_3reagents(input_1_temp, representative_2_temp, input_0_temp_C, self.reaction_param, "temp_AbC3.smi", workdir_temp, workdir_temp, "AbC")
        stenum.run_SATELLiTES_3reagents(input_0_temp_A, representative_2_temp, input_3_temp, self.reaction_param, "temp_AbC4.smi", workdir_temp, workdir_temp, "AbC")

        # aBC
        stenum.run_SATELLiTES_3reagents(representative_1_temp, input_0_temp_B, input_3_temp, self.reaction_param, "temp_aBC5.smi", workdir_temp, workdir_temp, "aBC")
        stenum.run_SATELLiTES_3reagents(representative_1_temp, input_2_temp, input_0_temp_C, self.reaction_param, "temp_aBC6.smi", workdir_temp, workdir_temp, "aBC")
        
        #############################~

        output_ABc1 = pd.read_csv(workdir_temp + "/temp_ABc1.smi", sep = "\t", names=["Smiles", "ID"])
        output_ABc2 = pd.read_csv(workdir_temp + "/temp_ABc2.smi", sep = "\t", names=["Smiles", "ID"])
        output_AbC1 = pd.read_csv(workdir_temp + "/temp_AbC3.smi", sep = "\t", names=["Smiles", "ID"])
        output_AbC2 = pd.read_csv(workdir_temp + "/temp_AbC4.smi", sep = "\t", names=["Smiles", "ID"])
        output_aBC1 = pd.read_csv(workdir_temp + "/temp_aBC5.smi", sep = "\t", names=["Smiles", "ID"])
        output_aBC2 = pd.read_csv(workdir_temp + "/temp_aBC6.smi", sep = "\t", names=["Smiles", "ID"])
        output_cpds = pd.concat([output_ABc1, output_ABc2, output_AbC1, output_AbC2, output_aBC1, output_aBC2]).reset_index(drop = True)

        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)
    
##############################################################################################################################################