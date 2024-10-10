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

@knext.node(name="3-reagents - Step #3", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_3reagents)
@knext.input_port(name="Chosen compounds from previous step", description="Input of chosen compounds from the previous step, using the node \"Chosen compounds from previous step\".", port_type=port_type_SATELLiTES)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents C", description="Input of reagents C with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated ABC compounds", description="Table of final and fully-enumerated ABC compounds.")

class Enumeration_3reagents_step3:
    """ **SATELLiTES - 3-reagents - Step #3 - **
    This node produces the enumeration of final and fully-enumerated ABC by the combination of reagents A/B/C with the chosen ABc/AbC/aBC compounds from the previous step. 

    As inputs tables, you need to provide:

    - The **chosen ABc/AbC/aBC compounds from the previous step** using the squared input port, obtained from the \"Chosen compounds from the previous step\" node.
    - The **reagents A table**, with the Smiles and ID.
    - The **reagents B table**, with the Smiles and ID.
    - The **reagents C table**, with the Smiles and ID.

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

        input_0_temp_AB_A = workdir_temp + "/input_0_temp_AB_A.smi"
        input_0_pandas_AB_A = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("ABc")]
        input_0_pandas_AB_A_IDs = input_0_pandas_AB_A.iloc[:,1].str.split("-", expand=True).iloc[:,1].to_list()
        with open(input_0_temp_AB_A, "w") as file_output_AB_A:
            for input_0_pandas_AB_A_ID in input_0_pandas_AB_A_IDs:
                file_output_AB_A.write("{0}\t{1}\n".format(input_1_pandas.loc[input_1_pandas.iloc[:,1] == input_0_pandas_AB_A_ID].iloc[0,0], input_0_pandas_AB_A_ID))

        input_0_temp_AB_B = workdir_temp + "/input_0_temp_AB_B.smi"
        input_0_pandas_AB_B = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("ABc")]
        input_0_pandas_AB_B_IDs = input_0_pandas_AB_B.iloc[:,1].str.split("-", expand=True).iloc[:,2].to_list()
        with open(input_0_temp_AB_B, "w") as file_output_AB_B:
            for input_0_pandas_AB_B_ID in input_0_pandas_AB_B_IDs:
                file_output_AB_B.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_AB_B_ID].iloc[0,0], input_0_pandas_AB_B_ID))

        #####~

        input_0_temp_AC_A = workdir_temp + "/input_0_temp_AC_A.smi"
        input_0_pandas_AC_A = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("AbC")]
        input_0_pandas_AC_A_IDs = input_0_pandas_AC_A.iloc[:,1].str.split("-", expand=True).iloc[:,1].to_list()
        with open(input_0_temp_AC_A, "w") as file_output_AC_A:
            for input_0_pandas_AC_A_ID in input_0_pandas_AC_A_IDs:
                file_output_AC_A.write("{0}\t{1}\n".format(input_1_pandas.loc[input_1_pandas.iloc[:,1] == input_0_pandas_AC_A_ID].iloc[0,0], input_0_pandas_AC_A_ID))

        input_0_temp_AC_C = workdir_temp + "/input_0_temp_AC_C.smi"
        input_0_pandas_AC_C = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("AbC")]
        input_0_pandas_AC_C_IDs = input_0_pandas_AC_C.iloc[:,1].str.split("-", expand=True).iloc[:,3].to_list()
        with open(input_0_temp_AC_C, "w") as file_output_AC_C:
            for input_0_pandas_AC_C_ID in input_0_pandas_AC_C_IDs:
                file_output_AC_C.write("{0}\t{1}\n".format(input_3_pandas.loc[input_3_pandas.iloc[:,1] == input_0_pandas_AC_C_ID].iloc[0,0], input_0_pandas_AC_C_ID))

        #####~

        input_0_temp_BC_B = workdir_temp + "/input_0_temp_BC_B.smi"
        input_0_pandas_BC_B = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("aBC")]
        input_0_pandas_BC_B_IDs = input_0_pandas_BC_B.iloc[:,1].str.split("-", expand=True).iloc[:,2].to_list()
        with open(input_0_temp_BC_B, "w") as file_output_BC_B:
            for input_0_pandas_BC_B_ID in input_0_pandas_BC_B_IDs:
                file_output_BC_B.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_BC_B_ID].iloc[0,0], input_0_pandas_BC_B_ID))

        input_0_temp_BC_C = workdir_temp + "/input_0_temp_BC_C.smi"
        input_0_pandas_BC_C = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("aBC")]
        input_0_pandas_BC_C_IDs = input_0_pandas_BC_C.iloc[:,1].str.split("-", expand=True).iloc[:,3].to_list()
        with open(input_0_temp_BC_C, "w") as file_output_BC_C:
            for input_0_pandas_BC_C_ID in input_0_pandas_BC_C_IDs:
                file_output_BC_C.write("{0}\t{1}\n".format(input_3_pandas.loc[input_3_pandas.iloc[:,1] == input_0_pandas_BC_C_ID].iloc[0,0], input_0_pandas_BC_C_ID))

        #############################~

        input_1_temp = workdir_temp + "/input_1_temp.smi"
        input_1_pandas.to_csv(input_1_temp, sep = "\t", header = None, index=None)

        input_2_temp = workdir_temp + "/input_2_temp.smi"
        input_2_pandas.to_csv(input_2_temp, sep = "\t", header = None, index=None)

        input_3_temp = workdir_temp + "/input_3_temp.smi"
        input_3_pandas.to_csv(input_3_temp, sep = "\t", header = None, index=None)

        #############################~

        stenum.run_SATELLiTES_3reagents(input_0_temp_AB_A, input_0_temp_AB_B, input_3_temp, self.reaction_param, "temp_ABC1.smi", workdir_temp, workdir_temp, "ABC")
        stenum.run_SATELLiTES_3reagents(input_0_temp_AC_A, input_2_temp, input_0_temp_AC_C, self.reaction_param, "temp_ABC2.smi", workdir_temp, workdir_temp, "ABC")
        stenum.run_SATELLiTES_3reagents(input_1_temp, input_0_temp_BC_B, input_0_temp_BC_C, self.reaction_param, "temp_ABC3.smi", workdir_temp, workdir_temp, "ABC")

        #############################~

        output_ABC1 = pd.read_csv(workdir_temp + "/temp_ABC1.smi", sep = "\t", names=["Smiles", "ID"])
        output_ABC2 = pd.read_csv(workdir_temp + "/temp_ABC2.smi", sep = "\t", names=["Smiles", "ID"])
        output_ABC3 = pd.read_csv(workdir_temp + "/temp_ABC3.smi", sep = "\t", names=["Smiles", "ID"])
        output_cpds = pd.concat([output_ABC1, output_ABC2, output_ABC3]).reset_index(drop = True)

        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)
    
##############################################################################################################################################