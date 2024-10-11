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


@knext.node(name="2-reagents - Step #2", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_2reagents)
@knext.input_port(name="Chosen compounds from previous step", description="Input of chosen compounds from the previous step, using the node \"Chosen compounds from previous step\".", port_type=port_type_SATELLiTES)
@knext.input_table(name="Reagents A", description="Input of reagents A.")
@knext.input_table(name="Reagents B", description="Input of reagents B.")
@knext.output_table(name="Enumerated AB compounds", description="Table of final and fully-enumerated AB compounds.")

class Enumeration_2reagents_step2:
    """ **SATELLiTES - 2-reagents - Step #2 - **
    This node produces the enumeration of final and fully-enumerated AB compounds by the combination of reagents A/B with the chosen Ab/aB compounds from the previous step. 

    As inputs tables, you need to provide:
    
    - The **chosen Ab/aB compounds from the previous step** using the squared input port, obtained from the \"Chosen compounds from the previous step\" node.
    - The **reagents A table**, with the Smiles and ID.
    - The **reagents B table**, with the Smiles and ID.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds AB** will be generated. These are the final and fully-assembled compounds. The selection of the best final compounds is determined by the method(s) of your choice.
    """

    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")

    def configure(self, configure_context, input_0, input_1, input_2):
        pass

    def execute(self, execute_context, input_0, input_1, input_2):
        input_0_pandas = input_0._model
        input_1_pandas = input_1.to_pandas()
        input_2_pandas = input_2.to_pandas()

        #############################~

        workdir = execute_context.get_workflow_data_area_dir()
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

        workdir_temp = workdir + "/tmp"
        os.makedirs(workdir_temp, exist_ok=True)

        #############################~

        input_0_temp_Ab = workdir_temp + "/input_0_temp_Ab1.smi"
        input_0_pandas_Ab = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("Ab")]
        input_0_pandas_Ab_IDs = input_0_pandas_Ab.iloc[:,1].str.split("-", expand=True).iloc[:,1].to_list()
        with open(input_0_temp_Ab, "w") as file_output_Ab:
            for input_0_pandas_Ab_ID in input_0_pandas_Ab_IDs:
                file_output_Ab.write("{0}\t{1}\n".format(input_1_pandas.loc[input_1_pandas.iloc[:,1] == input_0_pandas_Ab_ID].iloc[0,0], input_0_pandas_Ab_ID))

        input_0_temp_aB = workdir_temp + "/input_0_temp_aB2.smi"
        input_0_pandas_aB = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("aB")]
        input_0_pandas_aB_IDs = input_0_pandas_aB.iloc[:,1].str.split("-", expand=True).iloc[:,2].to_list()
        with open(input_0_temp_aB, "w") as file_output_aB:
            for input_0_pandas_aB_ID in input_0_pandas_aB_IDs:
                file_output_aB.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_aB_ID].iloc[0,0], input_0_pandas_aB_ID))

        #############################~

        input_1_temp = workdir_temp + "/input_1_temp.smi"
        input_1_pandas.to_csv(input_1_temp, sep = "\t", header = None, index=None)

        input_2_temp = workdir_temp + "/input_2_temp.smi"
        input_2_pandas.to_csv(input_2_temp, sep = "\t", header = None, index=None)

        #############################~

        stenum.run_SATELLiTES_2reagents(input_1_temp, input_0_temp_aB, self.reaction_param, "temp_AB1_step2.smi", workdir_temp, workdir_temp, "AB")
        stenum.run_SATELLiTES_2reagents(input_0_temp_Ab, input_2_temp, self.reaction_param, "temp_AB2_step2.smi", workdir_temp, workdir_temp, "AB")

        #############################~

        output_Ab = pd.read_csv(workdir_temp + "/temp_AB1_step2.smi", sep = "\t", names=["Smiles", "ID"])
        output_aB = pd.read_csv(workdir_temp + "/temp_AB2_step2.smi", sep = "\t", names=["Smiles", "ID"])
        output_cpds = pd.concat([output_Ab, output_aB]).reset_index(drop = True)

        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################