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

@knext.node(name="Representative - Custom", node_type=knext.NodeType.SOURCE, icon_path="icon.png", category=category_utilities)
@knext.input_table(name="Input Data", description="We read data from here")
@knext.output_table(name="Output Data", description="Whatever the node has produced")
class Representative_custom:
    """ **SATELLiTES - Custom representative**
    The **\"Representative\"** set of nodes is used to determine the reference reagent (named \"representative compound\" with its letter written in lower case) from a set of ligands.
    For example, for the reagent subset that will be used as **reagents A**, its representative compound will be specified as **representative a**.

    As an input table, you need to provide a **reagents table**, with the Smiles and ID. The easiest way is to use the \"File Reader\" node to import the .smi/.csv/.tsv/other file into Knime.

    For this specific node, the representative compound will be selected from an **user-specified SMILES** available in the input set. The representative compound will be referred as \"CustomX\" in the following steps.

    The next step will include the \"X-reagents - Step #X\" set of nodes.
    """

    smiles_param = knext.StringParameter("Custom SMILES", "Custom SMILES reagent that will be used as the representative.")

    def configure(self, configure_context, input_1):
        pass

    def execute(self, execute_context, input_1):
        input_1_pandas = input_1.to_pandas()

        #############################~

        workdir = execute_context.get_workflow_data_area_dir()
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

        workdir_temp = workdir + "/tmp"
        os.makedirs(workdir_temp, exist_ok=True)

        #############################~

        input_1_temp = workdir_temp + "/smallest_input.smi"
        input_1_pandas.to_csv(input_1_temp, sep = "\t", header = None, index=None)

        #############################~

        try:
            stenum.Chem.AddHs(stenum.Chem.MolFromSmiles(input_1_pandas.iloc[0,0]))
            output_temp = input_1_pandas.iloc[:,[0,1]]
        except:
            try:
                stenum.Chem.AddHs(stenum.Chem.MolFromSmiles(input_1_pandas.iloc[0,1]))
                output_temp = input_1_pandas.iloc[:,[1,0]]
            except:
                raise ValueError("Issue with input SMILES - The SMILES and their ID must be in the first 2 columns")

        #############################~

        output_temp.columns = ["SMILES", "ID"]
        output_temp["Representative"] = False
        output_temp.loc[len(output_temp)] = [self.smiles_param, "Custom", True]
        
        #############################~

        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_temp)

##############################################################################################################################################