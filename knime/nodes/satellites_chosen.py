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

@knext.node(name="Chosen compounds from previous step", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_utilities)
@knext.input_table(name="Best compounds", description="Best compounds")
@knext.output_port(name="Output Data", description="Whatever the node has produced", port_type=port_type_SATELLiTES)

class Enumeration_chosen:
    """ **SATELLiTES - Chosen compounds from previous step**
    The **\"Chosen compounds from previous step\"** creates the input that will be used for the **\"X-reagents\"** set of nodes starting steps #2.
    
    As an input table, you need to provide **at least the SATELLiTES ID of the chosen compounds in the first two columns**. It can also take as input a .smi/.csv/.tsv/other file using the \“File Reader\” node,as long as the ID is in the first 2 columns.

    The next step will include the \"X-reagents - Step #2/3\" set of nodes.
    """
        
    def configure(self, configure_context, input_1):
        pass

    def execute(self, execute_context, input_1):
        input_1_pandas = input_1.to_pandas()
        check = 0

        #############################~
        # Check if the IDs are in the first column or without a SMILES string

        if "Ab" in input_1_pandas.iloc[0,0] or "Ba" in input_1_pandas.iloc[0,0]:
            input_1_IDs = input_1_pandas.iloc[:,0]
        else:
            check += 1

        #####~

        if "Abc" in input_1_pandas.iloc[0,0] or "aBc" in input_1_pandas.iloc[0,0] or "abC" in input_1_pandas.iloc[0,0]:
            input_1_IDs = input_1_pandas.iloc[:,0]
        else:
            check += 1

        #####~

        if "ABc" in input_1_pandas.iloc[0,0] or "AbC" in input_1_pandas.iloc[0,0] or "aBC" in input_1_pandas.iloc[0,0]:
            input_1_IDs = input_1_pandas.iloc[:,0]
        else:
            check += 1


        #############################~
        # Check if the IDs are in the second column, before checking if there is a second column
        if input_1_pandas.shape[1] > 1:

            if "Ab" in input_1_pandas.iloc[0,1] or "Ba" in input_1_pandas.iloc[0,1]:
                input_1_IDs = input_1_pandas.iloc[:,1]
            else:
                check += 1

            #####~

            if "Abc" in input_1_pandas.iloc[0,1] or "aBc" in input_1_pandas.iloc[0,1] or "abC" in input_1_pandas.iloc[0,1]:
                input_1_IDs = input_1_pandas.iloc[:,1]
            else:
                check += 1

            #####~

            if "ABc" in input_1_pandas.iloc[0,1] or "AbC" in input_1_pandas.iloc[0,1] or "aBC" in input_1_pandas.iloc[0,1]:
                input_1_IDs = input_1_pandas.iloc[:,1]
            else:
                check += 1

        else:
            check += 3

        #############################~

        # Error if no IDs
        if check == 6:
            raise TypeError("You need to give SATELLiTES compound IDs, with tags at the start of the IDs.")

        #############################~

        # Check the location of SMILES if existing
        try:
            stenum.Chem.AddHs(stenum.Chem.MolFromSmiles(input_1_pandas.iloc[0,0]))
            input_1_SMILES = input_1_pandas.iloc[:,0]
        except:
            try:
                stenum.Chem.AddHs(stenum.Chem.MolFromSmiles(input_1_pandas.iloc[0,1]))
                input_1_SMILES = input_1_pandas.iloc[:,1]
            except:
                input_1_SMILES = ["C" for i in range(len(input_1_IDs))]
        
        #############################~

        output_pandas = pd.DataFrame.from_dict({"Smiles": input_1_SMILES, "ID":input_1_IDs})
        output_objspec = SatellitesObjectSpec(output_pandas)
        output_obj = SatellitesObject(output_objspec, output_pandas)

        #############################~

        return output_obj

##############################################################################################################################################

