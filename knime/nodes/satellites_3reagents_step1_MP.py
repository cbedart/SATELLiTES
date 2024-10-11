##############################################################################################################################################
# SATELLiTES - KNIME extension
# Copyright (C) Corentin BEDART, 2024
##############################################################################################################################################

import knime.extension as knext
from satellites_extension import *

import logging
import itertools
import glob 
import pickle
import os
import pandas as pd
import multiprocessing as mp
import shutil

import SATELLiTES.SATELLiTES_enumeration as stenum
import SATELLiTES.SATELLiTES_gui as stgui

##############################################################################################################################################

@knext.node(name="3-reagents - Step #1 - Multiprocessing", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_multiprocessing)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents C", description="Input of reagents C with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated Abc/aBc/abC compounds", description="Table of enumerated Abc/aBc/abC compounds.")

class Enumeration_3reagents_step1_MP:
    """ **SATELLiTES MP - 3-reagents - Step #1**
    This node is the **multiprocessing** alternative of the \"3-reagents - Step #1\" node.

    This node produces the enumeration of compounds Abc/aBc/abC by the combination of reagents A/B/C with the defined representatives reagents a/b/c. 

    As inputs tables, you need to provide:

    - The **reagents A table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents B table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents C table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds Abc, aBc, and abC** will be generated.

    The selection of these best compounds is determined by the method(s) of your choice. The next step will include the \"3-reagents - Step #2\" node, as well as the "Chosen compounds from previous step" node.
    """

    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")
    nb_cpu_max = mp.cpu_count()
    nb_cpu = knext.IntParameter("Number of CPUs", "Number of CPUs", default_value = 1, min_value = 1, max_value = nb_cpu_max)
    groups_of = knext.IntParameter("Split into groups of X compounds", "Split into groups of X compounds", default_value = 20, min_value = 1)

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

        stgui.CB_BB_splitter("temp", input_1_temp, "A", groupof=self.groups_of, gui=0)
        reagentA_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_1_temp")
        reagentA_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_1_temp/" + i for i in reagentA_list]

        stgui.CB_BB_splitter("temp", input_2_temp, "B", groupof=self.groups_of, gui=0)
        reagentB_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_2_temp")
        reagentB_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_2_temp/" + i for i in reagentB_list]

        stgui.CB_BB_splitter("temp", input_3_temp, "C", groupof=self.groups_of, gui=0)
        reagentC_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_3_temp")
        reagentC_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_3_temp/" + i for i in reagentC_list]

        #############################~

        arguments_list_Abc1 = []
        products_combinations_Abc1 = list(itertools.product(reagentA_list, [representative_2_temp], [representative_3_temp]))
        products_output_names_Abc1 = ["Abc1_{}_enum.smi".format(i) for i in range(1, len(products_combinations_Abc1) + 1)]
        for i in range(len(products_output_names_Abc1)):
            arguments_list_Abc1.append((products_combinations_Abc1[i][0], products_combinations_Abc1[i][1], products_combinations_Abc1[i][2], self.reaction_param, products_output_names_Abc1[i], workdir_temp, workdir_temp, "Abc"))

        arguments_list_aBc2 = []
        products_combinations_aBc2 = list(itertools.product([representative_1_temp], reagentB_list, [representative_3_temp]))
        products_output_names_aBc2 = ["aBc2_{}_enum.smi".format(i) for i in range(1, len(products_combinations_aBc2) + 1)]
        for i in range(len(products_output_names_aBc2)):
            arguments_list_aBc2.append((products_combinations_aBc2[i][0], products_combinations_aBc2[i][1], products_combinations_aBc2[i][2], self.reaction_param, products_output_names_aBc2[i], workdir_temp, workdir_temp, "aBc"))

        arguments_list_abC3 = []
        products_combinations_abC3 = list(itertools.product([representative_1_temp], [representative_2_temp], reagentC_list))
        products_output_names_abC3 = ["abC3_{}_enum.smi".format(i) for i in range(1, len(products_combinations_abC3) + 1)]
        for i in range(len(products_output_names_abC3)):
            arguments_list_abC3.append((products_combinations_abC3[i][0], products_combinations_abC3[i][1], products_combinations_abC3[i][2], self.reaction_param, products_output_names_abC3[i], workdir_temp, workdir_temp, "abC"))

        #############################~

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_Abc1):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_aBc2):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_abC3):
                print(result)

        #############################~

        files_Abc1 = glob.glob(workdir_temp + "/Abc1_*_enum.smi")
        files_aBc2 = glob.glob(workdir_temp + "/aBc2_*_enum.smi")
        files_abC3 = glob.glob(workdir_temp + "/abC3_*_enum.smi")
        all_files = files_Abc1 + files_aBc2 + files_abC3

        all_files_pd = [pd.read_csv(i, sep = "\t", names=["Smiles", "ID"]) for i in all_files]
        output_cpds = pd.concat(all_files_pd).reset_index(drop = True)

        #############################~

        shutil.rmtree("temp_results/", ignore_errors=True)
        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################