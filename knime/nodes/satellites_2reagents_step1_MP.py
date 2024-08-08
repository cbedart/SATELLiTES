##############################################################################################################################################
# SATELLiTES - KNIME extension
# Copyright (C) Corentin BEDART, 2024
##############################################################################################################################################

import knime.extension as knext
from satellites_extension import *

import logging
import itertools
import glob 

import os
import pandas as pd
import multiprocessing as mp
import shutil

import SATELLiTES.SATELLiTES_enumeration as stenum
import SATELLiTES.SATELLiTES_gui as stgui

##############################################################################################################################################

@knext.node(name="2-reagents - Step #1 - Multiprocessing", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_multiprocessing)
@knext.input_table(name="Reagents A", description="Reagents A")
@knext.input_table(name="Reagents B", description="Reagents B")
@knext.output_table(name="Enumerated Ab/aB compounds", description="Table of enumerated Ab/aB compounds.")

class Enumeration_2reagents_step1_MP:
    """ **SATELLiTES MP - 2-reagents - Step #1 - **
    This node is the **multiprocessing** alternative of the \"2-reagents - Step #1\" node.

    This node produces the enumeration of compounds Ab/aB by the combination of reagents A/B with the defined representatives reagents a/b. 

    As inputs tables, you need to provide:

    - The **reagents A table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.
    - The **reagents B table**, with the Smiles, ID, and the **Representative column** obtained with either the \"Smallest\", \"By ID\", or \"Custom\" representative nodes.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds Ab and aB** will be generated.
    The selection of these best compounds is determined by the method(s) of your choice. The next step will include the \"2-reagents - Step #2\" node, as well as the "Chosen compounds from previous step" node.
    """
        
    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")
    nb_cpu_max = mp.cpu_count()
    nb_cpu = knext.IntParameter("Number of CPUs", "Number of CPUs", default_value = 1, min_value = 1, max_value = nb_cpu_max)
    groups_of = knext.IntParameter("Split into groups of X compounds", "Split into groups of X compounds", default_value = 20, min_value = 1)
    
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

        stgui.CB_BB_splitter("temp", input_1_temp, "A", groupof=self.groups_of, gui=0)
        reagentA_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_1_temp")
        reagentA_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_1_temp/" + i for i in reagentA_list]

        stgui.CB_BB_splitter("temp", input_2_temp, "B", groupof=self.groups_of, gui=0)
        reagentB_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_2_temp")
        reagentB_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_2_temp/" + i for i in reagentB_list]

        #############################~

        arguments_list_Ab = []
        products_combinations_Ab = list(itertools.product(reagentA_list, [representative_2_temp]))
        products_output_names_Ab = ["Ab1_{}_enum.smi".format(i) for i in range(1, len(products_combinations_Ab) + 1)]
        for i in range(len(products_output_names_Ab)):
            arguments_list_Ab.append((products_combinations_Ab[i][0], products_combinations_Ab[i][1], self.reaction_param, products_output_names_Ab[i], workdir_temp, workdir_temp, "Ab"))

        arguments_list_aB = []
        products_combinations_aB = list(itertools.product([representative_1_temp], reagentB_list))
        products_output_names_aB = ["aB2_{}_enum.smi".format(i) for i in range(1, len(products_combinations_aB) + 1)]
        for i in range(len(products_output_names_Ab)):
            arguments_list_aB.append((products_combinations_aB[i][0], products_combinations_aB[i][1], self.reaction_param, products_output_names_aB[i], workdir_temp, workdir_temp, "aB"))

        #############################~

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_2reagents, arguments_list_Ab):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_2reagents, arguments_list_aB):
                print(result)

        #############################~

        files_Ab = glob.glob(workdir_temp + "/Ab1_*_enum.smi")
        files_aB = glob.glob(workdir_temp + "/aB2_*_enum.smi")
        all_files = files_Ab + files_aB

        all_files_pd = [pd.read_csv(i, sep = "\t", names=["Smiles", "ID"]) for i in all_files]
        output_cpds = pd.concat(all_files_pd).reset_index(drop = True)

        #############################~

        shutil.rmtree("temp_results/", ignore_errors=True)
        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################

