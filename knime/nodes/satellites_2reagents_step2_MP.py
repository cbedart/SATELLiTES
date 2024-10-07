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

@knext.node(name="2-reagents - Step #2 - Multiprocessing", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_multiprocessing)
@knext.input_port(name="Chosen compounds from previous step", description="Input of chosen compounds from the previous step, using the node \"Chosen compounds from previous step\".", port_type=port_type_SATELLiTES)
@knext.input_table(name="Reagents A", description="Input of reagents A.")
@knext.input_table(name="Reagents B", description="Input of reagents B.")
@knext.output_table(name="Enumerated AB compounds", description="Table of final and fully-enumerated AB compounds.")

class Enumeration_2reagents_step2_MP:
    """ **SATELLiTES MP - 2-reagents - Step #2**
    This node is the **multiprocessing** alternative of the \"2-reagents - Step #2\" node.

    This node produces the enumeration of final and fully-enumerated AB compounds by the combination of reagents A/B with the chosen Ab/aB compounds from the previous step. 

    As inputs tables, you need to provide:
    
    - The **chosen Ab/aB compounds from the previous step** using the squared input port, obtained from the \"Chosen compounds from the previous step\" node.
    - The **reagents A table**, with the Smiles and ID.
    - The **reagents B table**, with the Smiles and ID.

    The node will generate a new table with all the compounds enumerated. For this specific node, **compounds AB** will be generated. These are the final and fully-assembled compounds. The selection of the best final compounds is determined by the method(s) of your choice.
    """

    reaction_param = knext.StringParameter("Reaction", "Reaction SMARTS string")
    nb_cpu_max = mp.cpu_count()
    nb_cpu = knext.IntParameter("Number of CPUs", "Number of CPUs", default_value = 1, min_value = 1, max_value = nb_cpu_max)
    groups_of = knext.IntParameter("Split into groups of X compounds", "Split into groups of X compounds", default_value = 20, min_value = 1)
    
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

        stgui.CB_BB_splitter("temp", input_1_temp, "A", groupof=self.groups_of, gui=0)
        reagentA_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_1_temp")
        reagentA_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_1_temp/" + i for i in reagentA_list]

        stgui.CB_BB_splitter("temp", input_2_temp, "B", groupof=self.groups_of, gui=0)
        reagentB_list = os.listdir("./temp_results/building_blocks_forEnumeration/input_2_temp")
        reagentB_list = [workdir + "/temp_results/building_blocks_forEnumeration/input_2_temp/" + i for i in reagentB_list]

        #############################~

        arguments_list_Ab = []
        products_combinations_Ab = list(itertools.product(reagentA_list, [input_0_temp_aB]))
        products_output_names_Ab = ["AB1_{}_enum_step2.smi".format(i) for i in range(1, len(products_combinations_Ab) + 1)]
        for i in range(len(products_output_names_Ab)):
            arguments_list_Ab.append((products_combinations_Ab[i][0], products_combinations_Ab[i][1], self.reaction_param, products_output_names_Ab[i], workdir_temp, workdir_temp, "AB"))

        arguments_list_aB = []
        products_combinations_aB = list(itertools.product([input_0_temp_Ab], reagentB_list))
        products_output_names_aB = ["AB2_{}_enum_step2.smi".format(i) for i in range(1, len(products_combinations_aB) + 1)]
        for i in range(len(products_output_names_Ab)):
            arguments_list_aB.append((products_combinations_aB[i][0], products_combinations_aB[i][1], self.reaction_param, products_output_names_aB[i], workdir_temp, workdir_temp, "AB"))

        #############################~

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_2reagents, arguments_list_Ab):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_2reagents, arguments_list_aB):
                print(result)

        #############################~

        files_AB1 = glob.glob(workdir_temp + "/AB1_*_enum_step2.smi")
        files_AB2 = glob.glob(workdir_temp + "/AB2_*_enum_step2.smi")
        all_files = files_AB1 + files_AB2

        all_files_pd = [pd.read_csv(i, sep = "\t", names=["Smiles", "ID"]) for i in all_files]
        output_cpds = pd.concat(all_files_pd).reset_index(drop = True)

        #############################~

        shutil.rmtree("temp_results/", ignore_errors=True)
        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)

##############################################################################################################################################