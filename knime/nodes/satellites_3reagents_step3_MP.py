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

@knext.node(name="3-reagents - Step #3 - Multiprocessing", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category_multiprocessing)
@knext.input_port(name="Chosen compounds from previous step", description="Input of chosen compounds from the previous step, using the node \"Chosen compounds from previous step\".", port_type=port_type_SATELLiTES)
@knext.input_table(name="Reagents A", description="Input of reagents A with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents B", description="Input of reagents B with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.input_table(name="Reagents C", description="Input of reagents C with the reagent representative defined using one compatible nodes (either smallest, by ID, or custom).")
@knext.output_table(name="Enumerated ABC compounds", description="Table of final and fully-enumerated ABC compounds.")

class Enumeration_3reagents_step3_MP:
    """ **SATELLiTES MP - 3-reagents - Step #3 **
    This node is the **multiprocessing** alternative of the \"3-reagents - Step #3\" node.

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
    nb_cpu_max = mp.cpu_count()
    nb_cpu = knext.IntParameter("Number of CPUs", "Number of CPUs", default_value = 1, min_value = 1, max_value = nb_cpu_max)
    groups_of = knext.IntParameter("Split into groups of X compounds", "Split into groups of X compounds", default_value = 20, min_value = 1)

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
                file_output_AB_B.write("{0}\t{1}\n".format(input_1_pandas.loc[input_1_pandas.iloc[:,1] == input_0_pandas_AB_B_ID].iloc[0,0], input_0_pandas_AB_B_ID))

        #####~

        input_0_temp_AC_A = workdir_temp + "/input_0_temp_AC_A.smi"
        input_0_pandas_AC_A = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("AbC")]
        input_0_pandas_AC_A_IDs = input_0_pandas_AC_A.iloc[:,1].str.split("-", expand=True).iloc[:,1].to_list()
        with open(input_0_temp_AC_A, "w") as file_output_AC_A:
            for input_0_pandas_AC_A_ID in input_0_pandas_AC_A_IDs:
                file_output_AC_A.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_AC_A_ID].iloc[0,0], input_0_pandas_AC_A_ID))

        input_0_temp_AC_C = workdir_temp + "/input_0_temp_AC_C.smi"
        input_0_pandas_AC_C = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("AbC")]
        input_0_pandas_AC_C_IDs = input_0_pandas_AC_C.iloc[:,1].str.split("-", expand=True).iloc[:,3].to_list()
        with open(input_0_temp_AC_C, "w") as file_output_AC_C:
            for input_0_pandas_AC_C_ID in input_0_pandas_AC_C_IDs:
                file_output_AC_C.write("{0}\t{1}\n".format(input_2_pandas.loc[input_2_pandas.iloc[:,1] == input_0_pandas_AC_C_ID].iloc[0,0], input_0_pandas_AC_C_ID))

        #####~

        input_0_temp_BC_B = workdir_temp + "/input_0_temp_BC_B.smi"
        input_0_pandas_BC_B = input_0_pandas.loc[input_0_pandas.iloc[:,1].str.contains("aBC")]
        input_0_pandas_BC_B_IDs = input_0_pandas_BC_B.iloc[:,1].str.split("-", expand=True).iloc[:,2].to_list()
        with open(input_0_temp_BC_B, "w") as file_output_BC_B:
            for input_0_pandas_BC_B_ID in input_0_pandas_BC_B_IDs:
                file_output_BC_B.write("{0}\t{1}\n".format(input_3_pandas.loc[input_3_pandas.iloc[:,1] == input_0_pandas_BC_B_ID].iloc[0,0], input_0_pandas_BC_B_ID))

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

        arguments_list_ABC1 = []
        products_combinations_ABC1 = list(itertools.product([input_0_temp_AB_A], [input_0_temp_AB_B], reagentC_list))
        products_output_names_ABC1 = ["ABC1_{}_enum.smi".format(i) for i in range(1, len(products_combinations_ABC1) + 1)]
        for i in range(len(products_output_names_ABC1)):
            arguments_list_ABC1.append((products_combinations_ABC1[i][0], products_combinations_ABC1[i][1], products_combinations_ABC1[i][2], self.reaction_param, products_output_names_ABC1[i], workdir_temp, workdir_temp, "ABC"))

        arguments_list_ABC2 = []
        products_combinations_ABC2 = list(itertools.product([input_0_temp_AC_A], reagentB_list, [input_0_temp_AC_C]))
        products_output_names_ABC2 = ["ABC2_{}_enum.smi".format(i) for i in range(1, len(products_combinations_ABC2) + 1)]
        for i in range(len(products_output_names_ABC2)):
            arguments_list_ABC2.append((products_combinations_ABC2[i][0], products_combinations_ABC2[i][1], products_combinations_ABC2[i][2], self.reaction_param, products_output_names_ABC2[i], workdir_temp, workdir_temp, "ABC"))

        arguments_list_ABC3 = []
        products_combinations_ABC3 = list(itertools.product(reagentA_list, [input_0_temp_BC_B], [input_0_temp_BC_C]))
        products_output_names_ABC3 = ["ABC3_{}_enum.smi".format(i) for i in range(1, len(products_combinations_ABC3) + 1)]
        for i in range(len(products_output_names_ABC3)):
            arguments_list_ABC3.append((products_combinations_ABC3[i][0], products_combinations_ABC3[i][1], products_combinations_ABC3[i][2], self.reaction_param, products_output_names_ABC3[i], workdir_temp, workdir_temp, "ABC"))

        #############################~

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_ABC1):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_ABC2):
                print(result)

        with mp.Pool(self.nb_cpu) as pool:
            for result in pool.starmap(stenum.run_SATELLiTES_3reagents, arguments_list_ABC3):
                print(result)

        #############################~

        files_ABC1 = glob.glob(workdir_temp + "/ABC1_*_enum.smi")
        files_ABC2 = glob.glob(workdir_temp + "/ABC2_*_enum.smi")
        files_ABC3 = glob.glob(workdir_temp + "/ABC3_*_enum.smi")
        all_files = files_ABC1 + files_ABC2 + files_ABC3

        all_files_pd = [pd.read_csv(i, sep = "\t", names=["Smiles", "ID"]) for i in all_files]
        output_cpds = pd.concat(all_files_pd).reset_index(drop = True)

        #############################~

        shutil.rmtree("temp_results/", ignore_errors=True)
        shutil.rmtree("tmp/", ignore_errors=True)

        #############################~

        return knext.Table.from_pandas(output_cpds)
    
##############################################################################################################################################