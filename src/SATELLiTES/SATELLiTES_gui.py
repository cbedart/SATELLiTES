####################################################################################################
# SATELLiTES - version 1.0.10
# Copyright (C) 2024 Corentin BEDART
####################################################################################################

import multiprocessing as mp
import os
import sys
import itertools
import datetime
import glob
import time
import re
import subprocess
import numpy as np

from rdkit import Chem

from SATELLiTES.SATELLiTES_enumeration import *

####################################################################################################

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

#####

def path_completer(text, state):
    if '~' in text:
        text = os.path.expanduser('~')
    return [x for x in glob.glob(text+'*')][state]

#####

def CB_filetype_detection(first_line, second_line, filename):
    # Check the type of separator
    if len(first_line.strip().split("\t")) > 1 and len(second_line.strip().split("\t")) > 1:
        separator = "\t"
    elif len(first_line.strip().split(" ")) > 1 and len(second_line.strip().split(" ")) > 1:
        separator = " "
    elif len(first_line.strip().split(",")) > 1 and len(second_line.strip().split(",")) > 1:
        separator = ","
    else:
        exit("Wrong file format for {} - we need a tab, space or comma-separated file with the SMILES string in first or second position".format(filename))

    # Check if the SMILES is in the first or second position
    try:
        Chem.AddHs(Chem.MolFromSmiles(second_line.strip().split(separator)[0]))
        smiles_position = 0
        id_position = 1
    except:
        try:
            Chem.AddHs(Chem.MolFromSmiles(second_line.strip().split(separator)[1]))
            smiles_position = 1
            id_position = 0
        except:
            exit("Wrong file format for {} - we need a tab, space or comma-separated file with the SMILES string in first or second position".format(filename))

    # Check if the first line is a header or not
    try:
        Chem.AddHs(Chem.MolFromSmiles(first_line.strip().split(separator)[0]))
        skip_first = 0
    except:
        skip_first = 1

    return (separator, skip_first, smiles_position, id_position)

#####


def CB_smallest_finder(reagent_filename, job_name, tag, gui = 0):
    smallest_MW = 99999
    smallest_line = "".join(str(list(range(99))))

    with open(reagent_filename, "r") as BBs, open(job_name + "_results/" + reagent_filename.split("/")[-1].split(".")[0] + "_smallest.smi", "w") as BBs_smallest:
        BBs_lines = BBs.readlines()
        separator, skip_first, smiles_position, id_position = CB_filetype_detection(BBs_lines[0], BBs_lines[1], reagent_filename)

        for line in BBs_lines[skip_first:]:
            if len(line) > len(smallest_line) + 5:
                continue
            
            mw_temp = Chem.rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(line.strip().split(separator)[smiles_position]), onlyHeavy=True)
            if mw_temp < smallest_MW:
                smallest_MW = mw_temp
                smallest_line = line
        
        smallest_line_split = smallest_line.strip().split(separator)
        BBs_smallest.write("{}\t{}\n".format(smallest_line_split[smiles_position], smallest_line_split[id_position]))

    if gui == 0:
        print("The chosen " + tag + " representative :")
        print("MW = " + str(round(smallest_MW,2)) + " Da")
        print(smallest_line + "\n")
    else:
        gui.update_text("The chosen " + tag + " representative :\n")
        gui.update_text("MW = " + str(round(smallest_MW,2)) + " Da\n")
        gui.update_text(smallest_line + "\n")

#####

def CB_smallest_ID(reagent_filename, job_name, tag, id=0, gui = 0):

    reagent_name = reagent_filename.split("/")[-1].split(".")[0]

    if gui == 0:
        id_check = 0
        while id_check == 0:
            id = input("Your representative " + reagent_name + " as ID (Must be in your " + reagent_name +".smi file): ")
            if id == "":
                continue
            else:
                id_input_check = input("Is your ID " + id + " correct (Y/N): ")
                if id_input_check.upper() == "Y" or id_input_check == "":
                    smallest_line = ""
                    with open(reagent_filename, "r") as file_input:
                        lines = file_input.readlines()
                        separator, skip_first, smiles_position, id_position = CB_filetype_detection(lines[0], lines[1], reagent_filename)
                        for line in lines:
                            if re.search(id, line):
                                smallest_line = line

                        if smallest_line == "":
                            print("\n    Your ID " + id + " was not found in " + reagent_filename)
                            continue
                        else:
                            id_check = 1
                else:
                    continue
    else:
        # GUI mode
        # Need error statement
        smallest_line = ""
        with open(reagent_filename, "r") as file_input:
            lines = file_input.readlines()
            separator, skip_first, smiles_position, id_position = CB_filetype_detection(lines[0], lines[1], reagent_filename)
            for line in lines:
                if re.search(id, line):
                    smallest_line = line

            if smallest_line == "":
                print("\n    Your ID " + id + " was not found in " + reagent_filename)


    smiles_loaded = Chem.AddHs(Chem.MolFromSmiles(smallest_line.strip().split(separator)[smiles_position]))

    with open(job_name + "_results/" + reagent_name + "_representative.smi", "w") as BBs_smallest:
        BBs_smallest.write("{}\t{}\n".format(smallest_line.strip().split(separator)[smiles_position], smallest_line.strip().split(separator)[id_position]))

    if gui == 0:
        print("\nYour chosen " + reagent_name + " representative:")
        print("MW = " + str(round(Chem.rdMolDescriptors.CalcExactMolWt(smiles_loaded),2)) + " Da")
        print(smallest_line + "\n")
    else:
        gui.update_text("Your chosen " + reagent_name + " representative:\n")
        gui.update_text("MW = " + str(round(Chem.rdMolDescriptors.CalcExactMolWt(smiles_loaded),2)) + " Da\n")
        gui.update_text(smallest_line + "\n")

#####

def CB_smallest_smiles(reagent_filename, job_name, tag, smallest_smiles = 0, gui = 0):
    reagent_name = reagent_filename.split("/")[-1].split(".")[0]

    if gui == 0:
        smiles_check = 0
        while smiles_check == 0:
            smallest_smiles = input("Your representative " + reagent_name + " as Smiles string: ")
            smallest_name = input("Your representative " + reagent_name + " name (Default 'Custom" + reagent_name + "'): ")

            if smallest_name == "":
                smallest_name = "Custom" + tag

            if smallest_smiles == "":
                continue
            else:
                smallest_smiles_check = input("Is your " + smallest_name + " Smiles " + smallest_smiles + " correct (Y/N): ")
                if smallest_smiles_check.upper() == "Y" or smallest_smiles_check == "":
                    try:
                        smiles_loaded = Chem.AddHs(Chem.MolFromSmiles(smallest_smiles))
                        smiles_check = 1
                    except:
                        print("\n    RDKit error with your smiles " + smallest_smiles)
                else:
                    continue
    else:
        # Need an error statement in here if the input SMILES is wrong
        # GUI mode
        smallest_name = "Custom" + tag
        try:
            smiles_loaded = Chem.AddHs(Chem.MolFromSmiles(smallest_smiles))
            smiles_check = 1
        except:
            print("\n    RDKit error with your smiles " + smallest_smiles)

    with open(job_name + "_results/" + reagent_name + "_representative.smi", "w") as BBs_smallest:
        smallest_line = "{}\t{}\n".format(smallest_smiles, smallest_name)
        BBs_smallest.write(smallest_line)

    if gui == 0:
        print("\nYour chosen '" + smallest_name + "' " + reagent_name + " representative:")
        print("MW = " + str(round(Chem.rdMolDescriptors.CalcExactMolWt(smiles_loaded),2)) + " Da")
        print(smallest_line + "\n")
    else:
        gui.update_text("Your chosen '" + smallest_name + "' " + reagent_name + " representative:\n")
        gui.update_text("MW = " + str(round(Chem.rdMolDescriptors.CalcExactMolWt(smiles_loaded),2)) + " Da\n")
        gui.update_text(smallest_line + "\n")


#####

def CB_BB_splitter(job_name, reagent_filename, tag, groupof=400, gui=0):
    reagent_name = reagent_filename.split("/")[-1].split(".")[0]

    # Create main BB directory
    if not os.path.exists(job_name + "_results/building_blocks_forEnumeration"):
        os.makedirs(job_name + "_results/building_blocks_forEnumeration")

    # Create sub directory with reagent name
    if not os.path.exists(job_name + "_results/building_blocks_forEnumeration/" + reagent_name):
        os.makedirs(job_name + "_results/building_blocks_forEnumeration/" + reagent_name)

    # Read input smiles
    with open(reagent_filename, "r") as smiles_input:
        lines = smiles_input.readlines()
        separator, skip_first, smiles_position, id_position = CB_filetype_detection(lines[0], lines[1], reagent_filename)
        nb_lines = len(lines[skip_first:])
        lines = lines[skip_first:]

    splitA = list(range(0, nb_lines, groupof)) + [nb_lines] # ; splitA = splitA[0:nbgroups+1]

    for i in list(range(len(splitA) - 1)):
        with open("./" + job_name + "_results/building_blocks_forEnumeration/" + reagent_name + "/" + reagent_name + "_group" + str(i+1) + ".smi", "w") as write_output:
            for j in range(splitA[i], splitA[i+1]):
                write_output.write("{}\t{}\n".format(lines[j].strip().split(separator)[smiles_position], lines[j].strip().split(separator)[id_position]))
    
    if gui == 0:
        print(" - {2}{0}{4} from {1}: {2}{3} group(s){4} of maximum {5} building blocks\n".format(reagent_name, reagent_filename, color.BOLD + color.YELLOW, len(splitA)-1, color.END, groupof))
    else:
        gui.update_text(" - <b>{0}</b>: <b>{1} group(s)</b> of maximum {2} building blocks\n".format(reagent_name, len(splitA)-1, groupof))



#####

def CB_reaction_2R(job_name, reagents_list, reaction, reagents_smaller, tag, nb_cpu, output_path, app = 0, gui = 0):
    # Var init
    products_output_names = []

    arguments_list = []

    # Output library path
    # output_library_path = job_name + "_results/output_libraries_" + tag
    output_library_path = job_name + "_results/output_libraries"
    if not os.path.exists(output_library_path):
        os.makedirs(output_library_path)

    # Generation of combinations of reagentsA vs reagentsB
    products_combinations = list(itertools.product(reagents_list, reagents_smaller))

    # Output filenames - smi format for compounds
    for combination in products_combinations:
        products_output_names.append(tag + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + ".smi")

    # # Output filenames - csv format for statistics
    # for combination in products_combinations:
    #     products_output_stats.append(tag + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + "_" + tag + ".csv")

    # Generation of the big argument list to run python multiprocessing
    # Specific case if Ba = products_combinations must be inverted, as it's always AB in the reaction for now
    if tag == "Ab":
        for i in range(len(products_output_names)):
            arguments_list.append((products_combinations[i][0], products_combinations[i][1], reaction, products_output_names[i], output_library_path, output_path, tag))
    elif tag == "Ba":
        for i in range(len(products_output_names)):
            arguments_list.append((products_combinations[i][1], products_combinations[i][0], reaction, products_output_names[i], output_library_path, output_path, tag))
    
    # Logfile generation to benchmark, and see if there is some issues just in case 
    # with open(job_name + "_results/logfile.txt", "a") as logfile_in:
    #     logfile_in.write("Logfile for enumeration of " + tag + " products - Start : " + str(datetime.datetime.now()) + "\n")

    # Main multiprocessing loop
    with mp.Pool(nb_cpu) as pool:
        for result in pool.starmap(run_SATELLiTES_2reagents, arguments_list):
            if gui == 0:
                print(result)
            else:
                gui.update_text(result + "\n")

####################################################################################################

def CB_selection_step2_2R(reagentsA, reagentsB, job_name, pathA, pathB, gui = 0):
    if gui == 0:
        # # Tab completer
        # readline.set_completer_delims('\t')
        # readline.parse_and_bind("tab: complete")
        # readline.set_completer(path_completer)

        # step2_path_check = 0
        # while step2_path_check == 0:
        #     step2_filename = input(color.BOLD + color_ansi + name + " path (.smi file): " + color.END)
        #     if os.path.exists(step2_filename):
        #         if step2_filename[-4:] == ".smi":
        #             step2_path_check = 1
        #             cmd_copy = "cp " + step2_filename[0] + " " + job_name + "_results/."
        #             subprocess.run(cmd_copy, shell=True)
        #             output_path = ["./" + job_name + "_results/" + step2_filename[0].split("/")[-1]]
        #         else:
        #             print("\n    Not a .smi file")
        #     else:
        #         print("\n    File does not exist")

        exit("TO DO")

    else:
        # GUI mode ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # Errors to do, if does not exist (should not arrive) or in case of an header

        gui.update_text("Focused Enumeration Library - Extraction of the best selected reagents\n")

        #####

        output_path_A = job_name + "_results/selected_reagentsA.smi"
        output_path_B = job_name + "_results/selected_reagentsB.smi"

        #####

        selected_reagentsA = []
        selected_reagentsB = []

        with open(pathA, "r") as filein:
            for line in filein.readlines():
                line_splitted = line.rstrip().split("\t")
                id_splitted = line_splitted[-1].split("-")
                if id_splitted[0] == "Ab":
                    selected_reagentsA.append(id_splitted[1])
                if id_splitted[0] == "Ba":
                    selected_reagentsB.append(id_splitted[2])

        with open(pathB, "r") as filein:
            for line in filein.readlines():
                line_splitted = line.rstrip().split("\t")
                id_splitted = line_splitted[-1].split("-")
                if id_splitted[0] == "Ab":
                    selected_reagentsA.append(id_splitted[1])
                if id_splitted[0] == "Ba":
                    selected_reagentsB.append(id_splitted[2])

        selected_reagentsA_unique_iter = np.unique(selected_reagentsA).tolist()
        selected_reagentsB_unique_iter = np.unique(selected_reagentsB).tolist()

        #####

        with open(reagentsA, "r") as file_input, open(output_path_A, "w") as file_output:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsA_unique_iter:
                    file_output.write(line)
                    selected_reagentsA_unique_iter.remove(id_temp)
                    if len(selected_reagentsA_unique_iter) == 0:
                        break
        
        gui.update_text("Best reagents A copied successfully in the results directory - selected_reagentsA.smi\n")

        #####

        with open(reagentsB, "r") as file_input, open(output_path_B, "w") as file_output:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsB_unique_iter:
                    file_output.write(line)
                    selected_reagentsB_unique_iter.remove(id_temp)
                    if len(selected_reagentsB_unique_iter) == 0:
                        break

        gui.update_text("Best reagents B copied successfully in the results directory - selected_reagentsB.smi\n")

        #####

    return [output_path_A], [output_path_B]


####################################################################################################

def CB_run_2reagents(reagentsA, reagentsB, reaction, output_path, nb_cpu, modeA, modeB, inputA, inputB, window, app, job_name = "SATELLiTES"):
    start_time = time.time()

    window.update_test()
    # window.update_progress(0)
    window.update_progress(100)
    app.processEvents()

    if not os.path.exists(job_name + "_results"):
        os.makedirs(job_name + "_results")

    # Modes summary
    # 0 = MEL Smallest
    # 1 = MEL ID
    # 2 = MEL Smiles
    # 3 = FEL File

    # Reagents A
    if modeA == 0:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Smallest finder</div></b>\n", start=1)
        CB_smallest_finder(reagentsA, job_name, "A", gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_smallest.smi"]
        # window.update_progress(10)
        window.update_progress(100)

    if modeA == 1:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Find by ID</div></b>\n", start=1)
        CB_smallest_ID(reagentsA, job_name, "A", id=inputA, gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # window.update_progress(10)
        window.update_progress(100)

    if modeA == 2:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Custom SMILES</div></b>\n", start=1)
        CB_smallest_smiles(reagentsA, job_name, "A", smallest_smiles=inputA, gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # window.update_progress(10)
        window.update_progress(100)


    # Reagents B
    if modeB == 0:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Smallest finder</div></b>\n")
        CB_smallest_finder(reagentsB, job_name, "B", gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_smallest.smi"]

    if modeB == 1:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Find by ID</div></b>\n")
        CB_smallest_ID(reagentsB, job_name, "B", id=inputB, gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_representative.smi"]

    if modeB == 2:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Custom SMILES</div></b>\n")
        CB_smallest_smiles(reagentsB, job_name, "B", smallest_smiles=inputB, gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_representative.smi"]


    # If FEL
    if modeA == 3 and modeB == 3:
        window.update_text("<b><div style='color:green;'>Step #1 & #2 : Selected reagents A/B for the Focused Enumeration Library</div></b>\n", start=1)
        reagentA_smallest, reagentB_smallest = CB_selection_step2_2R(reagentsA, reagentsB, job_name, pathA=inputA, pathB=inputB, gui=window)


    # window.update_progress(20)
    window.update_progress(100)
    app.processEvents()

    # BB split
    window.update_text("<div>\t</div>\n<b><div>Step #3 : Splitting input building blocks</div></b>\n")
    CB_BB_splitter(job_name, reagentsA, "A", gui=window)
    # window.update_progress(30)
    window.update_progress(100)
    app.processEvents()

    CB_BB_splitter(job_name, reagentsB, "B", gui=window)
    # window.update_progress(40)
    window.update_progress(100)
    app.processEvents()

    # Enumeration lists
    reagentA_name = reagentsA.split("/")[-1].split(".")[0]
    reagentB_name = reagentsB.split("/")[-1].split(".")[0]
    
    reagentA_list = os.listdir("./" + job_name + "_results/building_blocks_forEnumeration/" + reagentA_name)
    reagentA_list = ["./" + job_name + "_results/building_blocks_forEnumeration/" + reagentA_name + "/" + i for i in reagentA_list]

    reagentB_list = os.listdir("./" + job_name + "_results/building_blocks_forEnumeration/" + reagentB_name)
    reagentB_list = ["./" + job_name + "_results/building_blocks_forEnumeration/" + reagentB_name + "/" + i for i in reagentB_list]

    ####################################################################################################

    window.update_text("<div>\t</div>\n<b><div style='color:green;'>Step #4 : Enumeration of reagents A with the representative reagent B (Ab)</div></b>\n")

    # Ab
    CB_reaction_2R(job_name, reagentA_list, reaction, reagentB_smallest, "Ab", nb_cpu, output_path, app, gui=window)
    # window.update_progress(70)
    window.update_progress(100)
    app.processEvents()

    ####################################################################################################
    window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #5 : Enumeration of reagents B with the representative reagent A (Ba)</div></b>\n")
    # Ba
    CB_reaction_2R(job_name, reagentB_list, reaction, reagentA_smallest, "Ba", nb_cpu, output_path, app, gui=window)
    window.update_progress(100) ; app.processEvents()

    ####################################################################################################
    end_time = time.time()
    total_time = round(end_time - start_time, 0)
    if total_time >= 60:
        total_time_output = "{} minute(s) and {} second(s)".format(int(total_time//60), int(round(total_time%60, 0)))
    else:
        total_time_output = "{} second(s)".format(int(total_time))

    window.update_text("<div>\t</div>\n<b><div>Enumerations ended successfully in " + total_time_output + " !</div></b>\n\n")
    
    return "ok"


####################################################################################################

def CB_run_3reagents(reagentsA, reagentsB, reagentsC, reaction, output_path, nb_cpu, modeA, modeB, modeC, inputA, inputB, inputC, step, step_path, window, app, job_name = "SATELLiTES"):
    start_time = time.time()

    window.update_test()
    # window.update_progress(0)
    window.update_progress(100)
    app.processEvents()

    if not os.path.exists(job_name + "_results"):
        os.makedirs(job_name + "_results")

    # Modes summary

    # Step #1 - Minimum Enumeration Library - Abc/aBc/abC
    # Step #2 - MEL #2 - ABc/AbC/aBC
    # Step #3 - Focused Enumeration Library

    # 0 = MEL Smallest
    # 1 = MEL ID
    # 2 = MEL Smiles

    # 3 = FEL File

    ##############################################################################################################################################
    # Reagents A
    if modeA == 0:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Smallest finder</div></b>\n", start=1)
        CB_smallest_finder(reagentsA, job_name, "A", gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_smallest.smi"]
        # window.update_progress(10)
        window.update_progress(100)

    if modeA == 1:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Find by ID</div></b>\n", start=1)
        CB_smallest_ID(reagentsA, job_name, "A", id=inputA, gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # window.update_progress(10)
        window.update_progress(100)

    if modeA == 2:
        window.update_text("<b><div style='color:green;'>Step #1 : Representative reagent A - Custom SMILES</div></b>\n", start=1)
        CB_smallest_smiles(reagentsA, job_name, "A", smallest_smiles=inputA, gui=window)
        reagentA_smallest = ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # window.update_progress(10)
        window.update_progress(100)

    ##############################################################################################################################################
    # Reagents B
    if modeB == 0:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Smallest finder</div></b>\n")
        CB_smallest_finder(reagentsB, job_name, "B", gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_smallest.smi"]

    if modeB == 1:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Find by ID</div></b>\n")
        CB_smallest_ID(reagentsB, job_name, "B", id=inputB, gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_representative.smi"]

    if modeB == 2:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent B - Custom SMILES</div></b>\n")
        CB_smallest_smiles(reagentsB, job_name, "B", smallest_smiles=inputB, gui=window)
        reagentB_smallest = ["./" + job_name + "_results/" + reagentsB.split("/")[-1].split(".")[0] + "_representative.smi"]

    ##############################################################################################################################################
    # Reagents C
    if modeC == 0:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent C - Smallest finder</div></b>\n")
        CB_smallest_finder(reagentsC, job_name, "B", gui=window)
        reagentC_smallest = ["./" + job_name + "_results/" + reagentsC.split("/")[-1].split(".")[0] + "_smallest.smi"]

    if modeC == 1:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent C - Find by ID</div></b>\n")
        CB_smallest_ID(reagentsC, job_name, "B", id=inputC, gui=window)
        reagentC_smallest = ["./" + job_name + "_results/" + reagentsC.split("/")[-1].split(".")[0] + "_representative.smi"]

    if modeC == 2:
        window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #2 : Representative reagent C - Custom SMILES</div></b>\n")
        CB_smallest_smiles(reagentsC, job_name, "B", smallest_smiles=inputC, gui=window)
        reagentC_smallest = ["./" + job_name + "_results/" + reagentsC.split("/")[-1].split(".")[0] + "_representative.smi"]

    ##############################################################################################################################################

    # If FEL
    # if modeA == 3 and modeB == 3:
    #     window.update_text("<b><div style='color:green;'>Step #1 & #2 : Selected reagents A/B/C for the Focused Enumeration Library</div></b>\n", start=1)
    #     reagentA_smallest, reagentB_smallest = CB_selection_step2(reagentsA, reagentsB, job_name, pathA=inputA, pathB=inputB, gui=window)

    # window.update_progress(20)
    window.update_progress(100)
    app.processEvents()


    ##############################################################################################################################################
    # If step #2
    if step == 2:
        window.update_text("3-reagents step2 - Extract bests A/B/C\n")
        selected_reagentA_list, selected_reagentB_list, selected_reagentC_list = CB_selection_step2_3R(reagentsA, reagentsB, reagentsC, job_name, step_path, gui=window)

    ##############################################################################################################################################
    # If step #3
    if step == 3:
        # COMBOS TO MANAGE

        os.makedirs(job_name + "_results/temp_combos", exist_ok = True)

        os.makedirs(job_name + "_results/temp_combos/AB", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/AB/A", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/AB/B", exist_ok = True)

        os.makedirs(job_name + "_results/temp_combos/AC", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/AC/A", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/AC/C", exist_ok = True)

        os.makedirs(job_name + "_results/temp_combos/BC", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/BC/B", exist_ok = True)
        os.makedirs(job_name + "_results/temp_combos/BC/C", exist_ok = True)

        reagentsAB_A_paths, reagentsAB_B_paths, reagentsAC_A_paths, reagentsAC_C_paths, reagentsBC_B_paths, reagentsBC_C_paths = CB_selection_step3_3R(reagentsA, reagentsB, reagentsC, job_name, step_path, gui=window)

    ##############################################################################################################################################
    # BB split
    window.update_text("<div>\t</div>\n<b><div>Step #3 : Splitting input building blocks</div></b>\n")
    CB_BB_splitter(job_name, reagentsA, "A", gui=window)
    # window.update_progress(30)
    window.update_progress(100)
    app.processEvents()

    CB_BB_splitter(job_name, reagentsB, "B", gui=window)
    # window.update_progress(40)
    window.update_progress(100)
    app.processEvents()

    CB_BB_splitter(job_name, reagentsC, "C", gui=window)
    # window.update_progress(50)
    window.update_progress(100)
    app.processEvents()

    ##############################################################################################################################################
    # Enumeration lists
    reagentA_name = reagentsA.split("/")[-1].split(".")[0]
    reagentA_list = os.listdir("./" + job_name + "_results/building_blocks_forEnumeration/" + reagentA_name)
    reagentA_list = ["./" + job_name + "_results/building_blocks_forEnumeration/" + reagentA_name + "/" + i for i in reagentA_list]

    #####

    reagentB_name = reagentsB.split("/")[-1].split(".")[0]
    reagentB_list = os.listdir("./" + job_name + "_results/building_blocks_forEnumeration/" + reagentB_name)
    reagentB_list = ["./" + job_name + "_results/building_blocks_forEnumeration/" + reagentB_name + "/" + i for i in reagentB_list]

    #####

    reagentC_name = reagentsC.split("/")[-1].split(".")[0]
    reagentC_list = os.listdir("./" + job_name + "_results/building_blocks_forEnumeration/" + reagentC_name)
    reagentC_list = ["./" + job_name + "_results/building_blocks_forEnumeration/" + reagentC_name + "/" + i for i in reagentC_list]

    ####################################################################################################

    window.update_text("<div>\t</div>\n<b><div style='color:green;'>Step #4 : Enumeration of reagents A with the representative reagent B (Ab)</div></b>\n")

    # Ab
    # CB_reaction(job_name, reagentA_list, reaction, reagentB_smallest, "Ab", nb_cpu, output_path, app, gui=window)
    # window.update_progress(70) ; app.processEvents()

    ####################################################################################################
    window.update_text("<div>\t</div>\n<b><div style='color:red;'>Step #5 : Enumeration of reagents B with the representative reagent A (Ba)</div></b>\n")
    # Ba
    # CB_reaction(job_name, reagentB_list, reaction, reagentA_smallest, "Ba", nb_cpu, output_path, app, gui=window)
    # window.update_progress(100) ; app.processEvents()

    if step == 1:
        # Abc
        CB_reaction_3R(job_name, reaction, reagentA_list, reagentB_smallest, reagentC_smallest, "Abc", nb_cpu, output_path, step,  app, gui=window)
        # aBc
        CB_reaction_3R(job_name, reaction, reagentA_smallest, reagentB_list, reagentC_smallest, "aBc", nb_cpu, output_path, step,  app, gui=window)
        # abC
        CB_reaction_3R(job_name, reaction, reagentA_smallest, reagentB_smallest, reagentC_list, "abC", nb_cpu, output_path, step,  app, gui=window)
        # Merge files
        CB_merge_libraries(job_name, "A")
        CB_merge_libraries(job_name, "B")
        CB_merge_libraries(job_name, "C")

    if step == 2:
        # ABc
        CB_reaction_3R(job_name, reaction, reagentA_list, selected_reagentB_list, reagentC_smallest, "ABc", nb_cpu, output_path, step,  app, gui=window)
        CB_reaction_3R(job_name, reaction, selected_reagentA_list, reagentB_list, reagentC_smallest, "ABc", nb_cpu, output_path, step,  app, gui=window)
        # AbC
        CB_reaction_3R(job_name, reaction, reagentA_list, reagentB_smallest, selected_reagentC_list, "AbC", nb_cpu, output_path, step, app, gui=window)
        CB_reaction_3R(job_name, reaction, selected_reagentA_list, reagentB_smallest, reagentC_list, "AbC", nb_cpu, output_path, step,  app, gui=window)
        # aBC
        CB_reaction_3R(job_name, reaction, reagentA_smallest, reagentB_list, selected_reagentC_list, "aBC", nb_cpu, output_path, step,  app, gui=window)
        CB_reaction_3R(job_name, reaction, reagentA_smallest, selected_reagentB_list, reagentC_list, "aBC", nb_cpu, output_path, step,  app, gui=window)
        # Merge files
        CB_merge_libraries(job_name, "AB")
        CB_merge_libraries(job_name, "AC")
        CB_merge_libraries(job_name, "BC")

    if step == 3:
        # Starting from AB
        CB_reaction_3R(job_name, reaction, reagentsAB_A_paths, reagentsAB_B_paths, reagentC_list, "AB", nb_cpu, output_path, step,  app, gui=window)
        # Starting from AC
        CB_reaction_3R(job_name, reaction, reagentsAC_A_paths, reagentB_list, reagentsAC_C_paths, "AC", nb_cpu, output_path, step,  app, gui=window)
        # Starting from BC
        CB_reaction_3R(job_name, reaction, reagentA_list, reagentsBC_B_paths, reagentsBC_C_paths, "BC", nb_cpu, output_path, step,  app, gui=window)
        # Merge files
        CB_merge_libraries(job_name, "ABC1")
        CB_merge_libraries(job_name, "ABC2")
        CB_merge_libraries(job_name, "ABC3")
        # Delete working directory
        os.rmdir(job_name + "_results/temp_combos")

    ####################################################################################################
    end_time = time.time()
    total_time = round(end_time - start_time, 0)
    if total_time >= 60:
        total_time_output = "{} minute(s) and {} second(s)".format(int(total_time//60), int(round(total_time%60, 0)))
    else:
        total_time_output = "{} second(s)".format(int(total_time))

    window.update_text("<div>\t</div>\n<b><div>Enumerations ended successfully in " + total_time_output + " !</div></b>\n\n")
    
    return "ok"

####################################################################################################

def CB_merge_libraries(job_name, tag):
    paths = glob.glob(job_name + "_results/output_libraries/" + tag + "*")
    with open(job_name + "_results/output_libraries/" + tag + ".smi", "w") as fileout:
        for path in paths:
            with open(path, "r") as filein:
                for line in filein.readlines():
                    fileout.write(line)
            os.remove(path)

####################################################################################################

def CB_reaction_3R(job_name, reaction, reagentsA, reagentsB, reagentsC, tag, nb_cpu, output_path, step, step_path, app = 0, gui = 0):
    # Var init
    products_output_names = []
    arguments_list = []

    # Output library path


    #######################################################
    # Generation of combinations of reagentsA vs reagentsB vs reagentsC

    ##### Step #1
    if step == 1:
        products_combinations = list(itertools.product(reagentsA, reagentsB, reagentsC))

        if tag == "Abc":
            tag_tosave = "A"
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + ".smi")
        if tag == "aBc":
            tag_tosave = "B"
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[1].split("/")[-1][:-4].split("group")[1] + ".smi")
        if tag == "abC":
            tag_tosave = "C"
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[2].split("/")[-1][:-4].split("group")[1] + ".smi")

    ##### Step #2
    if step == 2:
        products_combinations = list(itertools.product(reagentsA, reagentsB, reagentsC))
        if tag == "ABc":
            if reagentsA[0].split("/")[-1] == "selected_reagentsA.smi":
                tag_tosave = "AB1"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[1].split("/")[-1][:-4].split("group")[1] + ".smi")
            else:
                tag_tosave = "AB2"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + ".smi")
        if tag == "AbC":
            if reagentsA[0].split("/")[-1] == "selected_reagentsA.smi":
                tag_tosave = "AC1"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[2].split("/")[-1][:-4].split("group")[1] + ".smi")
            else:
                tag_tosave = "AC2"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + ".smi")
        if tag == "aBC":
            if reagentsB[0].split("/")[-1] == "selected_reagentsB.smi":
                tag_tosave = "BC1"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[2].split("/")[-1][:-4].split("group")[1] + ".smi")
            else:
                tag_tosave = "BC2"
                for combination in products_combinations:
                    products_output_names.append(tag_tosave + "_" + combination[1].split("/")[-1][:-4].split("group")[1] + ".smi")

    ##### Step #3
    if step == 3:
        products_combinations = []
        if tag == "AB":
            tag_tosave = "ABC1"
            for i in range(len(reagentsA)):
                products_combinations += list(itertools.product([reagentsA[i]], [reagentsB[i]], reagentsC))
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[0].split("/")[-1][:-4] + "_" + combination[2].split("/")[-1][:-4].split("group")[1] + ".smi")

        if tag == "AC":
            tag_tosave = "ABC2"
            for i in range(len(reagentsA)):
                products_combinations += list(itertools.product([reagentsA[i]], reagentsB, [reagentsC[i]]))
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[0].split("/")[-1][:-4] + "_" + combination[1].split("/")[-1][:-4].split("group")[1] + ".smi")

        if tag == "BC":
            tag_tosave = "ABC3"
            for i in range(len(reagentsB)):
                products_combinations += list(itertools.product(reagentsA, [reagentsB[i]], [reagentsC[i]]))
            for combination in products_combinations:
                products_output_names.append(tag_tosave + "_" + combination[1].split("/")[-1][:-4] + "_" + combination[0].split("/")[-1][:-4].split("group")[1] + ".smi")

        tag = "ABC" # Put back the proper tag


    # Generation of the big argument list to run python multiprocessing
    # output_library_path = job_name + "_results/output_libraries_" + tag
    output_library_path = job_name + "_results/output_libraries"
    os.makedirs(output_library_path, exist_ok = True)

    for i in range(len(products_output_names)):
        arguments_list.append((products_combinations[i][0], products_combinations[i][1], products_combinations[i][2], reaction, products_output_names[i], output_library_path, output_path, tag))

    # Logfile generation to benchmark, and see if there is some issues just in case 
    # with open(job_name + "_results/logfile.txt", "a") as logfile_in:
    #     logfile_in.write("Logfile for enumeration of " + tag + " products - Start : " + str(datetime.datetime.now()) + "\n")

    # Main multiprocessing loop
    with mp.Pool(nb_cpu) as pool:
        for result in pool.starmap(run_SATELLiTES_3reagents, arguments_list):
            if gui == 0:
                print(result)
            else:
                gui.update_text(result + "\n")

####################################################################################################

def CB_selection_step2_3R(reagentsA, reagentsB, reagentsC, job_name, step_path, gui = 0):
    if gui == 0:
        # # Tab completer
        # readline.set_completer_delims('\t')
        # readline.parse_and_bind("tab: complete")
        # readline.set_completer(path_completer)

        # step2_path_check = 0
        # while step2_path_check == 0:
        #     step2_filename = input(color.BOLD + color_ansi + name + " path (.smi file): " + color.END)
        #     if os.path.exists(step2_filename):
        #         if step2_filename[-4:] == ".smi":
        #             step2_path_check = 1
        #             cmd_copy = "cp " + step2_filename[0] + " " + job_name + "_results/."
        #             subprocess.run(cmd_copy, shell=True)
        #             output_path = ["./" + job_name + "_results/" + step2_filename[0].split("/")[-1]]
        #         else:
        #             print("\n    Not a .smi file")
        #     else:
        #         print("\n    File does not exist")

        exit("TO DO")

    else:
        # GUI mode ["./" + job_name + "_results/" + reagentsA.split("/")[-1].split(".")[0] + "_representative.smi"]
        # Errors to do, if does not exist (should not arrive) or in case of an header

        gui.update_text("3-reagents step #2 - Extraction of the best selected reagents\n")

        #####

        output_path_A = job_name + "_results/selected_reagentsA.smi"
        output_path_B = job_name + "_results/selected_reagentsB.smi"
        output_path_C = job_name + "_results/selected_reagentsC.smi"

        #####

        selected_reagentsA = []
        selected_reagentsB = []
        selected_reagentsC = []

        with open(step_path, "r") as filein:
            for line in filein.readlines():
                line_splitted = line.rstrip().split("\t")
                id_splitted = line_splitted[-1].split("-")
                if id_splitted[0] == "Abc":
                    selected_reagentsA.append(id_splitted[1])
                if id_splitted[0] == "aBc":
                    selected_reagentsB.append(id_splitted[2])
                if id_splitted[0] == "abC":
                    selected_reagentsC.append(id_splitted[3])

        # with open(pathB, "r") as filein:
        #     for line in filein.readlines():
        #         line_splitted = line.rstrip().split("\t")
        #         id_splitted = line_splitted[-1].split("-")
        #         if id_splitted[0] == "Abc":
        #             selected_reagentsA.append(id_splitted[1])
        #         if id_splitted[0] == "aBc":
        #             selected_reagentsB.append(id_splitted[2])
        #         if id_splitted[0] == "abC":
        #             selected_reagentsC.append(id_splitted[2])

        # with open(pathC, "r") as filein:
        #     for line in filein.readlines():
        #         line_splitted = line.rstrip().split("\t")
        #         id_splitted = line_splitted[-1].split("-")
        #         if id_splitted[0] == "Abc":
        #             selected_reagentsA.append(id_splitted[1])
        #         if id_splitted[0] == "aBc":
        #             selected_reagentsB.append(id_splitted[2])
        #         if id_splitted[0] == "abC":
        #             selected_reagentsC.append(id_splitted[2])

        selected_reagentsA_unique_iter = np.unique(selected_reagentsA).tolist()
        selected_reagentsB_unique_iter = np.unique(selected_reagentsB).tolist()
        selected_reagentsC_unique_iter = np.unique(selected_reagentsC).tolist()

        #####

        with open(reagentsA, "r") as file_input, open(output_path_A, "w") as file_output:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsA_unique_iter:
                    file_output.write(line)
                    selected_reagentsA_unique_iter.remove(id_temp)
                    if len(selected_reagentsA_unique_iter) == 0:
                        break
        
        gui.update_text("Best reagents A copied successfully in the results directory - selected_reagentsA.smi\n")

        #####

        with open(reagentsB, "r") as file_input, open(output_path_B, "w") as file_output:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsB_unique_iter:
                    file_output.write(line)
                    selected_reagentsB_unique_iter.remove(id_temp)
                    if len(selected_reagentsB_unique_iter) == 0:
                        break

        gui.update_text("Best reagents B copied successfully in the results directory - selected_reagentsB.smi\n")

        #####

        with open(reagentsC, "r") as file_input, open(output_path_C, "w") as file_output:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsC_unique_iter:
                    file_output.write(line)
                    selected_reagentsC_unique_iter.remove(id_temp)
                    if len(selected_reagentsC_unique_iter) == 0:
                        break

        gui.update_text("Best reagents C copied successfully in the results directory - selected_reagentsC.smi\n")

    return [output_path_A], [output_path_B], [output_path_C]

##################################################################################################################################################

def CB_selection_step3_3R(reagentsA, reagentsB, reagentsC, job_name, step_path, gui = 0):
    if gui == 0:
        # # Tab completer
        # readline.set_completer_delims('\t')
        # readline.parse_and_bind("tab: complete")
        # readline.set_completer(path_completer)

        # step2_path_check = 0
        # while step2_path_check == 0:
        #     step2_filename = input(color.BOLD + color_ansi + name + " path (.smi file): " + color.END)
        #     if os.path.exists(step2_filename):
        #         if step2_filename[-4:] == ".smi":
        #             step2_path_check = 1
        #             cmd_copy = "cp " + step2_filename[0] + " " + job_name + "_results/."
        #             subprocess.run(cmd_copy, shell=True)
        #             output_path = ["./" + job_name + "_results/" + step2_filename[0].split("/")[-1]]
        #         else:
        #             print("\n    Not a .smi file")
        #     else:
        #         print("\n    File does not exist")

        exit("TO DO")

    else:

        gui.update_text("3-reagents step #2 - Extraction of the best selected reagents\n")

        #####

        # output_path_A = job_name + "_results/selected_reagentsA.smi"
        # output_path_B = job_name + "_results/selected_reagentsB.smi"
        # output_path_C = job_name + "_results/selected_reagentsC.smi"

        #####

        selected_reagentsAB = []
        selected_reagentsAC = []
        selected_reagentsBC = []

        with open(step_path, "r") as filein:
            for line in filein.readlines():
                line_splitted = line.rstrip().split("\t")
                id_splitted = line_splitted[-1].split("-")
                if id_splitted[0] == "ABc":
                    selected_reagentsAB.append((id_splitted[1], id_splitted[2]))
                if id_splitted[0] == "AbC":
                    selected_reagentsAC.append((id_splitted[1], id_splitted[3]))
                if id_splitted[0] == "aBC":
                    selected_reagentsBC.append((id_splitted[2], id_splitted[3]))

        #####

        selected_reagentsAB_unique = np.unique(selected_reagentsAB, axis=0).tolist()
        selected_reagentsAC_unique = np.unique(selected_reagentsAC, axis=0).tolist()
        selected_reagentsBC_unique = np.unique(selected_reagentsBC, axis=0).tolist()

        #######################################################

        selected_reagentsAB_unique_A = [i[0] for i in selected_reagentsAB_unique]
        selected_reagentsAB_unique_A_iter = np.array(selected_reagentsAB_unique_A.copy())
        selected_reagentsAB_unique_B = [i[1] for i in selected_reagentsAB_unique]
        selected_reagentsAB_unique_B_iter = np.array(selected_reagentsAB_unique_B.copy())

        reagentsAB_A_paths = []
        reagentsAB_B_paths = []

        with open(reagentsA, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsAB_unique_A_iter:

                    temp_indexes = np.where(id_temp == np.array(selected_reagentsAB_unique_A))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/AB/A/{}.smi".format(temp_index + 1)
                        reagentsAB_A_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsAB_unique_A_iter = np.delete(selected_reagentsAB_unique_A_iter, np.where(id_temp == selected_reagentsAB_unique_A_iter))
                    if len(selected_reagentsAB_unique_A_iter) == 0:
                        break

        with open(reagentsB, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsAB_unique_B_iter:
                    temp_indexes = np.where(id_temp == np.array(selected_reagentsAB_unique_B))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/AB/B/{}.smi".format(temp_index + 1)
                        reagentsAB_B_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsAB_unique_B_iter = np.delete(selected_reagentsAB_unique_B_iter, np.where(id_temp == selected_reagentsAB_unique_B_iter))
                    if len(selected_reagentsAB_unique_B_iter) == 0:
                        break

        #######################################################

        selected_reagentsAC_unique_A = [i[0] for i in selected_reagentsAC_unique]
        selected_reagentsAC_unique_A_iter = np.array(selected_reagentsAC_unique_A.copy())
        selected_reagentsAC_unique_C = [i[1] for i in selected_reagentsAC_unique]
        selected_reagentsAC_unique_C_iter = np.array(selected_reagentsAC_unique_C.copy())

        reagentsAC_A_paths = []
        reagentsAC_C_paths = []

        with open(reagentsA, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsAC_unique_A_iter:

                    temp_indexes = np.where(id_temp == np.array(selected_reagentsAC_unique_A))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/AC/A/{}.smi".format(temp_index + 1)
                        reagentsAC_A_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsAC_unique_A_iter = np.delete(selected_reagentsAC_unique_A_iter, np.where(id_temp == selected_reagentsAC_unique_A_iter))
                    if len(selected_reagentsAC_unique_A_iter) == 0:
                        break

        with open(reagentsC, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsAC_unique_C_iter:
                    temp_indexes = np.where(id_temp == np.array(selected_reagentsAC_unique_C))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/AC/C/{}.smi".format(temp_index + 1)
                        reagentsAC_C_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsAC_unique_C_iter = np.delete(selected_reagentsAC_unique_C_iter, np.where(id_temp == selected_reagentsAC_unique_C_iter))
                    if len(selected_reagentsAC_unique_C_iter) == 0:
                        break

        #######################################################

        selected_reagentsBC_unique_B = [i[0] for i in selected_reagentsBC_unique]
        selected_reagentsBC_unique_B_iter = np.array(selected_reagentsBC_unique_B.copy())
        selected_reagentsBC_unique_C = [i[1] for i in selected_reagentsBC_unique]
        selected_reagentsBC_unique_C_iter = np.array(selected_reagentsBC_unique_C.copy())

        reagentsBC_B_paths = []
        reagentsBC_C_paths = []

        with open(reagentsB, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsBC_unique_B_iter:

                    temp_indexes = np.where(id_temp == np.array(selected_reagentsBC_unique_B))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/BC/B/{}.smi".format(temp_index + 1)
                        reagentsBC_B_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsBC_unique_B_iter = np.delete(selected_reagentsBC_unique_B_iter, np.where(id_temp == selected_reagentsBC_unique_B_iter))
                    if len(selected_reagentsBC_unique_B_iter) == 0:
                        break

        with open(reagentsC, "r") as file_input:
            for line in file_input:
                id_temp = line.rstrip().split("\t")[1]
                if id_temp in selected_reagentsBC_unique_C_iter:
                    temp_indexes = np.where(id_temp == np.array(selected_reagentsBC_unique_C))[0]
                    for temp_index in temp_indexes:
                        temp_path = job_name + "_results/temp_combos/BC/C/{}.smi".format(temp_index + 1)
                        reagentsBC_C_paths.append(temp_path)
                        with open(temp_path, "w") as file_output:
                            file_output.write(line)

                    selected_reagentsBC_unique_C_iter = np.delete(selected_reagentsBC_unique_C_iter, np.where(id_temp == selected_reagentsBC_unique_C_iter))
                    if len(selected_reagentsBC_unique_C_iter) == 0:
                        break

        #######################################################

    return reagentsAB_A_paths, reagentsAB_B_paths, reagentsAC_A_paths, reagentsAC_C_paths, reagentsBC_B_paths, reagentsBC_C_paths

##################################################################################################################################################


