####################################################################################################
# SATELLiTES - version 1.0.10
# Copyright (C) 2024 Corentin BEDART
####################################################################################################

#############################
#         Libraries         #
#############################

from rdkit import Chem
from rdkit.Chem import AllChem

import os
import itertools
import time

#############################
#       Main functions      #
#############################

def run_SATELLiTES_2reagents(input_BB_A, input_BB_B, reaction, input_output_name, output_library_path, output_path, tag):
    start_time = time.time()
    
    os.chdir(output_path)
    # Building blocks
    all_BBs_objects_list = []
    all_BBs_objects_list.append(extract_bbs(input_BB_A, "A")) #Reagents A with the tag A
    all_BBs_objects_list.append(extract_bbs(input_BB_B, "B")) #Reagents B with the tag B

    # Reaction
    reaction_rdkit = AllChem.ReactionFromSmarts(reaction)

    # Products
    Product_library = open(output_library_path + "/" + input_output_name, 'w')
    
    for BB_combo_objects in itertools.product(*all_BBs_objects_list):
        Build_Products_and_Outputs(BB_combo_objects, reaction_rdkit, Product_library, tag)


    Product_library.close()

    # Return statement
    end_time = time.time()

    return "Finished in {} seconds - {} vs {}".format(round(end_time - start_time, 2), input_BB_A.split("/")[-1].split(".")[0], input_BB_B.split("/")[-1].split(".")[0])

#############################

def run_SATELLiTES_3reagents(input_BB_A, input_BB_B, input_BB_C, reaction, input_output_name, output_library_path, output_path, tag):
    start_time = time.time()
    os.chdir(output_path)
    # Building blocks
    all_BBs_objects_list = []
    all_BBs_objects_list.append(extract_bbs(input_BB_A, "A")) #Reagents A with the tag A
    all_BBs_objects_list.append(extract_bbs(input_BB_B, "B")) #Reagents B with the tag B
    all_BBs_objects_list.append(extract_bbs(input_BB_C, "C")) #Reagents B with the tag C

    # Reaction
    reaction_rdkit = AllChem.ReactionFromSmarts(reaction)

    # Products
    Product_library = open(output_library_path + "/" + input_output_name, 'a')
    
    for BB_combo_objects in itertools.product(*all_BBs_objects_list):
        Build_Products_and_Outputs(BB_combo_objects, reaction_rdkit, Product_library, tag)

    Product_library.close()

    # Return statement
    end_time = time.time()

    return "Finished in {} seconds - {} vs {}".format(round(end_time - start_time, 2), input_BB_A.split("/")[-1].split(".")[0], input_BB_B.split("/")[-1].split(".")[0])

#############################
#      Building blocks      #
#############################

class BuildingBlock:
    def __init__(self, smiles, bbid, rxn_tag, mol):
        self.smiles = smiles
        self.bbid = bbid
        self.rxn_tag = rxn_tag
        self.mol = mol

#############################

def extract_bbs(file_name, type):
    reagent_bb_objects = []\
    
    with open(file_name, "r") as reagent_file:
        reagent_file_lines = reagent_file.readlines()

    for bb_line in reagent_file_lines:
        bb_char_list = bb_line.strip().split('\t')
        bb_smiles = bb_char_list[0]
        bb_ID = bb_char_list[1]
        bb_rxn_tag = type
        bb_mol = Chem.AddHs(Chem.MolFromSmiles(bb_smiles))

        reagent_bb_objects.append(BuildingBlock(bb_smiles, bb_ID, bb_rxn_tag, bb_mol))
    
    return reagent_bb_objects

#############################
#    Products and output    #
#############################

class Product_Identifiers:
    def __init__(self, mol, smiles, BB_IDs, BB_smiles):
        self.mol = mol
        self.smiles = smiles
        self.BB_IDs = BB_IDs
        self.BB_smiles = BB_smiles

    def BB_component_IDs(self):
        return "-".join(self.BB_IDs) 

#############################

def Determine_and_Run_Reactions(BB_combo, reaction_rdkit):
    products_list = []
    
    product_BB_mol_list = []
    product_BB_smi_list = []
    product_BB_ids_list = []

    for BB in BB_combo:
        product_BB_mol_list.append(BB.mol)
        product_BB_smi_list.append(BB.smiles)
        product_BB_ids_list.append(BB.bbid)

    try:
        product_mol = reaction_rdkit.RunReactants(tuple(product_BB_mol_list))[0][0]
        product_smiles = Chem.MolToSmiles(Chem.RemoveHs(product_mol))
    except:
        with open("errors.txt", "a") as error_file:
            error_file.write("error with {}\n".format(" - ".join(product_BB_ids_list)))
        print("Error with ", tuple(product_BB_smi_list))
        return 0

    products_list.append(Product_Identifiers(product_mol, product_smiles, product_BB_ids_list, product_BB_smi_list))

    return products_list

#############################

def Build_Products_and_Outputs(BB_combo, reaction_rdkit, Product_library, tag):
    All_Products = Determine_and_Run_Reactions(BB_combo, reaction_rdkit)

    if All_Products == 0:
        return 0
    
    for product in All_Products:
        Output_line = "{}\t{}-{}\n".format(product.smiles, tag, product.BB_component_IDs())
        # Output_line = "%s\t-%s\n" % (product.smiles, product.BB_component_IDs())
        Product_library.write(Output_line)

#############################

