#!/usr/bin/env python
# -*- tab-width:4;indent-tabs-mode:f;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vim: set ts=4 et sw=4:

# Written by Seth Vigneron, Shoichet Lab UCSF, 2021
# Edited by Andrii Kyrylchuk, 2022-2023
# Takes input files of smiles for each reagent Building Block (BB) type
# Loops through every combination of BBs (one from each list)
# Runs all the Reaction SMARTs that should be used for each specific combination of BBs as defined in the .json input file
# Builds the 2D structure of the products to be formed
# Checks to see if the product passes disqualification checks based on the disqualification critera defined in the .json input file
# If the product passes, writes it to the output file and then moves onto the next set of BBs
from ast import arg
import os
from random import randrange
import sys
import json
import itertools, math
import Bespoke_Disqualification_Criterea_AK
from multiprocessing import Process, Pool, Lock
from rdkit import Chem
from rdkit.Chem import AllChem


#####################################################################################################################
#####################################################################################################################

# Class Objects for Product Molecules, and Reactions
#####################################################################################################################


#######################################

class Product_Identifiers:

    # def __init__(self, mol, smiles, purch, tranche, BB_IDs, BB_smiles, code_descriptor):
    def __init__(self, mol, smiles, tranche, BB_IDs, BB_smiles, code_descriptor):
        self.mol = mol
        self.smiles = smiles
        # self.purch = purch
        self.tranche = tranche
        self.BB_IDs = BB_IDs
        # print("SELF BB_IDS", self.BB_IDs)
        self.BB_smiles = BB_smiles
        self.code_descriptor = code_descriptor
        # self.inchi = Chem.MolToInchi(mol)

    def BB_component_IDs(self):
        components_ID = ''
        for bb in self.BB_IDs:
            components_ID += bb
            components_ID += '_'
            # components_ID += bb
        return components_ID 

    def BB_component_Smiles(self):
        components_smiles = ''
        for bb in self.BB_smiles:
            components_smiles += bb
            components_smiles += '\t'
        return components_smiles


#######################################
# Taken from the Breakdown script
class Reaction:
    
    def __init__(self, Name, ID, Synthons, breakdownSMARTS, breakdownRxn, buildSMARTS, buildRxn, numReactants):
        self.Name = Name
        self.ID = ID
        self.Synthons = Synthons
        self.breakdownSMARTS = breakdownSMARTS
        self.breakdownRxn = breakdownRxn
        self.buildSMARTS = buildSMARTS
        self.buildRxn = buildRxn
        self.numReactants = numReactants

#######################################
class Restrictions_Identifiers:

    def __init__(self, mol, BB_combo):
        self.mol = mol
        self.BB_combo = BB_combo


#####################################################################################################################
#####################################################################################################################

# Functions
#####################################################################################################################

def read_BBs(bb_file, all_synthons_dict):
    bb_fname = os.path.basename(bb_file)
    bb_fname = bb_fname.replace(" ", "-")
    # Example: bb_50_synthon_30_ArylHalide-Aryl\ iodides
    vendor_id = bb_fname.split("_synthon_")[0]
    synthon_id, synthon_name = bb_fname.split("_synthon_")[-1].split("_")
    print("Processing " + vendor_id, synthon_id, synthon_name)
    bb_f = open(bb_file, 'r')
    all_synthons_dict[synthon_id] = []
    for line in bb_f.readlines():
        smiles, cid = line.strip("\n").split(" ")[:2]
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if not mol:
            continue
        bb_id = f"{vendor_id}_{cid}"
        bb_dict ={
                "bb_id" : bb_id,
                "smiles" : smiles,
                "mol" : mol,
                }
        all_synthons_dict[synthon_id].append(bb_dict)
    return all_synthons_dict
        
    
#######################################
# Get the built product molecules from a given combination of BBs 
def Build_Products_and_Outputs(BB_combo, Rxn, Bespoke_Lib_Map,  \
    parent_fp=None, tc=None):
        
    All_Products = Determine_and_Run_Reactions(BB_combo, Rxn, Bespoke_Lib_Map, \
        parent_fp, tc)
    if All_Products == None:
        return
    return All_Products

def process_data(data ,rxn, Bespoke_Lib_Map, result_file):
    product = Build_Products_and_Outputs(data, rxn, Bespoke_Lib_Map)
    out_f = open(result_file, "a")
    if product:
        if isinstance(product, bool):
            return
        Product_code = "{}_{}".format(product.code_descriptor, product.BB_component_IDs())[:-1]
        Output_line = "%s\t%s\t%s\t%s\n" % (Product_code, product.smiles, product.tranche, product.BB_component_Smiles())
        print(Output_line)
        # if product.code_descriptor in rxnid_prod_dict.keys():
        #     rxnid_prod_dict[product.code_descriptor].append(Output_line)
        # else:
        #     rxnid_prod_dict[product.code_descriptor] = [Output_line]
        out_f.write(Output_line)
        return product

    
    


#######################################
# Taken from the Breakdown script
def Extract_Rxns(rxnIDs_list, Analoging_Rxn_Map, rxn_list):
    rxns_list = []
    for rxn_ID in rxnIDs_list:
        if rxn_ID not in rxn_list:
            continue
        name = Analoging_Rxn_Map["Reactions"][rxn_ID]["Name"]
        synthonIDs = Analoging_Rxn_Map["Reactions"][rxn_ID]["SynthonIDs"]
        synthonIDs_list = synthonIDs.split(',')
        breakdown_SMARTS = Analoging_Rxn_Map["Reactions"][rxn_ID]["Breakdown_SMARTS"]
        breakdown_rxn = AllChem.ReactionFromSmarts(breakdown_SMARTS)
        build_SMARTS = Analoging_Rxn_Map["Reactions"][rxn_ID]["Build_SMARTS"]
        build_rxn = AllChem.ReactionFromSmarts(build_SMARTS)
        numReactants = Analoging_Rxn_Map["Reactions"][rxn_ID]["NumReactants"]
    
        rxns_list.append(Reaction(name, rxn_ID, synthonIDs_list, breakdown_SMARTS, breakdown_rxn, build_SMARTS, build_rxn, numReactants))
    print("RXN_LIST", rxns_list[0].__dict__)
    return rxns_list

#######################################
# Take the given combination of BBs and return a list of product molecules
def Determine_and_Run_Reactions(BB_combo, Rxn, Bespoke_Lib_Map, parent_fp=None, tc=None):
    # print(BB_combo) 
    BB_combo_rxn_tag = ''
    products_list = []
    
    product_BB_mol_list = []
    product_BB_smi_list = []
    product_BB_ids_list = []


    for BB in BB_combo:
        product_BB_mol_list.append(Chem.MolFromSmiles(BB['smiles']))
        product_BB_smi_list.append(BB['smiles'])
        product_BB_ids_list.append(BB['bb_id'])

    print(product_BB_mol_list)
    print(product_BB_smi_list)
    # print("++++++++++++++++++++++++++++Reaction", Rxn.__dict__)
    rxn = Rxn.buildRxn
    product_code_descriptor = Rxn.ID
    # print("PROD_BB_MOL", product_BB_mol_list)
    product = rxn.RunReactants(tuple(product_BB_mol_list))  # Run Rxn on reactants
    if len(product) == 0:
        return
    product_mol = product[0][0]
    try:
         # If product molecule passes disqualification criterea, collect its identifiers, save it to a list of products built from this set of BBs
        product_smiles = Chem.MolToSmiles(Chem.RemoveHs(product_mol))
    except Exception as err:
        return

    product_tranche = get_tranche(product_mol)
    product = Product_Identifiers(product_mol, product_smiles, product_tranche, product_BB_ids_list, product_BB_smi_list, product_code_descriptor)
    # Check if product molecule passes all the given disqualification criterea
    return Product_Char_and_Disqualification_Check(product, Bespoke_Lib_Map)

    


#######################################
# Check if Product is Disqualified Based on Critera Provided in the json file input
def Product_Char_and_Disqualification_Check(product_object, Bespoke_Lib_Map):
    restricted_obj = Restrictions_Identifiers(product_object.mol, product_object.BB_IDs)
    disqualification = False  # Default is the current product is not disqualified
    # For each restriction listed in the json file, run the associated restriction class function to check if the product molecule passess 
    
    for product_restriction in Bespoke_Lib_Map["Product_Restrictions_Map"].keys():
        restriction_function_call_string = product_restriction + "_check"
        current_restriction_function = getattr(Bespoke_Disqualification_Criterea_AK, restriction_function_call_string)
        disqualification = current_restriction_function(restricted_obj)
        if disqualification:
            print("==DISQUAL", product_restriction, Chem.MolToSmiles(Chem.RemoveHs(product_object.mol)))
            return disqualification
    return product_object # If product molecule passes all the criterea, return the disqualification = False

######################################



#######################################
# Get the product molecules tranche based on its Heavy Atom Count and LogP
def get_tranche(mol):
    Product_HAC = mol.GetNumHeavyAtoms()
    logP = Chem.Crippen.MolLogP(Chem.RemoveHs(mol))
    tranche_logP, sign = scale_logp_value(logP)

    if Product_HAC > 99:
        Product_HAC = 99

    tranche_code = "H%02d%s%03d" % (Product_HAC, sign, tranche_logP)

    return tranche_code

#######################################
# Round LogP to Tranche Format
def scale_logp_value(logp):
    if logp < -9.0:
        logp = -9.0
    elif logp > 9.0:
        logp = 9.0
    if logp < 0.0 or logp >= 5.0:
        logp = 100*int(logp)
    else:
        logp = 10*int(10*logp)

    sign = 'M' if logp < 0.0 else 'P'

    tranche_logP = abs(logp)

    return tranche_logP, sign



            

#####################################################################################################################
#####################################################################################################################


#####################################################################################################################


def main(args):



    # Get the Bespoke Library Map Json File
    # no purchasability checks should be provided
    param_json_file = open(args.enum,)
    Bespoke_Lib_Map = json.load(param_json_file)
    print("Bespoke_Lib_MAP", Bespoke_Lib_Map)
    Lib_ID = Bespoke_Lib_Map["Library_ID"] # Get the Bespoke Library ID

    # Get the reactions Json file
    analoging_reactions_json = open(args.rxnjson,)
    Analoging_Rxn_Map = json.load(analoging_reactions_json)
    #print("ANALOGING_RXN_MAP", Analoging_Rxn_Map)
    Analoging_Rxn_IDs = Analoging_Rxn_Map["Reaction_IDs"]
    #print("ANALOGING_RXN_IDS", Analoging_Rxn_IDs)
    # Json map should contain 2 keys:
    #       Analoging_Rxn_Map["Reaction_IDs"] gives all rxn IDs of analoging reactons stored in the Json map
    #       Analoging_Rxn_Map["Reactions"] gives a dictionary where the keys are the rxn IDs



    # Create a dictionary of SynthonID:RXN IDs
    Analoging_Rxn_ID_SynthonID_Map = {}
    for rxn in Analoging_Rxn_Map["Reactions"]:
        synthon_IDs = Analoging_Rxn_Map["Reactions"][rxn]["SynthonIDs"].replace(" ", "").split(",")
        for s_id in synthon_IDs:
            if Analoging_Rxn_ID_SynthonID_Map.get(s_id) == None:
                Analoging_Rxn_ID_SynthonID_Map[s_id] = [Analoging_Rxn_Map["Reactions"][rxn]["ID"]]
            else:
                Analoging_Rxn_ID_SynthonID_Map[s_id].append(Analoging_Rxn_Map["Reactions"][rxn]["ID"])
    #print(Analoging_Rxn_ID_SynthonID_Map)




    # For each of the disqualification criterea listed in json file, create a class variable that equals the specified cut off value
    for restriction_type in Bespoke_Lib_Map["Product_Restrictions_Map"].keys():
        setattr(Restrictions_Identifiers, restriction_type, Bespoke_Lib_Map["Product_Restrictions_Map"][restriction_type])
    
    # List of all synthon IDs with smiles, mol , and BB ID
    all_synthon_dict = {}

    # list of all applicable reactions
    rxn_list = []

    # Running through all synthons in breakdown_BB_filename
    for bb_file in args.input:
        read_BBs(bb_file, all_synthon_dict)
# Building Analoging reaction list, applicable for the current synthon list if a RXN Id was not specified
    if args.rxn == None:
        for rxn in Analoging_Rxn_Map["Reactions"]:
            synthon_IDs = Analoging_Rxn_Map["Reactions"][rxn]["SynthonIDs"].split(",")
            print(synthon_IDs)
            count = 0
            for s_id in synthon_IDs:
                if s_id in synth_ID_list:
                    count += 1
                    print(count)
            if count == len(synthon_IDs):
                rxn_list.append(Analoging_Rxn_Map["Reactions"][rxn]["ID"])
    else:
        rxn_list.append(Analoging_Rxn_Map["Reactions"][args.rxn]["ID"])
    print("==RXN LIST", rxn_list)

    print("+++++++++++++++++++++++++++++++++++++ Enumeration")


    # Create a list of Reaction Class objects for all analoging reactions in the Analoging Rxn Map
    Analoging_Rxns = Extract_Rxns(Analoging_Rxn_IDs, Analoging_Rxn_Map, rxn_list)
    print(Analoging_Rxns)
    print("ANALOGING_RXNS", Analoging_Rxns[0].Synthons)
   
    product_list = []
    # Loop through all reactions, then through all BBs at the corresponding positions
    all_products_list = []
    for curr_rxn in Analoging_Rxns:
        curr_list_synth = []
        for curr_synth_id in curr_rxn.Synthons:
            curr_synth_id = curr_synth_id.strip()
            print(curr_synth_id)
            curr_list_synth.append([info for info in all_synthon_dict[curr_synth_id]])
        # print("CURR_LIST_SYNTH", curr_list_synth)
        rxn_combos = itertools.product(*curr_list_synth)
        # Starting parallel execution!
        print("Start numeration...") 
        num_processes = args.proc if args.proc else 5
        chunk_size = 50000
        rxnid_prod_dict = {}
      
        out_file_path = f"{args.output}_{curr_rxn.ID}.out"
        out_f = open(out_file_path, "w")
        with Pool(processes=num_processes) as pool:
            results = pool.starmap(process_data, [(x, curr_rxn, Bespoke_Lib_Map, out_file_path) for x in rxn_combos], chunksize=chunk_size)
        out_f.close() 






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Enumerate 2 or more BB libraries",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", type=str, nargs='+', help="List of snput BB files. BB file name should be vendor_synthon_ID")
    parser.add_argument("-r", "--rxn", type=str, help="Reaction ID if needed")
    parser.add_argument("-o", "--output", type=str, help="Output files prefix name")
    parser.add_argument("-j", "--json", type=str, help="Json file containing in/ex rules for each synthon")
    parser.add_argument("-e", "--enum", type=str, help="Bespoke Library parameters json, including filters for analog molecules")
    parser.add_argument("-a", "--rxnjson", type=str, help="Json file with the list of reactions")
    parser.add_argument("-p", "--proc", type=int, help="Number of cores to use")


    args = parser.parse_args()
    main(args)


