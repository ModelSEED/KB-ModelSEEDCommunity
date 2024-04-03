from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import jinja2
import pandas as pd
from optlang.symbolics import Zero, add
import cobra
from cobrakbase.core.kbasefba import FBAModel
from modelseedpy import AnnotationOntology, MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection,MSModelReport
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from kbbasemodules.basemodelingmodule import BaseModelingModule
from modelseedpy.community.mscommunity import MSCommunity

logger = logging.getLogger(__name__)

class MSCommunityModule(BaseModelingModule):
    def __init__(self,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModelingModule.__init__(self,"ModelSEEDCommunity",config,module_dir=module_dir,working_dir=working_dir,token=token,clients=clients,callback=callback)
        self.core_template = None
        self.util = None
        self.gs_template = None
        self.version = "0.1.1.msr"
        self.module_dir = module_dir
        self.native_ontology = False
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
    def build_community_model(self,params):
        self.initialize_call("build_community_model",params,True)
        self.validate_args(params,["workspace","model_refs","fbamodel_output_id"],{
            "abundances":{},
            "mixed_bag_model":False,
            "return_model_object":False,
            "save_report_to_kbase":True,
            "save_model_to_kbase":True,
            "save_gapfilling_fba_to_kbase":True
        })
        mdlutl = {}
        output = {}
        self.build_dataframe_report()
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if params["return_model_object"]:
            output["model_obj"] = mdlutl
        return output
    
    def edit_update_community_model(self,params):
        self.initialize_call("edit_update_community_model",params,True)
        self.validate_args(params,["workspace","community_model_ref"],{
            "return_model_object":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "save_models_to_kbase":True,
            "save_gapfilling_fba_to_kbase":True
        }) 
        output = {}
        self.build_dataframe_report()
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = result_table.to_json()
        if params["return_model_object"]:
            output["model_obj"] = params["model_obj"]
        return output

    def run_community_fba(self,params):
        self.initialize_call("run_community_fba",params,True)
        self.validate_args(params,["workspace","community_model_ref","fba_output_id"],{
            "media_id":"KBaseMedia/AuxoMedia",
            "target_reaction":"bio1",
            "expression_ref":None,
            "expression_condition":None,
            "feature_ko_list":[],
		    "reaction_ko_list":[],
            "prfba":True,
            "objective_fraction":0.9,
            "expression_coef":-10,
            "min_probability":0.05,
            "kinetics_coef":400,
            "exchange_coef":1,
            "max_c_uptake":None,
            "max_n_uptake":None,
            "max_p_uptake":None,
            "max_s_uptake":None,
            "max_o_uptake":None,
            "predict_community_composition":False,
            "return_fba_object":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "save_fba_to_kbase":True,
            "abundances":None
        })
        #Getting model
        mdlutl = self.get_model(params["community_model_ref"])
        #Setting objective
        mdlutl.model.objective = params["target_reaction"]
        #Setting media
        media = self.get_media(params["media_id"],None)
        mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        #Setting elemental uptake constraint
        exchange_string = "None"
        elemental_hash = {}
        if params["max_c_uptake"]:
            elemental_hash["C"] = params["max_c_uptake"]
            exchange_string += "C: "+str(params["max_c_uptake"])+"; "
        if params["max_n_uptake"]:
            elemental_hash["N"] = params["max_n_uptake"]
            exchange_string += "N: "+str(params["max_n_uptake"])+"; "
        if params["max_p_uptake"]:
            elemental_hash["P"] = params["max_p_uptake"]
            exchange_string += "P: "+str(params["max_p_uptake"])+"; "
        if params["max_s_uptake"]:
            elemental_hash["S"] = params["max_s_uptake"]
            exchange_string += "S: "+str(params["max_s_uptake"])+"; "
        if params["max_o_uptake"]:
            elemental_hash["O"] = params["max_o_uptake"]
            exchange_string += "O: "+str(params["max_o_uptake"])+"; "
        if len(elemental_hash) > 0:
            mdlutl.pkgmgr.getpkg("ElementUptakePkg").build_package(elemental_hash)
        #Dealing with abundances
        abundances = None
        names = None
        if params["abundances"]:
            lines = params["abundances"].split("\n")
            names = []
            abundances = {}
            for line in lines:
                array = line.split("=")
                if len(array) == 2:
                    names.append(array[0])
                    abundances[array[0]] = float(array[1])
        
        #Adding commkinetic constraints
        attributes = self.mdlutl.get_attributes()
        if "community" not in attributes and not params["abundances"]:
                print("Must provide abundances if this is an old community model")
                return {}
        comm_model = MSCommunity(
            model=mdlutl,
            names=names,
            abundances=abundances
        )
        mdlutl.pkgmgr.getpkg("CommKineticPkg").build_package(float(params["kinetics_coef"]), comm_model)
        output = {}
        #Optimizing and constraining objective to fraction of optimum
        output["max_growth"] = mdlutl.model.slim_optimize()
        if "C" in elemental_hash:
            output["carbon_uptake"] = mdlutl.getpkg("ElementUptakePkg").variables["elements"]["C"].primal
        if str(output["max_growth"]) != "nan":
            return output
        mdlutl.pkgmgr.getpkg("ObjConstPkg").build_package(output["max_growth"] * float(params["objective_fraction"]), None)
        #Creating minimal probability objective
        coef = {}
        for rxn in comm_model.model.reactions:
            if "rxn" == rxn.id[0:3]:
                currrxn = comm_model.model.reactions.get_by_id(rxn.id)
                if bool(params["prfba"]):
                    coef.update(
                        {
                            currrxn.forward_variable: max(
                                float(params["min_probability"]), (1 - float(rxn.probability))
                            )
                        }
                    )
                    coef.update(
                        {
                            currrxn.reverse_variable: max(
                                float(params["min_probability"]), (1 - float(rxn.probability))
                            )
                        }
                    )
                else:
                    coef.update({currrxn.forward_variable: 1})
                    coef.update({currrxn.reverse_variable: 1})
            elif "EX_" == rxn.id[0:3]:
                currrxn = comm_model.model.reactions.get_by_id(rxn.id)
                coef.update({currrxn.forward_variable: float(params["exchange_coef"])})
                coef.update({currrxn.reverse_variable: float(params["exchange_coef"])})
        #Adding expression data to minimum probability objective
        """
        for clade in feature_entries:
            total = 0
            for ftr in feature_entries[clade]:
                total += feature_entries[clade][ftr][condition]
            for ftr in feature_entries[clade]:
                feature_entries[clade][ftr][condition] = feature_entries[clade][ftr][condition]/total
        for rxn in mdlutl.model.reactions:
            highest_exp = 0
            for gene in rxn.genes:
                array = gene.id.split("_")
                array.pop()
                clade = "_".join(array)
                if gene.id in feature_entries[clade] and condition in feature_entries[clade][gene.id]:
                    if feature_entries[clade][gene.id][condition] > highest_exp:
                        highest_exp = feature_entries[clade][gene.id][condition]
            if highest_exp > min_expression: 
                coef.update(
                    {
                        rxn.forward_variable: float(params["expression_coef"]) * highest_exp
                    }
                )
                coef.update(
                    {
                        rxn.reverse_variable: float(params["expression_coef"]) * highest_exp
                    }
                )"""
        #Setting the objective
        mdlutl.model.objective = mdlutl.model.problem.Objective(Zero, direction="min")
        mdlutl.model.objective.set_linear_coefficients(coef)    
        #Solving the LP
        solution = mdlutl.model.optimize()
        fba_obj = self.save_solution_as_fba(solution,mdlutl,media,params["fba_output_id"],workspace=params["workspace"],fbamodel_ref=params("fbamodel_id"))
        context = {
            "overview": {
                'Model ID':mdlutl.wsid,
                'Media ID':media.id,
                'Target reaction':params["target_reaction"],
                'Gene knockouts':"; ".join(params["feature_ko_list"]),
                'Reaction knockouts':"; ".join(params["reaction_ko_list"]),
                'Exchange limits':exchange_string,
                'Kinetics coefficient':params["kinetics_coef"],
                'Objective fraction':params["objective_fraction"],
                'prFBA':params["prfba"],
                'Max growth':fba_obj.objective_value
            },
            "reactions": [],
            "exchange": [],
            "interaction": []
        }
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.module_dir+"/data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml']))
        html = env.get_template("FBAReportTemplate.html").render(context)
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Saving the output
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if bool(params["return_fba_object"]):
            output["fba_obj"] = fba_obj
        return output

    def build_dataframe_report(self,table,model_objs=None):        
        context = {
            "initial_model":table.iloc[0]["Model"]
        }
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.module_dir+"/data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml']))
        html = env.get_template("ReportTemplate.html").render(context)
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        if model_objs:
            for model in model_objs:
                msmodrep = MSModelReport(model)
                msmodrep.build_report(self.working_dir+"/html/"+model.wsid+"-recon.html")
                msmodrep.build_multitab_report(self.working_dir+"/html/"+model.wsid+"-full.html") 
        print("Output dir:",self.working_dir+"/html/index.html")
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Creating data table file
        json_str = '{"data":'+table.to_json(orient='records')+'}'
        with open(self.working_dir+"/html/data.json", 'w') as f:
            f.write(json_str)