# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
sys.path.append("/deps/KBBaseModules/")
import json
from os.path import exists
from ModelSEEDCommunity.mscommunitymodule import MSCommunityModule
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
    level=logging.INFO) 
#END_HEADER


class ModelSEEDCommunity:
    '''
    Module Name:
    ModelSEEDCommunity

    Module Description:
    A KBase module: ModelSEEDCommunity
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/ModelSEED/KB-ModelSEEDCommunity.git"
    GIT_COMMIT_HASH = "c92006ee9d33a83738ba50bf61cc010da036d4d8"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config["ATP_media_workspace"] = "94026"
        if "appdev" in self.config["workspace-url"]:
            self.config["ATP_media_workspace"] = "68393" 
        config["version"] = self.VERSION
        clients = {
            "Workspace":Workspace(self.config["workspace-url"], token=os.environ['KB_AUTH_TOKEN']),
            "KBaseReport":KBaseReport(os.environ['SDK_CALLBACK_URL'],token=os.environ['KB_AUTH_TOKEN']),
            "DataFileUtil":DataFileUtil(os.environ['SDK_CALLBACK_URL'],token=os.environ['KB_AUTH_TOKEN'])
        }
        self.msmodule = MSCommunityModule(config,"/kb/module",config['scratch'],os.environ['KB_AUTH_TOKEN'],clients,os.environ['SDK_CALLBACK_URL'])
        #END_CONSTRUCTOR
        pass


    def build_community_model(self, ctx, params):
        """
        This function combines multiple individual metabolic models into a community metabolic model
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN build_community_model
        output = self.msmodule.build_community_model(params)
        #END build_community_model

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method build_community_model return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def edit_update_community_model(self, ctx, params):
        """
        This function combines multiple individual metabolic models into a community metabolic model
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN edit_update_community_model
        output = self.msmodule.edit_update_community_model(params)
        #END edit_update_community_model

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method edit_update_community_model return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_community_flux_balance_analysis(self, ctx, params):
        """
        This function combines multiple individual metabolic models into a community metabolic model
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_community_flux_balance_analysis
        output = self.msmodule.run_community_flux_balance_analysis(params)
        #END run_community_flux_balance_analysis

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_community_flux_balance_analysis return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
