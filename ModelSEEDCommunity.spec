/*
A KBase module: ModelSEEDCommunity
*/

module ModelSEEDCommunity {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This function combines multiple individual metabolic models into a community metabolic model
    */
    funcdef build_community_model(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

	/*
        This function combines multiple individual metabolic models into a community metabolic model
    */
    funcdef edit_update_community_model(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

	/*
        This function combines multiple individual metabolic models into a community metabolic model
    */
    funcdef run_community_flux_balance_analysis(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;
};
