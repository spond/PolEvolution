LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("GrabBag");

present_name_patterns = {"0" : "^([0-9a-z]+)\\|([0-9]{8})", 
                         "1" : "^([0-9a-z]+)\\|([0-9]{4})-([0-9]{2})-([0-9]{2})",
                         "2" : "^([0-9a-z]+)_([0-9][0-9]?)-([0-9][0-9]?)-([0-9]{4})"};

SetDialogPrompt ("Load the raw pol sequence file");
DataSet       pol = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter raw_filter = CreateFilter (pol, 1);

for (seq_id = 0; seq_id < raw_filter.species; seq_id += 1) {
    GetString (seq_name, raw_filter, seq_id);
    pattern_type = matchStringToSetOfPatterns (seq_name, present_name_patterns);
    if (pattern_type < 0) { continue; }
    match = seq_name $ present_name_patterns[pattern_type];
    if (pattern_type == 0) {
        new_name = seq_name[match[2]][match[3]] + "|" + 
                                         add_zero_if_needed(seq_name[match[4]+4][match[4] + 7]) + 
                                         "-" + add_zero_if_needed(seq_name[match[4]][match[4] + 1]) + 
                                         "-" + add_zero_if_needed(seq_name[match[4]+2][match[4] + 3]);
                                         
    } else {
        if (pattern_type == 1) {
            new_name = seq_name[match[2]][match[3]] + "|" + 
                                             add_zero_if_needed(seq_name[match[4]][match[5]]) + 
                                             "-" + add_zero_if_needed(seq_name[match[6]][match[7]]) + 
                                             "-" + add_zero_if_needed(seq_name[match[8]][match[9]]);
        
        } else {
            new_name = seq_name[match[2]][match[3]] + "|" + 
                                             add_zero_if_needed(seq_name[match[8]][match[9]]) + 
                                             "-" + add_zero_if_needed(seq_name[match[4]][match[5]]) + 
                                             "-" + add_zero_if_needed(seq_name[match[6]][match[7]]);
        }
    }
    SetParameter (raw_filter, seq_id, new_name);
}

fprintf ("../clean_pol.fas", CLEAR_FILE, raw_filter);
fprintf (stdout, "\n\n");

function add_zero_if_needed (s) {
    if (Abs(s) == 1) {
        return "0" + s;
    }
    return s;
}