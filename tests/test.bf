models = {{"JC69","F81","TN93","GTR"}};


_model = 0;
_codon_options   = {{"No","Yes"}};    
_feature_options = {{"Neither", "HLA", "DRAM"}};
_feature_file = {{"","../data/hla.csv","../data/arv.csv"}};
_group_options = {{"No","Yes"}};

modelName = {4,1};

for (_do_codon = 0; _do_codon < 2; _do_codon += 1) {
    if (_do_codon) {
        modelName[0] = "codon";
    } else {
        modelName[0] = "nuc";
    }
    for (_do_feature = 0; _do_feature < 3; _do_feature += 1) {
        modelName [1] =  _feature_options [_do_feature];
        
        for (_do_group = 0; _do_group < 2; _do_group += 1) {
            if (_do_group) {
                modelName [2] = "group";
            } else {
                modelName [2] = "";
            }
            modelName[3] = models[_model];
            jmn = Join ("_", modelName);
            
            fprintf (stdout, "\n------------------------------------\n[TESTING `jmn`]\n------------------------------------\n");
            
            _options = {"00": _feature_options [_do_feature],
                                                     //"01": "../data/hla.csv",
                                                     "02": _codon_options [_do_codon],
                                                     "03": models[_model],
                                                     "04": "../data/clean_pol.fas",
                                                     "05": _group_options[_do_group],
                                                     //"06": "../data/group.csv",
                                                     //"07": "Dually Infected",
                                                     "08": "../results/`jmn`.csv",
                                                     "09": "../results/`jmn`.json"};
        
            if (Abs (_feature_file[_do_feature])) {
                _options ["01"] = _feature_file[_do_feature];
            }
            
            if (_do_group) {
                if (_do_feature == 2) {
                    _options ["06"] = "../data/group_treated.csv";
                    _options ["07"] = "Treated";
                
                } else {
                    _options ["06"] = "../data/group.csv";
                    _options ["07"] = "Dually Infected";
                }
            }
        
            ExecuteAFile ("../HBL/RateAnalysis.bf", _options);
            _model += 1;
            _model = _model % 4;
        }
    }
}


