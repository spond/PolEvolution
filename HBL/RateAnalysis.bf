
days_per_month = {{31,28,31,30,31,30,31,31,30,31,30,31}};

LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("LocalMGREV");
LoadFunctionLibrary ("CF3x4");
LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("dSdNTreeTools");

doHLA    = 0;
doDRM    = 0;
do_codon = 0;

ASSERTION_BEHAVIOR = 1;

ChoiceList (phenoClass,"Include HLA or DRM information",1,SKIP_NONE,
                "Neither", "No HLA or DRM information",
                "HLA", "Provide a .csv file with HLA haplotypes for each patient ID",
                "DRAM", "Provide a .csv file with ARV lists to determine relevant DRM for each patient ID"
            );
            
assert (phenoClass >= 0);
 

if (phenoClass == 1) {
    doHLA = 1;
    ExecuteAFile ("../CTL/CTLParser.bf");
    featureLabel = "HLA-targeted epitopes";
} else {
    if (phenoClass == 2) {
        doDRM = 1;
        ExecuteAFile ("../DRAM/DRAMParser.bf");
        featureLabel = "DRAM sites";
    }
}

ChoiceList (do_codon,"Perform codon analyses",1,SKIP_NONE,
                "No", "Only perform nucleotide rate analyses [FASTER]",
                "Yes", "Perform codon-based rate analyses (dN/dS) in addition to nucleotide rate analyses"
            );

assert (do_codon >= 0);

if (do_codon)  {
    dNdS_stensils = ComputeScalingStencils (0);
}
	
/* ________________________________________________________________________________________________*/


function InitializeDistancesFilter (filterID) {
	ExecuteCommands ("HarvestFrequencies (_dNucFreq,`filterID`,1,1,0)");
	
	_d_fR = _dNucFreq[0]+_dNucFreq[2];
	_d_fY = _dNucFreq[1]+_dNucFreq[3];
	
	if (_dNucFreq[0] == 0 || _dNucFreq[1] == 0 || _dNucFreq[2] == 0 || _dNucFreq[3] == 0) {
		_useK2P = 1;
	}
	else {
		_d_TN_K1 = 2*_dNucFreq[0]*_dNucFreq[2]/_d_fR;
		_d_TN_K2 = 2*_dNucFreq[1]*_dNucFreq[3]/_d_fY;
		_d_TN_K3 = 2*(_d_fR*_d_fY-_dNucFreq[0]*_dNucFreq[2]*_d_fY/_d_fR-_dNucFreq[1]*_dNucFreq[3]*_d_fR/_d_fY);
		_useK2P = 0;
	}
	
	
	summingVector = {4,1}["1"];
	return 0;
}

/* ________________________________________________________________________________________________*/

whichModel                    = "TN93";
doIG                          = 0;

ChoiceList (whichModel,"Use this nucleotide model",1,SKIP_NONE,
                "JC69", "Jukes Cantor 69 (all equal rates, all equal frequencies)",
                "F81", "Felsenstein 81 (all equal rates, empirical frequencies)",
                "TN93", "Tamura and Nei 93 (two transition rates, shared transversion rate, empirical frequencies)",
                "GTR", "General Time Reversible (6 substitution rates, empirical frequencies)"
                //,"GTR+V", "General Time Reversible with rate variation(6 substitution rates, empirical frequencies, General Discrete 3 rate distribution)"
            );
            
assert (whichModel >= 0);

whichModel = SELECTION_STRINGS;
if (whichModel == "GTR+V") {
    whichModel = "GTR";
    doIG = 1;
}

ACCEPT_ROOTED_TREES       = 1;
COUNT_GAPS_IN_FREQUENCIES = 0;
TRY_NUMERIC_SEQUENCE_MATCH = 1;

if (whichModel != "GTR") {
    global TRSV = 1;
    global TRST_CT = 1;    
    TN93_matrix = {{*,TRSV*t,t,TRSV*t}{TRSV*t,*,TRSV*t,TRST_CT*t}{t,TRSV*t,*,TRSV*t}{TRSV*t,TRST_CT*t,TRSV*t,*}};
    if (whichModel == "F81" || whichModel == "JC69") {
        TRSV := 1;
        TRST_CT := 1;
    }
} else {
    global nAC = 1;
    global nAT = 1;
    global nCG = 1;
    global nCT = 1;
    global nGT = 1;
    if (doIG) {
        GLOBAL_FPRINTF_REDIRECT = "/dev/null";
        LoadFunctionLibrary ("defineGamma", {"0":"General Discrete", "1": "3"});
        GLOBAL_FPRINTF_REDIRECT = "";
        TN93_matrix = {{*,nAC*t*c,t*c,nAT*t*c}{nAC*t*c,*,nCG*t*c,nCT*t*c}{t*c,nCG*t*c,*,nGT*t*c}{nAT*t*c,nCT*t*c,nGT*t*c,*}};
    } else {
        TN93_matrix = {{*,nAC*t,t,nAT*t}{nAC*t,*,nCG*t,nCT*t}{t,nCG*t,*,nGT*t}{nAT*t,nCT*t,nGT*t,*}};
    }
}


SetDialogPrompt ("Load the combined HIV polymerase sequence file");
DataSet       pol = ReadDataFile (PROMPT_FOR_FILE);


ChoiceList (groupAnalysis,"Split subjects into two groups based on external information (e.g. dual infection, treatment status)",1,SKIP_NONE,
                        "No", "Perform a joint analysis on all subjects",
                        "Yes", "Split the subjects into two groups, perform two-group analyses, and compare the rates between groups");
                        
assert (groupAnalysis >= 0);                        

allowed_ids = {};

if (groupAnalysis) {
    SetDialogPrompt ("Partitioning list of the form: pid, [0/1]:");
    id_list   = (ReadCSVTableText ("", 1))[1];
    byValue = {};
    for (k = 0; k < Abs (id_list); k+=1) {
        rec = id_list[k];
        allowed_ids[rec[0]] = 1 + rec[1];
        byValue[allowed_ids[rec[0]]] = 1;
    }
    
    ASSERTION_BEHAVIOR = 0;
    
    assert (Abs (byValue) > 1, "[RateAnalysis ERROR: all subjects belong to the same group]");
    
    fprintf (stdout, "\nDescriptive label for 'positive' patients, e.g. Treated:");
    fscanf (stdin, "String", groupLabel);
} 


SetDialogPrompt ("Save the resulting .CSV file to:");
fprintf   (PROMPT_FOR_FILE, CLEAR_FILE);
_masterCSVOutput = LAST_FILE_PATH;

SetDialogPrompt ("Save a .JSON file with divergence plots for all subjects (suitable for plotting using D3) to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
_masterJSONoutput = LAST_FILE_PATH;


fprintf (stdout, "\n\n");

DataSetFilter raw_filter = CreateFilter (pol, 1);
HarvestFrequencies(nuc_freqs,raw_filter,1,1,1);
HarvestFrequencies(nuc_3x4, raw_filter,3,1,1);
InitializeDistancesFilter ("raw_filter");


if (whichModel == "TN93") {
    AT := AC;
    CG := AC;
    GT := AC;
} else {
    if (whichModel == "F81" || whichModel == "JC69") {
        AT := 1;
        CG := 1;
        GT := 1;
        CT := 1;
        AC := 1;        
        if (whichModel == "JC69") {
            nuc_freqs = {4,1}["0.25"];
        }
    } 
}

Model TN93 = (TN93_matrix,nuc_freqs);

if (do_codon) {
    PopulateModelMatrix ("MG", nuc_3x4);
    codonFreqs = BuildCodonFrequencies (CF3x4 (nuc_3x4, GeneticCodeExclusions));
    Model MG94 = (MG, codonFreqs,0);
}

binByID = {};

for (seq_id = 0; seq_id < raw_filter.species; seq_id += 1) {
    GetString (seq_name, raw_filter, seq_id);
    id_date = splitOnRegExp (seq_name, "\\|");
    id      = id_date[0];
    date    = id_date[1];
    if (Abs(binByID [id]) == 0) {
        binByID [id] = {};
    }
    binByID [id] + {"id": seq_id, "date": date};
}

        
pids_with_data = {};

nucTreePrefix       = "T_";
nucFilterPrefix     = "nucData_";
nucRatePrefix       = "nucRate_";

doFeature = doHLA || doDRM;

featureAnalysisDoneByID = {};

if (doFeature) {
    nucTreePrefixAnnotated   = "Annotated_T_";
    nucFilterPrefixAnnotated = "Annotated_nucData_";
    nucRatePrefixAnnotated  = "Annotated_nucRate_";
}


maxD = 0; maxT = 0;
distancesByPID = {};
clockReject = {4,1};


fprintf (_masterCSVOutput, Join(",",{{"PID","SAMPLES","ML_NUCLEOTIDE","AVETN93","RSLVTN93", "p_mol_clock", "FOLLOWUP"}}));

if (groupAnalysis) {
    fprintf (_masterCSVOutput, ",", groupLabel);
}

if (do_codon) {
    fprintf (_masterCSVOutput, ",MG94_SYN,OMEGA,OMEGA_TEST_P");
}

if (doFeature) {
    fprintf (_masterCSVOutput, ",", Join(",",{{"ML_NUCLEOTIDE_FEATURE","FEATURE_RATE_FASTER","FEATURE_RATE_SLOWER"}}));
    if (do_codon) {
        fprintf (_masterCSVOutput, ",OMEGA_FEATURE");
    }
}

fprintf (_masterCSVOutput, "\n");

hla_diff_nuc = {2,1};
binByID ["handle_a_patient"][""];

if (doFeature) {
    fprintf (stdout, "\n##\t`featureLabel` nucleotide rate is significantly faster in ", hla_diff_nuc[0], " subjects");
    fprintf (stdout, "\n##\t`featureLabel` nucleotide rate is significantly slower in ", hla_diff_nuc[1], " subjects\n");
} else {
    fprintf (stdout, "\n");
}

if (do_codon) {
    if (doFeature) {
        fprintf (stdout, "##\tSubjects for whom dN/dS differs between `featureLabel` and NOT `featureLabel` = ", pos_selection_overall, "\n");   
    } else {
        fprintf (stdout, "##\tSubjects with overall positive selection = ", pos_selection_overall, "\n");
    }
}

/*
R_file = "../plot.R";
fprintf (R_file, CLEAR_FILE, KEEP_OPEN, 
"plot (0,0,type=\"l\",xaxp = c(0,",(maxT+1)$1,",",(maxT+1)$1,"), xlab=\"Time from baseline, years\",ylab=\"TN93 divergence from baseline, substitutions/site\",xlim=c(0,", maxT,
"), ylim = c (0,", Min(0.05,maxD), "));");
*/

fprintf (_masterJSONoutput, distancesByPID);

//distancesByPID ["plotter"][""];

runAGlobalModelFit (pids_with_data,doFeature);
if (do_codon) {
    runAGlobalModelFitCodon (pids_with_data,doFeature);
}

/*
maxTS = "" + (maxT+1);

fprintf (R_file,
"
abline(h=0.01,lty=\"dashed\");
polygon(x=c(0,`maxTS`,`maxTS`,0),y=c(0,",(global_fits["mi"])[0]*maxT,",",(global_fits["mi"])[2]*maxT,",0),col=rgb(0,0,0,0.25),lty=\"blank\");
lines(c(0,`maxTS`),c(0,",(global_fits["mi"])[1]*maxT,"),lwd=2);
legend(",maxT*0.6,",",Min(0.05,maxD)*0.8,",legend=c(\"Monoinfected\",\"Dually infected\"),fil=c(rgb(0,0,1,1),rgb(1,0,0,1)));
",
CLOSE_FILE);
*/

//fprintf (stdout, clockReject);

//------------------------------------------------------------------------------------

function plotter (key,value) {
    is_di = allowed_ids[key]>1;
    fprintf (R_file, "\n#",key,"\n\npoints (c(", Join(",",value[-1][0]),"),c(",Join(",",value[-1][1]),"),type = \"l\",col=rgb(",is_di,",0,",1-is_di,",0.75));");
    return 0;
}

//------------------------------------------------------------------------------------

function rateToBL (model,stencil) {
    if (Abs(stencil)) {
    	BRANCH_LENGTH_STENCIL = stencil;
    	pn = "synRate";
    	divider = 3;
    } else {
        ExecuteCommands ("GetString(pn,`model`,0)");    
    	divider = 1;
    }
    ExecuteCommands (pn + "=1");
    ExecuteCommands ("GetString(bl,`model`,-1)");
    bl = bl ^ {{"\\*c",""}};
    
    if (Abs(stencil)) {
    	BRANCH_LENGTH_STENCIL = 0;
    }
    
    return Eval (bl)/divider;
}

//------------------------------------------------------------------------------------
function reportTest (LA,LN,DF, oneSided) {
    LRT = 2*(LA-LN);
    p   = 1-CChi2 (LRT, DF);
    if (oneSided) { p = p / 2; }
    
    return " p-value = " + Format (p,6,3) + " [LRT = " + Format (LRT, 5, 2) + "]";
}   

//------------------------------------------------------------------------------------
function reportMLE (cmx, units) {
    return Format (cmx[1],10,4) + ", 95% CI [" + Format (cmx[0],8,3) + "," + Format (cmx[2],8,3) + "] " + units;
}

//------------------------------------------------------------------------------------

function runAGlobalModelFit (pids, useAnnotation) {

    SetParameter (STATUS_BAR_STATUS_STRING, "[RateAnalysis: Working on the joint nucleotide analysis]                                            ", 0);

    //LoadFunctionLibrary ("defineGamma", {"0": "General Discrete", "1" : "4"});
    
    if (useAnnotation) {
        global overall_nuc_rateHLA = 1;
        global overall_nuc_rateHLA = 1;    
    } 
    
    global overall_nuc_rate = 1;
    global overall_nuc_rate_di = 1;
    
    pat_count = Abs (pids);
    lkFuncComponents = {pat_count*2+2*Abs(featureAnalysisDoneByID),1};
    shift = 0;
    
    if (Abs (featureAnalysisDoneByID) == 0) {
        useAnnotation = 0;
    }
    
    if (useAnnotation) {
         
         for (_id = 0; _id < Abs(featureAnalysisDoneByID); _id += 1) {
            _pid = featureAnalysisDoneByID["INDEXORDER"][_id];
            lkFuncComponents [_id*2]   = nucFilterPrefixAnnotated + _pid;
            lkFuncComponents [_id*2+1] = nucTreePrefixAnnotated   + _pid;
            overall_nuc_rate += Eval (nucRatePrefixAnnotated + _pid);
            if (allowed_ids[_pid] > 1) {
                ExecuteCommands (nucRatePrefixAnnotated + _pid+ ":=overall_nuc_rate_diHLA");
            } else {        
                ExecuteCommands (nucRatePrefixAnnotated + _pid + ":=overall_nuc_rateHLA");
            }
            treeName = nucTreePrefixAnnotated + _pid;
            ExecuteCommands ("ReplicateConstraint (\"this1.?.?:="+nucRatePrefixAnnotated + _pid+"*this2.?.?__\",`treeName`,`treeName`)");
        }
        shift = 2*Abs(featureAnalysisDoneByID);   
    }
        
    for (_id = 0; _id < pat_count; _id += 1) {
        lkFuncComponents [shift + _id*2]   = nucFilterPrefix + pids[_id];
        lkFuncComponents [shift + _id*2+1] = nucTreePrefix   + pids[_id];
        overall_nuc_rate += Eval (nucRatePrefix + pids[_id]);
        if (allowed_ids[pids[_id]] > 1) {
            ExecuteCommands (nucRatePrefix + pids[_id] + ":=overall_nuc_rate_di");
        } else {        
            ExecuteCommands (nucRatePrefix + pids[_id] + ":=overall_nuc_rate");
        }
        treeName = nucTreePrefix + pids[_id];
        ExecuteCommands ("ReplicateConstraint (\"this1.?.?:="+nucRatePrefix + pids[_id]+"*this2.?.?__\",`treeName`,`treeName`)");
    }
    
    ExecuteCommands ("LikelihoodFunction nuc_joint = (" + Join (",", lkFuncComponents) + ");");
    Optimize (nuc_res, nuc_joint);
    
    COVARIANCE_PARAMETER = "overall_nuc_rate";
    COVARIANCE_PRECISION = 0.95;
    CovarianceMatrix (cmx,nuc_joint);
    COVARIANCE_PARAMETER = "overall_nuc_rate_di";
    CovarianceMatrix (cmxdi,nuc_joint);
    
    if (useAnnotation) {
        COVARIANCE_PARAMETER = "overall_nuc_rateHLA";
        COVARIANCE_PRECISION = 0.95;
        CovarianceMatrix (cmxHLA,nuc_joint);
        COVARIANCE_PARAMETER = "overall_nuc_rate_diHLA";
        CovarianceMatrix (cmxdiHLA,nuc_joint);       
    }

         
    conversion_factor = rateToBL("TN93",0);
    if (groupAnalysis == 0) {
        xtra = "";
        if (useAnnotation) {
            xtra = "[NOT `featureLabel`]";
        }
        fprintf (stdout, "\n##\tGlobal nucleotide substitution rate `xtra` = ", reportMLE (cmx*conversion_factor, "substitutions/site/year"), "\n");    
    } else {
        fprintf (stdout, "\n##\tGlobal nucleotide substitution rate [not `groupLabel`] = ", reportMLE (cmx*conversion_factor, "substitutions/site/year"), "\n");    
        fprintf (stdout, "##\tGlobal nucleotide substitution rate [`groupLabel`] = ", reportMLE (cmxdi*conversion_factor, "substitutions/site/year"), "\n");    
    }
    
    if (useAnnotation) {
        if (groupAnalysis == 0) {
             fprintf (stdout, "##\tGlobal nucleotide substitution rate [`featureLabel`] = ",  reportMLE (cmxHLA*conversion_factor, "substitutions/site/year"), "\n");          
        } else {
            fprintf (stdout, "##\tGlobal nucleotide substitution rate [`featureLabel` not `groupLabel`] = ",  reportMLE (cmxHLA*conversion_factor, "substitutions/site/year"), "\n");
            fprintf (stdout, "##\tGlobal nucleotide substitution rate [`featureLabel` `groupLabel`] = ",  reportMLE (cmxdiHLA*conversion_factor, "substitutions/site/year"), "\n");          
        }
        
        overall_nuc_rate    := overall_nuc_rateHLA;
        Optimize (nuc_res_constrained, nuc_joint);
        fprintf (stdout, "##\tGlobal nucleotide substitution rate is different between `featureLabel` and NOT `featureLabel` ", reportTest(nuc_res[1][0],nuc_res_constrained[1][0],1,0 ), "\n");
        overall_nuc_rate    = overall_nuc_rateHLA;
        if (groupAnalysis) {
            overall_nuc_rate_di := overall_nuc_rate_diHLA;
            Optimize (nuc_res_constrained, nuc_joint);
            fprintf (stdout, "##\tGlobal nucleotide substitution rate is different between `featureLabel` and NOT `featureLabel` in the `groupLabel` group, ", reportTest (nuc_res[1][0],nuc_res_constrained[1][0],1,0 ), "\n");
            overall_nuc_rate_di = overall_nuc_rate_diHLA;
        }
    } else {    
    
        if (groupAnalysis) {
            overall_nuc_rate := overall_nuc_rate_di;
            if (useAnnotation) {
                overall_nuc_rateHLA := overall_nuc_rate_diHLA;
            }
            Optimize (nuc_res_constrained, nuc_joint);

            fprintf (stdout, "##\tGlobal nucleotide substitution rate is different between `groupLabel` and NOT `groupLabel`, ", reportTest(nuc_res[1][0],nuc_res_constrained[1][0],1 + useAnnotation,0), "\n");
        }
    }
    
}

//------------------------------------------------------------------------------------

function runAGlobalModelFitCodon (pids, useAnnotation) {
    
    SetParameter (STATUS_BAR_STATUS_STRING, "[RateAnalysis: Working on the joint codon analysis]                                            ", 0);
    global overall_dNdS_rate    = 1;
    global overall_dNdS_rate_di = 1;
    
    pat_count = Abs (pids);
    lkFuncComponents = {pat_count*2+2*Abs(featureAnalysisDoneByID),1};
    shift = 0;
    
    if (Abs (featureAnalysisDoneByID) == 0) {
        useAnnotation = 0;
    }
    
    if (useAnnotation) {
         global overall_dNdS_rateHLA    = 1;
         global overall_dNdS_rate_diHLA = 1;
         for (_id = 0; _id < Abs(featureAnalysisDoneByID); _id += 1) {
            _pid = featureAnalysisDoneByID["INDEXORDER"][_id];
            codon_rate      = "codon_rate_" + _pid;
            omega           = "`codon_rate`_dNdSHLA";
            lkFuncComponents [_id*2]   = "codon_filter_" + _pid + "HLA";
            lkFuncComponents [_id*2+1] = nucTreePrefix   + _pid + "_codonHLA";
            if (allowed_ids[_pid] > 1) {
                ExecuteCommands ("`omega`:=overall_dNdS_rate_diHLA");
            } else {        
                ExecuteCommands ("`omega`:=overall_dNdS_rateHLA");
            }
            treeName = lkFuncComponents [_id*2+1];
        }
        shift = 2*Abs(featureAnalysisDoneByID);   
    } 
    
    for (_id = 0; _id < pat_count; _id += 1) {
        codon_rate      = "codon_rate_" + pids[_id];
        omega           = "`codon_rate`_dNdS";
        lkFuncComponents [shift + _id*2]   = "codon_filter_" + pids[_id];
        lkFuncComponents [shift + _id*2+1] = nucTreePrefix   + pids[_id] + "_codon";
        overall_nuc_rate += Eval (nucRatePrefix + pids[_id]);
        if (allowed_ids[pids[_id]] > 1) {
            ExecuteCommands ("`omega`:=overall_dNdS_rate_di");
        } else {        
            ExecuteCommands ("`omega`:=overall_dNdS_rate");
        }
    }
    
    
    ExecuteCommands ("LikelihoodFunction codon_join = (" + Join (",", lkFuncComponents) + ");");
 
    Optimize (codon_res, codon_join);

    COVARIANCE_PARAMETER = "overall_dNdS_rate";
    COVARIANCE_PRECISION = 0.95;
    CovarianceMatrix (cmx,codon_join);
    COVARIANCE_PARAMETER = "overall_dNdS_rate_di";
    CovarianceMatrix (cmxdi,codon_join);
    
    if (useAnnotation) {
        COVARIANCE_PARAMETER = "overall_dNdS_rateHLA";
        COVARIANCE_PRECISION = 0.95;
        CovarianceMatrix (cmxHLA,codon_join);
        COVARIANCE_PARAMETER = "overall_dNdS_rate_diHLA";
        CovarianceMatrix (cmxdiHLA,codon_join);       
    }
     
    if (groupAnalysis == 0) {
        fprintf (stdout, "\n##\tGlobal dN/dS ", reportMLE (cmx, ""), "\n");    
    } else {
        fprintf (stdout, "\n##\tGlobal dN/dS [not `groupLabel`] = ",  reportMLE (cmx, ""), "\n");
        fprintf (stdout, "##\tGlobal dN/dS  [`groupLabel`] = ",  reportMLE (cmxdi, ""), "\n");    
    }

    
    if (useAnnotation) {
        if (groupAnalysis == 0) {
             fprintf (stdout, "##\tGlobal dN/dS [`featureLabel`] = ",  reportMLE (cmxHLA, ""), "\n");    
        } else {
            fprintf (stdout, "##\tGlobal dN/dS [`featureLabel` not `groupLabel`] = ", reportMLE (cmxHLA, ""), "\n");    
            fprintf (stdout, "##\tGlobal dN/dS  [`featureLabel` `groupLabel`] = ", reportMLE (cmxdiHLA, ""), "\n");    
        }
        
    
        overall_dNdS_rate    := overall_dNdS_rateHLA;
        Optimize (codon_res_constrained, codon_join);
        fprintf (stdout, "##\tdN/dS is different between `featureLabel` and NOT `featureLabel`, ", reportTest (codon_res[1][0],codon_res_constrained[1][0],1,0), "\n");
        overall_dNdS_rate    = overall_dNdS_rateHLA;
        
        if (groupAnalysis) {
            overall_dNdS_rateHLA := overall_dNdS_rate_diHLA;
            Optimize (codon_res_constrained, codon_join);
            fprintf (stdout, "##\tdN/dS is different between `featureLabel` and NOT `featureLabel` in the `groupLabel` group, ", reportTest(codon_res[1][0], codon_res_constrained[1][0],1,0), "\n");
            overall_dNdS_rateHLA = overall_dNdS_rate_diHLA;
        }
    } else {
        if (groupAnalysis) {
            overall_dNdS_rate := overall_dNdS_rate_di;
            if (useAnnotation) {
                overall_dNdS_rateHLA := overall_dNdS_rate_diHLA;
            }
            Optimize (codon_res_constr, codon_join);

            fprintf (stdout, "##\tGlobal dN/dS is different between `groupLabel` and NOT `groupLabel`, ", reportTest(codon_res[1][0],codon_res_constr[1][0],1 + useAnnotation,0), "\n");
        }
    
    }
    
    return 0;
}

//------------------------------------------------------------------------------------

function handle_a_patient (key, value) {


    runFeatureStuff = 0;
    if (doFeature) {
        selector = feature_filter_by_pid [key];
        if (Abs (selector) == 0) {
            fprintf (stdout, "\n[RateAnalysis: Missing `featureLabel` data for subject `key`]\n", 0);
        }    
        else {
            featureAnalysisDoneByID [key] = 1;
            runFeatureStuff = 1;
        }
    }

    SetParameter (STATUS_BAR_STATUS_STRING, "[RateAnalysis: Working on subject `key` with " + Abs (value) + " samples]                                            ", 0);

    sample_count = Abs(value);
    if (allowed_ids[key] == 0 && groupAnalysis == 1) {
        return 0;
    }
    is_di = (allowed_ids[key]>1);
    
    if (sample_count == 1) {
        fprintf (MESSAGES_LOG, "[RateAnalysis: WARNING: single sequence for ", key, "; no analyses will be run for this individual]\n");
        return 0;
    }
    
    filter_string = {sample_count,1};
    dates_to_sort = {sample_count,3};
    
    for (seq_id = 0; seq_id < sample_count; seq_id += 1) {
        dates_to_sort[seq_id][0] = seq_id;
        dates_to_sort[seq_id][1] = numericDate((value[seq_id])["date"]);
    }
    

    pids_with_data + key;
    
    dates_to_sort = dates_to_sort % 1;
    
    for (seq_id = 1; seq_id < sample_count; seq_id += 1) {
       dates_to_sort[seq_id][2] = compute_date_difference((value[dates_to_sort[seq_id-1][0]])["date"], (value[dates_to_sort[seq_id][0]])["date"]);
    }
    
    for (seq_id = 0; seq_id < sample_count; seq_id += 1) {
        filter_string[sample_count-1-seq_id] = (value[dates_to_sort[seq_id][0]])["id"];
    }
    
    filter_string = Join(",",filter_string);
    
    filterName       = nucFilterPrefix + key;
    filterNameAnnotated    = nucFilterPrefixAnnotated + key;
    treeName         = nucTreePrefix + key;
    global_rate      = nucRatePrefix + key;
    
    if (runFeatureStuff) {
        treeNameHLA         = nucTreePrefixAnnotated + key;
        global_rateHLA      = nucRatePrefixAnnotated + key;
        ExecuteCommands      ("global `global_rateHLA` = 1");
    }
 
    codonFilterName  = "codon_filter_" + key;
    codon_rate      = "codon_rate_" + key;
    omega           = "`codon_rate`_dNdS";
    
    ExecuteCommands ("global `global_rate` = 1");
    ExecuteCommands ("global `codon_rate` = 1");
    ExecuteCommands ("global `codon_rate`_dNdS = 1");
    
    if (runFeatureStuff) {
        omegaHLA = omega+"HLA";
        ExecuteCommands ("global `omegaHLA` = 1");
    }
    
    if (runFeatureStuff) {
        selector = feature_filter_by_pid [key];
        ExecuteCommands ("DataSetFilter `filterName` = CreateFilter (pol, 1,selector[siteIndex]==0, filter_string)");
        ExecuteCommands ("DataSetFilter `filterNameAnnotated` = CreateFilter (pol, 1,selector[siteIndex]>0, filter_string)");
    } else {
        ExecuteCommands ("DataSetFilter `filterName` = CreateFilter (pol, 1,, filter_string)");
    }
    
    DataSetFilter computeDistanceFrom = CreateFilter (pol, 1, , filter_string);
    
    if (runFeatureStuff) {
        ExecuteCommands ("DataSetFilter `codonFilterName` = CreateFilter (`filterName`, 3,,, GeneticCodeExclusions)");
        ExecuteCommands ("DataSetFilter `codonFilterName`HLA = CreateFilter (`filterNameAnnotated`, 3,,, GeneticCodeExclusions)");
    } else {
        ExecuteCommands ("DataSetFilter `codonFilterName` = CreateFilter (pol, 3,, filter_string, GeneticCodeExclusions)");
    }
    
    UseModel        (TN93);
    treeString      = construct_tree_string (dates_to_sort[-1][2]);
    ExecuteCommands ("Tree `treeName` = `treeString`");
    ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rate`*this2.?.?__\",`treeName`,`treeName`)");
    if (runFeatureStuff) {
        ExecuteCommands ("Tree `treeNameHLA` = `treeString`");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rateHLA`*this2.?.?__\",`treeNameHLA`,`treeNameHLA`)");
        ExecuteCommands ("LikelihoodFunction localLF = (`filterName`,`treeName`, `filterNameAnnotated`, `treeNameHLA`);");
    } else {
        ExecuteCommands ("LikelihoodFunction localLF = (`filterName`,`treeName`);");
    }
    
    if (do_codon) {
        UseModel        (MG94);
        ExecuteCommands ("Tree `treeName`_codon = `treeString`");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.nonSynRate:=`omega`*this2.?.synRate\",`treeName`_codon,`treeName`_codon)");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=`codon_rate`*this2.?.t__\",`treeName`_codon,`treeName`)");
        if (runFeatureStuff) {
            ExecuteCommands ("Tree `treeName`_codonHLA = `treeString`");
            ExecuteCommands ("ReplicateConstraint (\"this1.?.nonSynRate:=`omegaHLA`*this2.?.synRate\",`treeName`_codonHLA,`treeName`_codonHLA)");
            ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=`codon_rate`*this2.?.t__\",`treeName`_codonHLA,`treeNameHLA`)");
            ExecuteCommands ("LikelihoodFunction localLF_codon = (`codonFilterName`,`treeName`_codon,`codonFilterName`HLA,`treeName`_codonHLA);");
        } else {
            ExecuteCommands ("LikelihoodFunction localLF_codon = (`codonFilterName`,`treeName`_codon);");
        }
    }
    
    UseModel        (TN93);

   
    if (whichModel == "TN93") {
        AC              = 1;
        CT              = 1;
        TRST_CT         = 1;
        TRSV            = 1;    
    } else {
        if (whichModel == "GTR") {
            AC = 1; AT = 1; CG = 1; CT = 1; GT = 1;
            nAC = 1; nAT = 1; nCG = 1; nCT = 1; nGT = 1;
            if (doIG) {
                PS_1=1/3;
                PS_2 =1/2;
                RS_1=0.5;
                RS_3=2;
            }
        }
    }
    
    if (runFeatureStuff) {
        Eval ("`global_rate` = 0.001");
        global hla_multiplier = 0.1;
        hla_multiplier :< 1;
        ExecuteCommands ("`global_rateHLA` := hla_multiplier*`global_rate`");
        Optimize (resHLA, localLF);
        Eval ("`global_rateHLA` = 0.001");
        ExecuteCommands ("`global_rate` := hla_multiplier*`global_rateHLA`");        
        Optimize (resHLA2, localLF);
        ExecuteCommands ("`global_rate` = `global_rateHLA`");
    }
    
    Optimize        (res,localLF);
    
    if (runFeatureStuff) {           
    
        p_hla = 0.5-0.5*CChi2 (2*(res[1][0]-resHLA[1][0]),1);
        if (p_hla < 0.05) {
            hla_diff_nuc[0] += 1;
        }
        p_hla2 = 0.5-0.5*CChi2 (2*(res[1][0]-resHLA2[1][0]),1);
        if (p_hla2 < 0.05) {
            hla_diff_nuc[1] += 1;
        }
    } 
    
    
    if (do_codon) {
        if (runFeatureStuff) {
            ExecuteCommands ("`omegaHLA` := hla_multiplier*`omega`");
            Optimize        (res_codon_null, localLF_codon);
            ExecuteCommands ("`omegaHLA` = `omega`");           
        } else {
            Eval ("`omega` = 0.1");
            ExecuteCommands ("`omega` :< 1");
            Optimize        (res_codon_null, localLF_codon);
            ExecuteCommands ("`omega` :< 10000");
        }
        
        USE_LAST_RESULTS = 1;
        Optimize        (res_codon, localLF_codon);
        USE_LAST_RESULTS = 0;
        
        pps = 0.5-0.5*CChi2(2*(res_codon[1][0]-res_codon_null[1][0]),1);
        if (pps <= 0.05) {
            pos_selection_overall += 1;
        }
    }
    
    time = 0;
 
    conversion_factor   = rateToBL("TN93",0);
    if (do_codon) {
        conversion_factorMG = rateToBL("MG94",dNdS_stensils["Syn"]);
    }
    
    
    if (sample_count > 2) {
        ExecuteCommands ("`global_rate`:=`global_rate`__");
        ExecuteCommands ("Tree `treeName`_local = " + construct_tree_string (dates_to_sort[-1][2]));
        ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rate`*this2.?.?__\",`treeName`_local,`treeName`_local)");
        bNames = Eval   ("BranchName(`treeName`,-1)");
        
        df   = 0;
        time = 0;
       
        for (bID = Columns (bNames)-2; bID >= 0 ; bID += -1) {
            if ((bNames[bID]$"^[0-9]+$")[0] == 0) {
                continue;
            }
            ExecuteCommands (treeName + "_local." + bNames[bID] + ".t= " + treeName_local + "." + bNames[bID] + ".t");
            df += 1;
            time += dates_to_sort[df][2]/2;
        }
        clockReject [is_di] +=1;
        
        ExecuteCommands ("DataSetFilter `filterName`_local = CreateFilter (pol, 1,, filter_string)");
        ExecuteCommands ("LikelihoodFunction localLF2 = (`filterName`_local,`treeName`_local);");
        
        Optimize (res_local, localLF2);
        idx = 1;
        div = 0;
        distances = {sample_count, 2 + runFeatureStuff};
        time = 0;
        
        for (bID = Columns (bNames)-2; bID > 0 ; bID += -1) {
            if ((bNames[bID]$"^[0-9]+$")[0] == 0) {
                continue;
            }
            time += dates_to_sort[idx][2];
            bll   = Eval("BranchLength (`treeName`_local,bID)");
            div  += bll;
            distances [idx][0] = time/365;
            distances [idx][1] = div;
            
            idx += 1;
        }
        
        p = 1-CChi2(2*(res_local[1][0] - res[1][0]), 1 );
        
    } else {
        distances = {1, 2};
        distances [idx][0] = dates_to_sort[1][2]/365;
        distances [idx][1] = +Eval("BranchLength (`treeName`,0)");
        p = 1;
    }
    if (p <= 0.05) {
        clockReject[2 + is_di] += 1;
    }
    distancesByPID [key] = distances;
    
        
    maxT = Max(maxT,Max(distances[-1][0],0));
    maxD = Max(maxD,Max(distances[-1][1],0));
    
    
    GetDataInfo (siteDifferenceCount, computeDistanceFrom, 0, sample_count-1, AVERAGE_AMBIGUITIES);
    GetDataInfo (siteDifferenceCountR, computeDistanceFrom, 0, sample_count-1, RESOLVE_AMBIGUITIES);
    followup = +dates_to_sort[-1][2];
    

    fprintf (_masterCSVOutput, key, ",", sample_count, ",", Format(Eval (global_rate + "*conversion_factor"),12,6), ",", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCount)/followup*365,12,6), ",", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCountR)/followup*365,12,6), ",", 
              Format(p,8,4), ",", followup);
        
    if (groupAnalysis) {
        fprintf (_masterCSVOutput, ",", is_di);
    }    
        
    if (do_codon) {
        fprintf (_masterCSVOutput, ",", Eval("`codon_rate`*conversion_factorMG"), ",", Format(Eval (omega),5,2), ",", Format(pps,4,2));
    }
    
    if (runFeatureStuff) {
        fprintf (_masterCSVOutput, ",", Format(Eval (global_rateHLA + "*conversion_factor"),12,6), ",", Format(p_hla,8,4),",", Format(p_hla2,8,4));
        if (do_codon) {
            fprintf (_masterCSVOutput, ",", Format(omegaHLA,5,2));
        }        
    }
    fprintf (_masterCSVOutput, "\n");
    
    return 0;
}

function construct_tree_string (diffs) {
    diffs = diffs * (1/365);
    seq_count = Rows(diffs);
    if (seq_count == 2) {
        return "(1:"+diffs[1]+",2:0)";
    }
    
    baseString = "((1:0):"+diffs[1]+",2:0)";
    for (_seq = 2; _seq < seq_count; _seq += 1) {
        baseString = "(" + baseString + ":" + diffs[_seq] + "," + (_seq+1 ) + ":0)";
    }
    //fprintf (stdout, baseString, "\n");
    
    return baseString;
}

function numericDate (date1) {
    yyyymmdd = splitOnRegExp (date1,"-");
    return 10000*(0 + yyyymmdd[0]) + 100*(0 + yyyymmdd[1]) + yyyymmdd[2];
}

function compute_date_difference (date1, date2) {
    yyyymmdd1 = splitOnRegExp (date1,"-");
    yyyymmdd2 = splitOnRegExp (date2,"-");
    
    date_array1 = {{0 + yyyymmdd1[0], 0 + yyyymmdd1[1], 0 + yyyymmdd1[2]}};
    date_array2 = {{0 + yyyymmdd2[0], 0 + yyyymmdd2[1], 0 + yyyymmdd2[2]}};

    d = 0;
    for (i = 0; i < 3; i += 1) {
        if (date_array1 [i] < date_array2[i]) {
            d = -1;
            break;
        }
        if (date_array1 [i] > date_array2[i]) {
            d = 1;
            break;
        }
    }
    
    if (d == 0) {
        return 0;
    }
    
    if (d == 1) {
        d = date_array2;
        date_array2 = date_array1;
        date_array1 = d;
    }
    
    diff = 0;
    
    for (y = date_array1[0] + 1; y < date_array2[0]; y += 1) {
        diff += 365 + is_leap_year (y);
    }
    
    diff += up_to_date (date_array2) + 365*(date_array1[0]!=date_array2[0]) - up_to_date(date_array1) + is_leap_year(date_array1[0]);
    
    return diff;
}

function up_to_date (ymd) {
    days = 0;
    if (ymd[1] > 1) {
        days += +days_per_month[{{0,0}}][{{0,ymd[1]-2}}];
    }  
    days += ymd[2];
    if (ymd[2] >= 3) {
        days += is_leap_year(ymd[0]);
    }
    return days;
}

function is_leap_year (year) {
    if (year % 4 == 0) {
        if (year % 100 == 0) {
            return (year%400 == 0);
        }
        return 1;
    }
    return 0;
}

