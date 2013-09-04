days_per_month = {{31,28,31,30,31,30,31,31,30,31,30,31}};

LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("LocalMGREV");
LoadFunctionLibrary ("CF3x4");
LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("dSdNTreeTools");
LoadFunctionLibrary ("ReadDelimitedFiles");

doHLA = 1;
do_codon = 1;

if (doHLA) {
    ExecuteAFile ("../CTL/LANLParser.bf");
}

dNdS_stensils = ComputeScalingStencils (0);
		
/* ________________________________________________________________________________________________*/


function InitializeDistancesFilter (filterID)
{
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

ACCEPT_ROOTED_TREES       = 1;
COUNT_GAPS_IN_FREQUENCIES = 0;
TRY_NUMERIC_SEQUENCE_MATCH = 1;

if (whichModel != "REV") {
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
        LoadFunctionLibrary ("defineGamma", {"0":"Gamma", "1": "4"});
        TN93_matrix = {{*,nAC*t*c,t*c,nAT*t*c}{nAC*t*c,*,nCG*t*c,nCT*t*c}{t*c,nCG*t*c,*,nGT*t*c}{nAT*t*c,nCT*t*c,nGT*t*c,*}};
    } else {
        TN93_matrix = {{*,nAC*t,t,nAT*t}{nAC*t,*,nCG*t,nCT*t}{t,nCG*t,*,nGT*t}{nAT*t,nCT*t,nGT*t,*}};
    }
}


SetDialogPrompt ("Load the  pol sequence file");
DataSet       pol = ReadDataFile (PROMPT_FOR_FILE);

SetDialogPrompt ("PID filter");
id_list   = (ReadCSVTable ("", 1))[1];

allowed_ids = {};
for (k = 0; k < Rows (id_list); k+=1) {
    allowed_ids["0" + id_list[k][0]] = 1 + id_list[k][1];
}

fprintf (stdout, "\n\n");

DataSetFilter raw_filter = CreateFilter (pol, 1);
HarvestFrequencies(nuc_freqs,raw_filter,1,1,1);
HarvestFrequencies(nuc_3x4, raw_filter,3,1,1);

InitializeDistancesFilter ("raw_filter");

PopulateModelMatrix ("MG", nuc_3x4);

if (whichModel == "TN93") {
    AT := AC;
    CG := AC;
    GT := AC;
} else {
    if (whichModel == "F81") {
        AT := 1;
        CG := 1;
        GT := 1;
        CT := 1;
        AC := 1;        
    }
}

codonFreqs = BuildCodonFrequencies (CF3x4 (nuc_3x4, GeneticCodeExclusions));

Model MG94 = (MG, codonFreqs,0);

if (whichModel == "JC69") {
    nuc_freqs = {4,1}["0.25"];
}

Model TN93 = (TN93_matrix,nuc_freqs);


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


TN93_vectors 
     = {"MI RESOLVE" : {},
        "MI AVERAGE" : {},
        "DI RESOLVE" : {},
        "DI AVERAGE" : {}};
        
pids_with_data = {};

nucTreePrefix = "T_";
nucFilterPrefix = "nucData_";
nucRatePrefix   = "nucRate_";

nucTreePrefixHLA = "hla_T_";
nucFilterPrefixHLA = "hla_nucData_";
nucRatePrefixHLA   = "hla_nucRate_";

fprintf (stdout, "\n");

maxD = 0; maxT = 0;
distancesByPID = {};
clockReject = {4,1};

if (doHLA) {
    fprintf (stdout, "\n", Join("\t",{{"PID","Samples","MLTN93_HLA","MLTN93_NOHLA","AVETN93","RSLVTN93", "p", "MG94_syn","FOLLOWUP", "DI", "P_HLA_FASTER", "P_HLA_SLOWER"}}));
} else{
    fprintf (stdout, "\n", Join("\t",{{"PID","Samples","MLTN93","AVETN93","RSLVTN93", "p", "MG94_syn","FOLLOWUP", "DI"}}));
}
if (do_codon) {
    fprintf (stdout, "\tdN/dS\tpps");
}
fprintf (stdout, "\n");

hla_diff_nuc = {2,1};
binByID ["handle_a_patient"][""];

if (doHLA) {
    fprintf (stdout, "\nHLA rate faster = ", hla_diff_nuc[0], "\n");
    fprintf (stdout, "\nHLA rate slower = ", hla_diff_nuc[1], "\n");
}

if (do_codon) {
    fprintf (stdout, "\nSamples with overall positive selection = ", pos_selection_overall, "\n");
}

//fprintf (stdout, maxD, "\n", maxT, "\n");

R_file = "../plot.R";
fprintf (R_file, CLEAR_FILE, KEEP_OPEN, 
"plot (0,0,type=\"l\",xaxp = c(0,",(maxT+1)$1,",",(maxT+1)$1,"), xlab=\"Time from baseline, years\",ylab=\"TN93 divergence from baseline, substitutions/site\",xlim=c(0,", maxT,
"), ylim = c (0,", Min(0.05,maxD), "));");


distancesByPID ["plotter"][""];

global_fits = runAGlobalModelFit (pids_with_data,doHLA);
runAGlobalModelFitCodon (pids_with_data);

maxTS = "" + (maxT+1);

fprintf (R_file,
"
abline(h=0.01,lty=\"dashed\");
polygon(x=c(0,`maxTS`,`maxTS`,0),y=c(0,",(global_fits["mi"])[0]*maxT,",",(global_fits["mi"])[2]*maxT,",0),col=rgb(0,0,0,0.25),lty=\"blank\");
lines(c(0,`maxTS`),c(0,",(global_fits["mi"])[1]*maxT,"),lwd=2);
legend(",maxT*0.6,",",Min(0.05,maxD)*0.8,",legend=c(\"Monoinfected\",\"Dually infected\"),fil=c(rgb(0,0,1,1),rgb(1,0,0,1)));
",
CLOSE_FILE);

fprintf (stdout, clockReject);

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

function runAGlobalModelFit (pids, useHLA) {

    LoadFunctionLibrary ("defineGamma", {"0": "General Discrete", "1" : "4"});
    
    if (useHLA) {
        global overall_nuc_rateHLA = 1;
        global overall_nuc_rateHLA = 1;    
    } 
    
    global overall_nuc_rate = 1;
    global overall_nuc_rate_di = 1;
    
    pat_count = Abs (pids);
    lkFuncComponents = {pat_count*(2+2*useHLA),1};
    shift = 0;
    
    if (useHLA) {
         for (_id = 0; _id < pat_count; _id += 1) {
            lkFuncComponents [_id*2]   = nucFilterPrefixHLA + pids[_id];
            lkFuncComponents [_id*2+1] = nucTreePrefixHLA   + pids[_id];
            overall_nuc_rate += Eval (nucRatePrefixHLA + pids[_id]);
            if (allowed_ids[pids[_id]] > 1) {
                ExecuteCommands (nucRatePrefixHLA + pids[_id] + ":=overall_nuc_rate_diHLA");
            } else {        
                ExecuteCommands (nucRatePrefixHLA + pids[_id] + ":=overall_nuc_rateHLA");
            }
            treeName = nucTreePrefixHLA + pids[_id];
            ExecuteCommands ("ReplicateConstraint (\"this1.?.?:="+nucRatePrefixHLA + pids[_id]+"*this2.?.?__\",`treeName`,`treeName`)");
        }
        shift = 2*pat_count;   
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
    
    overall_nuc_rate = overall_nuc_rate / pat_count * (1 + useHLA > 0);
    fprintf (stdout, "\nMean rate = ", overall_nuc_rate, "\n");
    ExecuteCommands ("LikelihoodFunction nuc_joint = (" + Join (",", lkFuncComponents) + ");");
    Optimize (nuc_res, nuc_joint);
    
    Export (lfExport, nuc_joint);
    filename = "../byPID/joint.lf";
    fprintf (filename,  CLEAR_FILE, lfExport);

    COVARIANCE_PARAMETER = "overall_nuc_rate";
    COVARIANCE_PRECISION = 0.95;
    CovarianceMatrix (cmx,nuc_joint);
    COVARIANCE_PARAMETER = "overall_nuc_rate_di";
    CovarianceMatrix (cmxdi,nuc_joint);
    
    if (useHLA) {
        COVARIANCE_PARAMETER = "overall_nuc_rateHLA";
        COVARIANCE_PRECISION = 0.95;
        CovarianceMatrix (cmxHLA,nuc_joint);
        COVARIANCE_PARAMETER = "overall_nuc_rate_diHLA";
        CovarianceMatrix (cmxdiHLA,nuc_joint);       
    }

    //GetInformation (cInfo, c);
     
    fprintf (stdout, "Log(L) = ", nuc_res[1][0], "\n");
    
    conversion_factor = rateToBL("TN93",0);
    fprintf (stdout, "\nMonoinfected rate = ", conversion_factor*overall_nuc_rate, " [", conversion_factor*cmx[0], " - ", conversion_factor*cmx[2], "]");
    fprintf (stdout, "\nDually infected rate = ", conversion_factor*overall_nuc_rate_di, " [", conversion_factor*cmxdi[0], " - ", conversion_factor*cmxdi[2], "]\n");
    
    if (useHLA) {
        fprintf (stdout, "\nMonoinfected rate (HLA) = ", conversion_factor*overall_nuc_rateHLA, " [", conversion_factor*cmxHLA[0], " - ", conversion_factor*cmxHLA[2], "]");
        fprintf (stdout, "\nDually infected rate (HLA) = ", conversion_factor*overall_nuc_rate_diHLA, " [", conversion_factor*cmxdiHLA[0], " - ", conversion_factor*cmxdiHLA[2], "]\n");
        
        overall_nuc_rate    := overall_nuc_rateHLA;
        //overall_nuc_rate_di := overall_nuc_rate_diHLA;
        Optimize (nuc_res_constrained, nuc_joint);
        fprintf (stdout, "HLA == non-HLA (mono): p = ", 1-CChi2(2*(nuc_res[1][0]-nuc_res_constrained[1][0]),1 ), "\n");
        overall_nuc_rate    = overall_nuc_rateHLA;
        overall_nuc_rate_di := overall_nuc_rate_diHLA;
        Optimize (nuc_res_constrained, nuc_joint);
        fprintf (stdout, "HLA == non-HLA (dual): p = ", 1-CChi2(2*(nuc_res[1][0]-nuc_res_constrained[1][0]),1 ), "\n");
        overall_nuc_rate_di = overall_nuc_rate_diHLA;
    }
    
    overall_nuc_rate := overall_nuc_rate_di;
    if (useHLA) {
        overall_nuc_rateHLA := overall_nuc_rate_diHLA;
    }
    
    Optimize (nuc_res_constrained, nuc_joint);
    //fprintf (stdout, nuc_joint, "\n");
    
    fprintf (stdout, "MI == DI: p = ", 1-CChi2(2*(nuc_res[1][0]-nuc_res_constrained[1][0]),1 + useHLA), "\n");
   
    return {"mi":cmx*conversion_factor,"di":cmxdi*conversion_factor};
}

//------------------------------------------------------------------------------------

function runAGlobalModelFitCodon (pids) {
    
    global overall_dNdS_rate    = 1;
    global overall_dNdS_rate_di = 1;
    global overall_dNdS_rateHLA    = 1;
    global overall_dNdS_rate_diHLA = 1;
    
    pat_count = Abs (pids);
    lkFuncComponents = {pat_count*4,1};
    shift = 0;
    
     for (_id = 0; _id < pat_count; _id += 1) {
        codon_rate      = "codon_rate_" + pids[_id];
        omega           = "`codon_rate`_dNdSHLA";
        lkFuncComponents [_id*2]   = "codon_filter_" + pids[_id] + "HLA";
        lkFuncComponents [_id*2+1] = nucTreePrefix   + pids[_id] + "_codonHLA";
        if (allowed_ids[pids[_id]] > 1) {
            ExecuteCommands ("`omega`:=overall_dNdS_rate_diHLA");
        } else {        
            ExecuteCommands ("`omega`:=overall_dNdS_rateHLA");
        }
        treeName = lkFuncComponents [_id*2+1];
    }
    shift = 2*pat_count;   

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
    Export (lfExport, codon_join);
    filename = "../byPID/joint_codon.lf";
    fprintf (filename,  CLEAR_FILE, lfExport);
    Optimize (codon_res, codon_join);
    //fprintf (stdout, codon_join, "\n");

    COVARIANCE_PARAMETER = "overall_dNdS_rate";
    COVARIANCE_PRECISION = 0.95;
    CovarianceMatrix (cmx,codon_join);
    COVARIANCE_PARAMETER = "overall_dNdS_rate_di";
    CovarianceMatrix (cmxdi,codon_join);
    
    COVARIANCE_PARAMETER = "overall_dNdS_rateHLA";
    COVARIANCE_PRECISION = 0.95;
    CovarianceMatrix (cmxHLA,codon_join);
    COVARIANCE_PARAMETER = "overall_dNdS_rate_diHLA";
    CovarianceMatrix (cmxdiHLA,codon_join);       
     
    fprintf (stdout, "Log(L) = ", codon_res[1][0], "\n");
    
    fprintf (stdout, "\nMonoinfected dNdS = ", overall_dNdS_rate, " [", cmx[0], " - ", cmx[2], "]");
    fprintf (stdout, "\nDually infected dNdS = ", overall_dNdS_rate_di, " [", cmxdi[0], " - ", cmxdi[2], "]\n");
    
    fprintf (stdout, "\nMonoinfected dNdS (HLA) = ", overall_dNdS_rateHLA, " [", cmxHLA[0], " - ", cmxHLA[2], "]");
    fprintf (stdout, "\nDually infected dNdS (HLA) = ", overall_dNdS_rate_diHLA, " [", cmxdiHLA[0], " - ", cmxdiHLA[2], "]\n");
    
    overall_dNdS_rate    := overall_dNdS_rateHLA;
    Optimize (codon_res_constrained, codon_join);
    fprintf (stdout, "HLA == non-HLA (mono): p = ", 1-CChi2(2*(codon_res[1][0]-codon_res_constrained[1][0]),1 ), "\n");
    overall_dNdS_rate    = overall_dNdS_rateHLA;
    overall_dNdS_rateHLA := overall_dNdS_rate_diHLA;
    Optimize (codon_res_constrained, codon_join);
    fprintf (stdout, "HLA == non-HLA (dual): p = ", 1-CChi2(2*(codon_res[1][0]-codon_res_constrained[1][0]),1 ), "\n");
    
    return 0;
}

//------------------------------------------------------------------------------------

function handle_a_patient (key, value) {
    sample_count = Abs(value);
    if (allowed_ids[key] == 0) {
        return 0;
    }
    is_di = (allowed_ids[key]>1);
    
    if (sample_count == 1) {
        fprintf (stdout, "WARNING: single sequence for ", key, "\n");
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
    //fprintf (stdout, filter_string, "\n");
    
    filterName       = nucFilterPrefix + key;
    filterNameHLA    = nucFilterPrefixHLA + key;
    treeName         = nucTreePrefix + key;
    global_rate      = nucRatePrefix + key;
    
    if (doHLA) {
        treeNameHLA         = nucTreePrefixHLA + key;
        global_rateHLA      = nucRatePrefixHLA + key;
        ExecuteCommands      ("global `global_rateHLA` = 1");
    }
 
    codonFilterName  = "codon_filter_" + key;
    codon_rate      = "codon_rate_" + key;
    omega           = "`codon_rate`_dNdS";
    
    ExecuteCommands ("global `global_rate` = 1");
    ExecuteCommands ("global `codon_rate` = 1");
    ExecuteCommands ("global `codon_rate`_dNdS = 1");
    
    if (doHLA) {
        omegaHLA = omega+"HLA";
        ExecuteCommands ("global `omegaHLA` = 1");
    }
    
    if (doHLA) {
        selector = hla_filter_by_pid [key];
        assert (Abs (selector) > 0, "Missing HLA information for " +  key);    
        //fprintf (stdout, key, ":", selector, "\n");
        ExecuteCommands ("DataSetFilter `filterName` = CreateFilter (pol, 1,selector[siteIndex]==0, filter_string)");
        ExecuteCommands ("DataSetFilter `filterNameHLA` = CreateFilter (pol, 1,selector[siteIndex]>0, filter_string)");
    } else {
        ExecuteCommands ("DataSetFilter `filterName` = CreateFilter (pol, 1,, filter_string)");
    }
    
    DataSetFilter computeDistanceFrom = CreateFilter (pol, 1, , filter_string);
    
    if (doHLA) {
        ExecuteCommands ("DataSetFilter `codonFilterName` = CreateFilter (`filterName`, 3,,, GeneticCodeExclusions)");
        ExecuteCommands ("DataSetFilter `codonFilterName`HLA = CreateFilter (`filterNameHLA`, 3,,, GeneticCodeExclusions)");
    } else {
        ExecuteCommands ("DataSetFilter `codonFilterName` = CreateFilter (pol, 3,, filter_string, GeneticCodeExclusions)");
    }
    
    UseModel        (TN93);
    treeString      = construct_tree_string (dates_to_sort[-1][2]);
    ExecuteCommands ("Tree `treeName` = `treeString`");
    ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rate`*this2.?.?__\",`treeName`,`treeName`)");
    if (doHLA) {
        ExecuteCommands ("Tree `treeNameHLA` = `treeString`");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rateHLA`*this2.?.?__\",`treeNameHLA`,`treeNameHLA`)");
        ExecuteCommands ("LikelihoodFunction localLF = (`filterName`,`treeName`, `filterNameHLA`, `treeNameHLA`);");
    } else {
        ExecuteCommands ("LikelihoodFunction localLF = (`filterName`,`treeName`);");
    }
    
    if (do_codon) {
        UseModel        (MG94);
        ExecuteCommands ("Tree `treeName`_codon = `treeString`");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.nonSynRate:=`omega`*this2.?.synRate\",`treeName`_codon,`treeName`_codon)");
        ExecuteCommands ("ReplicateConstraint (\"this1.?.synRate:=`codon_rate`*this2.?.t__\",`treeName`_codon,`treeName`)");
        if (doHLA) {
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
        if (whichModel == "REV") {
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
    
    if (doHLA) {
        Eval ("`global_rate` = 0.001");
        global hla_multiplier = 0.1;
        hla_multiplier :< 1;
        ExecuteCommands ("`global_rateHLA` := hla_multiplier*`global_rate`");
        Optimize (resHLA, localLF);
        //fprintf (stdout, Eval (global_rate), ":", Eval (global_rateHLA), ":", resHLA2[1][0], "\n");
        Eval ("`global_rateHLA` = 0.001");
        ExecuteCommands ("`global_rate` := hla_multiplier*`global_rateHLA`");        
        Optimize (resHLA2, localLF);
        //fprintf (stdout, Eval (global_rate), ":", Eval (global_rateHLA), ":", resHLA[1][0], "\n");
        ExecuteCommands ("`global_rate` = `global_rateHLA`");
    }
    
    Optimize        (res,localLF);
    
    
    if (doHLA) {           
    
        p_hla = 0.5-0.5*CChi2 (2*(res[1][0]-resHLA[1][0]),1);
        if (p_hla < 0.05) {
            hla_diff_nuc[0] += 1;
        }
        p_hla2 = 0.5-0.5*CChi2 (2*(res[1][0]-resHLA2[1][0]),1);
        if (p_hla2 < 0.05) {
            hla_diff_nuc[1] += 1;
        }
    } 
    
    //fprintf (stdout, Eval (global_rate), ":", Eval (global_rateHLA), "\n");
   
    Export (lfExport, localLF);
    filename = "../byPID/" + key + ".lf";
    fprintf (filename,  CLEAR_FILE, lfExport);

    if (do_codon) {
        if (doHLA) {
            ExecuteCommands ("`omegaHLA` := hla_multiplier*`omega`");
            Optimize        (res_codon_null, localLF_codon);
            ExecuteCommands ("`omegaHLA` = `omega`");
            
        } else {
            Eval ("`omega` = 0.1");
            ExecuteCommands ("`omega` :< 1");
            Optimize        (res_codon_null, localLF_codon);
            ExecuteCommands ("`omega` :< 10000");
        }
        
        Optimize        (res_codon, localLF_codon);
        pps = 0.5-0.5*CChi2(2*(res_codon[1][0]-res_codon_null[1][0]),1);
        //fprintf (stdout, Eval (omega), ":", Eval (omegaHLA), "\n");
        if (pps <= 0.05) {
            pos_selection_overall += 1;
        }
        Export (lfExport, localLF_codon);
        filename = "../byPID/" + key + "_codon.lf";
        fprintf (filename,  CLEAR_FILE, lfExport);
        
    }
    
    time = 0;
 
    conversion_factor   = rateToBL("TN93",0);
    if (do_codon) {
        conversion_factorMG = rateToBL("MG94",dNdS_stensils["Syn"]);
    }
    
    //fprintf (stdout, "\n\n", conversion_factorMG, "\n\n");
    
    
    
    if (sample_count > 2) {
        ExecuteCommands ("`global_rate`:=`global_rate`__");
        ExecuteCommands ("Tree `treeName`_local = " + construct_tree_string (dates_to_sort[-1][2]));
        ExecuteCommands ("ReplicateConstraint (\"this1.?.?:=`global_rate`*this2.?.?__\",`treeName`_local,`treeName`_local)");
        bNames = Eval ("BranchName(`treeName`,-1)");
        
        
        ExecuteCommands ("global `global_rate`_slope=`global_rate`;`global_rate`_slope:>-100;");
        ExecuteCommands ("global `global_rate`_intercept=0;");
        
       df   = 0;
       time = 0;
       
        for (bID = Columns (bNames)-2; bID >= 0 ; bID += -1) {
            if ((bNames[bID]$"^[0-9]+$")[0] == 0) {
                continue;
            }
            ExecuteCommands (treeName + "_local." + bNames[bID] + ".t= " + treeName_local + "." + bNames[bID] + ".t");
            //fprintf (stdout, treeName + "_local." + bNames[bID] + ".t:=`global_rate`_slope*" + time/365 + "\n");
            df += 1;
            time += dates_to_sort[df][2]/2;
            /*ExecuteCommands (treeName + "_local." + bNames[bID] + ".t:=`global_rate`_slope*" + time/365 + "+`global_rate`_intercept;");
            time += dates_to_sort[df][2]/2;*/
        }
        clockReject [is_di] +=1;
        
        ExecuteCommands ("DataSetFilter `filterName`_local = CreateFilter (pol, 1,, filter_string)");
        ExecuteCommands ("LikelihoodFunction localLF2 = (`filterName`_local,`treeName`_local);");
        
        Optimize (res_local, localLF2);
        idx = 1;
        div = 0;
        distances = {sample_count, 2 + doHLA};
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
            
            //fprintf (stdout, "\n", bNames[bID], ":", bll, " : ", div, ":", time/365);
            idx += 1;
        }
        Export (lfExport, localLF2);
        filename = "../byPID/" + key + "_local.lf";
        fprintf (filename,  CLEAR_FILE, lfExport);
        //fprintf (stdout, "\n", key, "\n", distances, "\n", localLF2);
        // p = 1-CChi2(2*(res_local[1][0] - res[1][0]), df -1 );
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
    
    /*if (Max(distances[-1][1],0) > 0.01) {
        fprintf (stdout, distances);
    }*/
    
    /*bNames = Eval ("BranchName(`treeName`,-1)");
    fprintf (stdout, Eval ("Format(`treeName`, 1,1)"), "\n");
    for (bID = 0; bID < Columns (bNames); bID += 1) {
        fprintf (stdout, bNames[bID], " : ", Eval ("`treeName`." + bNames[bID] + ".t"), "\n");
    }*/
    
    //assert (sample_count < 3);
    
    GetDataInfo (siteDifferenceCount, computeDistanceFrom, 0, sample_count-1, AVERAGE_AMBIGUITIES);
    GetDataInfo (siteDifferenceCountR, computeDistanceFrom, 0, sample_count-1, RESOLVE_AMBIGUITIES);
    followup = +dates_to_sort[-1][2];
    
   
    //if (p <= 0.05 && !is_di) {
        if (doHLA) {
            fprintf (stdout, key, "\t", sample_count, "\t", Format(Eval (global_rateHLA + "*conversion_factor"),12,6), "\t", Format(Eval (global_rate + "*conversion_factor"),12,6), "\t", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCount)/followup*365,12,6), "\t", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCountR)/followup*365,12,6), "\t", 
            Format(p,8,4), "\t", Eval("`codon_rate`*conversion_factorMG"), "\t", followup, "\t", is_di, "\t", Format(p_hla,8,4),"\t", Format(p_hla2,8,4));
        
        } else {
            fprintf (stdout, key, "\t", sample_count, "\t", Format(Eval (global_rate + "*conversion_factor"),12,6), "\t", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCount)/followup*365,12,6), "\t", Format(ComputeDistanceFormulaFromDiffMx(siteDifferenceCountR)/followup*365,12,6), "\t", 
            Format(p,8,4), "\t", Eval("`codon_rate`*conversion_factorMG"), "\t", followup, "\t", is_di);
        }
        
        if (do_codon) {
            fprintf (stdout, "\t", Format(Eval (omega),5,2), "\t", Format(pps,4,2));
        }
        fprintf (stdout, "\n");
    //}
    
    file_name = "../byPID/" + key + ".fas";
    ExecuteCommands ("fprintf(file_name,CLEAR_FILE,`filterName`);");

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

