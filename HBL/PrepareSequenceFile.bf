LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("HXB2Mapper");

_hxb_alignOptions_codon["SEQ_ALIGN_NO_TP"] = 1;
alignOptions["SEQ_ALIGN_NO_TP"] = 1;

ChoiceList (whichGene,"Which gene is being analyzed",1,SKIP_NONE,
                "pol", "polymerase",
                "gag", "gag/matrix",
                "nef", "nef"
            );
            
assert (whichGene >= 0);            
gene = SELECTION_STRINGS; 
fprintf (stdout, gene, "\n");
ref_seqs = {"pol": 15,
            "gag" :4,
            "nef":3
            };
            
SetDialogPrompt ("Load the raw sequence file");
DataSet       pol = ReadDataFile (PROMPT_FOR_FILE);
_inPath = LAST_FILE_PATH;
DataSetFilter raw_filter = CreateFilter (pol, 1);

realign = 0;

for (seq_id = 0; seq_id < raw_filter.species; seq_id += 1) {
    GetString   (seq_name, raw_filter, seq_id);  
    GetDataInfo (seq_value, raw_filter, seq_id);
      
    assert ((seq_name $ "^[^\\|]+\\|[0-9]{4}\\-[0-9]{2}\\-[0-9]{2}")[0]==0,
             "Sequence name `seq_name` is not of the form PatientID|YYYY-MM-DD");
          
    if (realign == 0) {         
        SetParameter (STATUS_BAR_STATUS_STRING, "[PrepareSequenceFile: Checking the alignment with HXB2_`gene` for `seq_name`]                                            ", 0);
        mp = mapSequenceToHXB2Aux (seq_value^{{"\\-",""}}, RefSeqs[13], 2); 
        if (mp [0] != 0) {
            fprintf (stdout, "\n[PrepareSequenceFile: Sequence `seq_name` is not aligned with the beginning of HXB2 `gene`]");   
            realign = 1;
        } else {
            for (k = 1; k < Rows (mp); k+=1) {
                if (mp [k] - mp[k-1] != 1) {
                    fprintf (stdout, "\n[PrepareSequenceFile: Sequence `seq_name` has indels relative to HXB2 `gene`]");
                    realign = 1;                
                    break;
                }
            }
        }   
    }
}

if (realign) {
    SetDialogPrompt ("Save sequences realigned to HBX2 `gene` to");
    fprintf (PROMPT_FOR_FILE,CLEAR_FILE);
    _outputPath = LAST_FILE_PATH;
    _path = _outputPath;
    checkReferenceSequenceForStopCodons = 1;
    doCodonAlignment        = 1;
    refSeq = ref_seqs[gene];
    LoadFunctionLibrary ("SeqAlignmentNucShared", {"0": "HXB2_`gene`", "1": _inPath, "2" : "No", "3": _outputPath});
    DataSet       ds = ReadDataFile (_outputPath);
    DataSetFilter filtered_data = CreateFilter (ds, 1, "", speciesIndex >= 1);
    DataSetFilter reference_seq = CreateFilter (ds, 1, "", speciesIndex == 0);
    
    GetDataInfo (refSeq, reference_seq, 0);
    gaps = (refSeq || "\\-" );
    if (gaps[0] >= 0) {
        fprintf (stdout, "\n[PrepareSequenceFile: IMPORTANT! At least one of the sequences has an INSERTION relative to HXB2. TO USE HLA and DRAM annotated data with these sequences, please use the xx_clipped file]\n");
        filtering_vector = {Abs(refSeq), 1};
        for (i = 0; i < Rows (gaps); i+=2) {
            filtering_vector[gaps[i]] = 1;
        }   
    } else {
        filtering_vector = 0;
    }
    
    
    DATA_FILE_PRINT_FORMAT = 9;
    Export (fs, filtered_data);
    DataSet       no_ref = ReadFromString (fs);
    DataSetFilter allGaps = CreateFilter (no_ref, 1, "/^\-+$/");
    if (Abs (filtering_vector)) {
        DataSetFilter clipped = CreateFilter (no_ref, 1, filtering_vector[siteIndex] == 0);
    }
    
    if (allGaps.sites != filtered_data.sites) {
        k = filtered_data.sites-1; 
        i = Columns (allGaps.site_map)-1;
            
        while (allGaps.site_map[i] == k) {
            i = i-1;
            k = k-1;
        }
        DataSetFilter allGaps = CreateFilter (no_ref, 1, siteIndex <= k);
        Export (fs, allGaps);
        
        if (Abs (filtering_vector)) {
            DataSetFilter clipped = CreateFilter (allGaps, 1, filtering_vector[siteIndex] == 0 && siteIndex <= k);
            Export (fsc, clipped);
        }    
        
    } else {
        if (Abs (filtering_vector)) {
            DataSetFilter clipped = CreateFilter (no_ref, 1, filtering_vector[siteIndex] == 0);
            Export (fsc, clipped);
        }    
    }
    
    fprintf (_outputPath, CLEAR_FILE, fs);
    if (Abs (filtering_vector)) {
        _outputPath += "_clipped";
        fprintf (_outputPath, CLEAR_FILE, fsc);
    }
    
    
} else {
    _path = _inPath;
}

fprintf (stdout, "\n[PrepareSequenceFile: Use `_path` with RateAnalysis.bf]\n");

