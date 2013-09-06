array_dimension = 1000;
min_score       = 10;

LoadFunctionLibrary         ("ReadDelimitedFiles");
needMapFor = {};


lfunction readDRAM (filePath, array_dimension, offset, mappingArray) {
    fscanf (filePath, REWIND, "Lines", allLines);
    dnames = splitStringByTab (allLines[0]);
    coordToName     = {};
    for (d = 2; d < Abs (dnames); d+=1) {
        (^mappingArray) [dnames[d]] = {1, array_dimension*3} ["0"];
        coordToName [d] = dnames[d];
    }
    
    
    for (l = 1; l < Columns (allLines); l += 1) {
        locs = splitStringByTab (allLines[l]);
        residue = ((-1) + locs[0] + offset) * 3;
        for (d = 2; d < Abs (locs); d+=1) {
            if (0 + locs[d] >= min_score) {
                ((^mappingArray)[coordToName[d]])[residue] += 1;
                ((^mappingArray)[coordToName[d]])[residue+1] += 1;
                ((^mappingArray)[coordToName[d]])[residue+2] += 1;
            }
        }
    }
    
    return 0;
}

dram_array = {};

fscanf ("files.lst", "Lines", allFiles);

for (k = 0; k < Columns (allFiles); k += 1) {
    bits = splitOnRegExp(allFiles[k], ",");
    assert (Abs (bits) == 2, "Incorrect DRAM file list line specification: '" + allFiles[k] + "'");
    readDRAM (bits[0], array_dimension, 0 + bits[1], &dram_array);
}

SetDialogPrompt ("Load a .csv file with DRAM information");
dram_data = ReadCSVTableText ("",0);

feature_filter_by_pid = {};

for (k = 0; k < Abs (dram_data); k+=1) {
    id = (dram_data[k])[0];
    patient_selector = {1,array_dimension*3}["0"];
    for (d = 1; d < Abs (dram_data[k]); d += 1) {
        drug_name = ((dram_data[k])[d] && 1) ^ {{"\\ ",""}};
        if (drug_name == "RTV") {
            continue;
        }
        dd = dram_array[drug_name];
        assert (Abs (dd), drug_name + " is an invalid ARV med name");
        patient_selector += dd;
    }
    feature_filter_by_pid[id] = patient_selector;
}


