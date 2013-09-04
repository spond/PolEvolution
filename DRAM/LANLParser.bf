array_dimension = 550;


LoadFunctionLibrary         ("ReadDelimitedFiles");

SetDialogPrompt ("Load a .csv file with CTL haplotype information");
hla_data = ReadCSVTableText ("",1);

needMapFor = {};

superTypeLookup = {"A2 supertype" : {{"A*0201", "A*0202", "A*0203", "A*0206", "A*6802"}},
                   "A3 supertype" : {{"A*0301", "A*1101", "A*3101", "A*3301", "A*6801"}},
                   "A28 supertype" : {{"A*28"}},
                   };


lfunction simpleCTL (in) {
    if (in [0] == "Cw") {
        in [0] = "C";
    }
    if (Abs (in[1]) > 2) {
        in[1] = (in[1])[0][1];
    } else {
        if (Abs (in[1]) == 1) {
            in[1] = "0" + in[1];
        }
    }
    return in[0] + "*" + in[1];
}

lfunction standardizeCTL (in) {
    result = {};
    entry = extractSubexpressions (in, "([A-Za-z]+)\\*?([0-9]+)$", 0, 0);
    if (Abs (entry) == 0) {
        lookup = (^"superTypeLookup")[in];
        if (Abs(lookup)) {
            //fprintf (stdout, "\n\n\n", in, "\n", lookup, "\n\n\n");
            for (k = 0; k < Columns (lookup); k+=1) {
                result + (simpleCTL(splitOnRegExp(lookup[k],"\\*")));
            }
        } else {
            (^"needMapFor") [in] = 1;
        }
    } else {
        result [0] = simpleCTL (entry);
    }
    return result;
}

lfunction readSVG (filePath, array_dimension, mappingArray) {
    fscanf (filePath, REWIND, "Lines", allLines);
    // ctl_search?results=Search;protein=RT;inrange=18-26;hla=B*0801
    for (l = 0; l < Columns (allLines); l += 1) {
        entry = extractSubexpressions (allLines[l], "ctl_search\\?results=Search;protein=([A-Za-z]+);inrange=([0-9]+)\\-([0-9]+);hla=([^\\\"]+)", 0, "");
        
        if (Abs (entry)) {
            if (entry[0] == "RT") {
                residue_shift = 98;
            } else {
                residue_shift = (-1);
            }
            from = 3*(residue_shift + (0+entry[1]));
            to = 3*(residue_shift + (0+entry[2]))-1;
            if (to - from >= 20) {
                handled = standardizeCTL (entry[3]);
                
                for (k = 0; k < Abs (handled); k+=1) {
                    ctl = handled[k];
                    if (Abs((^mappingArray)[ctl]) == 0) {
                        (^mappingArray)[ctl] = {1,array_dimension*3}["0"];
                    }
                    ((^mappingArray)[ctl]) += {1,array_dimension*3} ["_MATRIX_ELEMENT_COLUMN_>="+from+"&&_MATRIX_ELEMENT_COLUMN_<="+to]; 
               }
                
            }
        }
    }
    
    return 0;
}

ctl_array = {};


fscanf ("files.lst", "Lines", allFiles);

for (k = 0; k < Columns (allFiles); k += 1) {
    readSVG (allFiles[k], array_dimension, &ctl_array);
}


hla_filter_by_pid = {};
hla_data = hla_data [1];

letter_by_index = {"4":"A", "5": "A", "6": "B", "7": "B", "8": "C", "9": "C"};

for (k = 0; k < Abs (hla_data); k+=1) {
    id = (hla_data[k])[12];
    patient_selector = {1,array_dimension*3}["0"];
    for (c = 4; c < 10; c+=1) {
        possibles = splitOnRegExp ((hla_data[k])[c], "/");
        for (p = 0; p < Abs (possibles); p += 1) {
            key = letter_by_index[c] + possibles[p];
            if (Abs (ctl_array[key]) > 0) {
                patient_selector += ctl_array[key];
            }
        }
    }
    hla_filter_by_pid[id] = patient_selector;
}

fprintf (stdout, hla_filter_by_pid, "\n");

