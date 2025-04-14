import pandas as pd
import yaml, sys, os, re

def create_config(input_csv_fhs):
    config = {}
    for input_csv in input_csv_fhs:
        dataset = {}
        #read in primer sheet
        df = pd.read_csv(input_csv)
        #filter to samples with indexes
        df2 = df[(pd.notna(df["fwd_index"])) & (pd.notna(df["rev_index"]))]
        for row in range(df2.shape[0]):            
            sec_str_primer = df2.iloc[row]["Forward_Primer_2ndRd_Sequence"]
		
            if "Index" in df2.iloc[row]["fwd_index"]:
                Index = "Index_primer"
            elif "S5" in df2.iloc[row]["fwd_index"]:
                Index = "Index_primer"    
            else:
                Index = "None"
                
            rev_primer = df2.iloc[row]["Reverse_Primer_2ndRd_Sequence"]
            fwd_index = df2.iloc[row]["fwd_index"]
            rev_index = df2.iloc[row]["rev_index"]
									
            dataset[df2.iloc[row]["Sample"]] = {
                "sec_str_primer": sec_str_primer,
                "index_type": Index,
                "fwd_index": fwd_index,
                "rev_index": rev_index,
                "rev_primer": rev_primer
                }

        config[os.path.basename(input_csv).split('.')[0]] = dataset
    sys.stdout.write(yaml.dump(config))

create_config(sys.argv[1:])
