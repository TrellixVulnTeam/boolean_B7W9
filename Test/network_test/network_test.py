import pandas as pd
import sys
import os
from io import StringIO
from pprint import pprint
import re

sys.path.insert(0, "/Users/Oliver/Code/")
import bone


def get_file(file: str):
    file_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(file_dir, file)
    return file_path


def bone_bv():
    # create the survival and expression file with the GEO class
    gse40060 = bone.GEO(accessionID="GSE40060", download_soft=False)
    gse40060_survival = gse40060.survival()  # if unspecified, first GPL is used
    gse40060_expr = gse40060.expr(get_genes=True)
    gse40060_expr = gse40060_expr.fillna(0)  # expr file cannot contain any NA values

    # initialize bone with expression and survival file
    my_bone = bone.BoNE(expr=gse40060_expr, survival=gse40060_survival)

    # define variables for bone
    survival_name = "c source_name_ch1"
    # first group used (in this case 'endogenous') is always control
    gse40060_groups = {"E": "endogenous", "O": "overexpressed"}
    alz_weights = {
        -3: ["SVOP", "CACNG3", "PCYOX1L", "BEX1", "TUBB3", "NRN1", "GAP43", "RGS4"],
        1: ["BGN", "EHD2", "FCGRT", "NT5DC2", "ITGB5", "PDGFRB", "GPR4", "LAMB2"],
    }
    my_bone.init(
        survival_col=survival_name,
        gene_weights=alz_weights,
        groups=gse40060_groups,
    )
    bv = my_bone.bv()
    return bv


bv = bone_bv()
# bv = pd.read_csv(get_file("bv.csv"), index_col=(0, 1))

sm = bone.Stepminer(bv)
output = sm.network
print(output)
