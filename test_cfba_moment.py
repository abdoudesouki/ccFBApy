#6/2/2024
test_cfba_moment():

import os
os.chdir(os.path.expanduser("D:\\Dropbox\\py\\Bioinfo_py\\ccFBApy"))
import mod2irrev

import cobra
model = cobra.io.read_sbml_model('E:\\Cobra\\iAF1260.xml')

mod2=mod2irrev(model)
import pandas as pd
mw=pd.read_csv("D:\\Dropbox\\Sybil\\CCFBA\\MOMENT\\mw.csv")
kcat=pd.read_csv("D:\\Dropbox\\Sybil\\CCFBA\\MOMENT\\kcat.csv")

sol=cfba_moment(model,mod2,kcat,MW=mw,verbose=3,RHS=0.27,medval=3600*22.6)
