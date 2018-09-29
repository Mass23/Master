import numpy as np
import pandas as pd

data = pd.read_csv("test.csv")
data = data.as_matrix(["gene_id","mkr","polymorphism","substitution","gc_content","uncertainty","length","identity"])
print(data)
