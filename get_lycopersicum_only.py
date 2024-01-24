import pandas as pd
ndf = pd.read_csv("/home/adejoro/Desktop/from-server/filtered/compare.csv")
print(ndf)

ndf_e = ndf[(ndf['empty'] == 1)]
print(ndf_e)