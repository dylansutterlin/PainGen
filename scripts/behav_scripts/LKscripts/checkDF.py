import pandas as pd

dir = '/home/dsutterlin/projects/genPain/'
df = pd.read_csv(dir + 'results/' + 'SCEBLmri_Finaldata_TxT_N36.csv')
print(df.shape, df.head)
print(df.describe())