import pandas as pd

selected_attrs = [
    'AMT_CREDIT',
    'CNT_FAM_MEMBERS',
    'ENTRANCES_AVG',
    'BASEMENTAREA_AVG',
    'APARTMENTS_AVG',
    'OBS_30_CNT_SOCIAL_CIRCLE',
    'AMT_GOODS_PRICE',
    'FLOORSMAX_AVG',
    'FLOORSMIN_MEDI',
    'LIVINGAREA_AVG',
    'OBS_60_CNT_SOCIAL_CIRCLE',
    'CNT_CHILDREN',
    'LIVINGAPARTMENTS_AVG',
    'AMT_ANNUITY',
    'ELEVATORS_AVG',
    'COMMONAREA_AVG',
]

n = 60000
filepath = "data/train.csv"
output_filepath = f"data/loan.csv"

df = pd.read_csv(filepath, nrows=n)
df = df[selected_attrs]
df.fillna(df.mean(), inplace=True)
        
for col in df:
    print(df[col].describe())
    # normalize each attr to [-1, 1]
    df[col].clip(upper=df[col].quantile(q=0.95), inplace=True)
    df[col] -= (df[col].min() + df[col].max()) / 2
    df[col] /= df[col].max() + 1e-5
    print("after:", df[col].min(), df[col].max(), df[col].mean())
    print(df[col].describe())
    print("="*20)

df.to_csv(output_filepath, index=False, header=False)
