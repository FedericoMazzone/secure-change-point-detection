import pandas as pd

n = 100000
filepath = "data/train.csv"
output_filepath = f"data/taxi.csv"

time_attrs = ['pickup_datetime', 'dropoff_datetime']
num_attrs = ['passenger_count', 'pickup_longitude', 'pickup_latitude',  'dropoff_longitude',
                'dropoff_latitude', 'trip_duration']

df = pd.read_csv(filepath, nrows=n+500)
df = df[df['trip_duration'] < 5000]
df = df.head(n=n)
for col in df:
    if col in time_attrs:
        df[col] = pd.to_datetime(df[col])
        print("before:", df[col].min(), df[col].max(), df[col].mean())
        df[col] = df[col].astype(int) / 1e9
        # only consider the time in a week
        df[col] = df[col].mod(60 * 60 * 24 * 7)
    if col in num_attrs or col in time_attrs:
        # normalize each attr to [-1, 1]
        df[col] -= (df[col].min() + df[col].max()) / 2
        df[col] /= df[col].max()
        print("after:", df[col].min(), df[col].max(), df[col].mean())

# store only columns in time_attrs and num_attrs, and without header
df = df[time_attrs + num_attrs]
df.to_csv(output_filepath, index=False, header=False)