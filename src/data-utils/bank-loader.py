import pandas as pd
from sklearn.preprocessing import StandardScaler

# Step 1: Load the dataset
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00222/bank-additional.zip"
file_name = "bank-additional/bank-additional-full.csv"
data_zip = "bank_additional.zip"
data_csv = "data/bank.csv"

# Download the dataset
import requests, zipfile, io
response = requests.get(url, verify=False)
with zipfile.ZipFile(io.BytesIO(response.content)) as z:
    z.extract(file_name)

# Read the dataset
data = pd.read_csv(file_name, sep=';')

# Step 2: Select significant columns for clustering
columns_to_keep = ['age', 'duration', 'campaign', 'previous', 'emp.var.rate', 'cons.conf.idx', 'nr.employed']
data = data[columns_to_keep]

# Remove records with outliers below the 5th percentile and above the 95th percentile
data = data[(data >= data.quantile(0.05)) & (data <= data.quantile(0.95))].dropna()
# print number of records after removing outliers
print(f"Number of records after removing outliers: {len(data)}")

# Step 3: Normalize the data
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data)

# Step 4: Save the processed data to a CSV
final_data = pd.DataFrame(data_scaled, columns=columns_to_keep)
final_data.to_csv(data_csv, index=False, header=False)

print(f"Preprocessed dataset saved to {data_csv}")
