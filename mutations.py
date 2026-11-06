import pandas as pd

# 1. Load the tabular mutation data
df = pd.read_csv("data.txt")
print("Original Dataset Structure:\n", df.head())

# 2. Extract basic metrics
total_mutations = len(df)
print(f"\nTotal mutations tracked: {total_mutations}")

# 3. Filter data: Keep only 'High' impact mutations
high_impact_df = df[df['Impact'] == 'High']

# 4. Sort the filtered data by Chromosome Position
sorted_high_impact = high_impact_df.sort_values(by='Position')

print("\nFiltered & Sorted High Impact Mutations:\n", sorted_high_impact)

# 5. Export the cleaned sheet to a new file
sorted_high_impact.to_csv("cleaned_high_impact_mutations.csv", index=False)
