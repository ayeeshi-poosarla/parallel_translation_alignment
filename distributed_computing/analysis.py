import csv
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from io import StringIO

# Function to load data from a string
def load_data_from_string(data_string):
    data = {}
    reader = csv.reader(StringIO(data_string))
    next(reader)  # Skip header
    for row in reader:
        data[row[0]] = row[1]  # {Reference_Gene_Name: Protein_Gene_Name}
    return data

# Ground truth and predictions as strings
ground_truth_data = """Reference_Gene_Name,Protein_Gene_Name
Gene1,Gene1_Protein
Gene2,Gene2_Protein
Gene3,Gene3_Protein
Gene4,Gene4_Protein
Gene5,Gene5_Protein
Gene6,Gene6_Protein
Gene7,Gene7_Protein
Gene8,Gene8_Protein
Gene9,Gene9_Protein
Gene10,Gene10_Protein
Gene11,Gene11_Protein
Gene12,Gene12_Protein
Gene13,Gene13_Protein
Gene14,Gene14_Protein
"""

predictions_data = """Reference_Gene_Name,Protein_Gene_Name
Gene1,Gene1_Protein
Gene2,Gene2_Protein
Gene3,Gene3_Protein
Gene4,Gene4_Protein
Gene5,Gene4_Protein
Gene6,Gene6_Protein
Gene7,Gene7_Protein
Gene8,Gene4_Protein
Gene9,Gene9_Protein
Gene10,Gene10_Protein
Gene11,Gene11_Protein
Gene12,Gene12_Protein
Gene13,Gene13_Protein
Gene14,Gene4_Protein
"""

# Load the data
ground_truth = load_data_from_string(ground_truth_data)
predictions = load_data_from_string(predictions_data)

# Align ground truth and predictions
y_true = []
y_pred = []

for gene in ground_truth:
    y_true.append(ground_truth[gene])
    y_pred.append(predictions.get(gene, "NO_MATCH"))  # Default to "NO_MATCH" if no prediction

# Calculate metrics
accuracy = accuracy_score(y_true, y_pred)
precision = precision_score(y_true, y_pred, average='micro', zero_division=0)  # 'micro' accounts for all labels
recall = recall_score(y_true, y_pred, average='micro', zero_division=0)
f1 = f1_score(y_true, y_pred, average='micro', zero_division=0)

# Print metrics
print(f"Accuracy: {accuracy}")
print(f"Precision: {precision}")
print(f"Recall: {recall}")
print(f"F1-Score: {f1}")
