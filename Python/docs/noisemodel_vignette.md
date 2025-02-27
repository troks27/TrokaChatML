# TrokaChat Noise Model Package - Jupyter Notebook Vignette

This vignette demonstrates how to use the TrokaChat Noise Model package for analyzing cellular communication data and identifying significant pathway interactions.

## Introduction

The TrokaChat Noise Model uses machine learning to distinguish true biological signals from random noise in cell-cell communication data. It builds a predictive model on null distribution data, then applies this model to real experimental data to identify statistically significant interactions.

## 1. Setup and Configuration

First, import the necessary modules:

```python
from trokachat_noise_model import data_processing as dp
from trokachat_noise_model import model_training as mt
from trokachat_noise_model import evaluation as ev
from trokachat_noise_model import plotting as pl
from trokachat_noise_model import config
from sklearn.model_selection import train_test_split
```

Check the current configuration:

```python
print("Current configuration:")
print("1. filename:", config.filename)
print("2. real_data_filename:", config.real_data_filename)
print("3. output_model_path:", config.output_model_path)
print("4. output_predictions_path:", config.output_predictions_path)
```

You can update the configuration variables as needed:

```python
# Update config variables with your data
config.filename = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/NULLDIST_NG_vs_NG_allpathways.h5"
config.real_data_filename = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/NG DEGs_allpathways.xlsx"
config.output_model_path = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/Noise Model/noisemodel_NG.pkl"
config.output_predictions_path = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/Noise Model/noise_model_predictions_NG.csv"

print("\nUpdated configuration:")
print("1. filename:", config.filename)
print("2. real_data_filename:", config.real_data_filename)
print("3. output_model_path:", config.output_model_path)
print("4. output_predictions_path:", config.output_predictions_path)
```

## 2. Loading and Preparing Data

List available groups in the HDF5 file and split them into training and testing sets:

```python
# List groups in the HDF5 file
groups = dp.list_h5_groups(config.filename)

# Split datasets into training and testing groups
train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)
```

## 3. Model Training

Train a LightGBM model using K-fold cross-validation with Bayesian hyperparameter optimization:

```python
# Perform K-Fold CV with Bayesian optimization
final_model, best_params = mt.perform_kfold_cv(
    train_groups, 
    config.filename, 
    n_trials=50,  # Number of hyperparameter combinations to try
    n_jobs=None,  # Use all available CPU cores
    batch_size=100  # Process data in batches to manage memory
)

print("Best hyperparameters:", best_params)
```

Save the trained model:

```python
# Save the trained model
mt.save_model(final_model, config.output_model_path)
```

## 4. Model Evaluation on Test Data

Preprocess the test data and evaluate the model performance:

```python
# Preprocess test data
batch_size = 100
X_test, y_test = dp.preprocess_test_data(test_groups, config.filename, batch_size, n_jobs=9)

# Get statistics about the data
stdev_logTF_communication_score = y_test.std()
print("Standard deviation of logTF communication score:", stdev_logTF_communication_score)

range_logTF_communication_score = y_test.max() - y_test.min()
print("Range of logTF communication score:", range_logTF_communication_score)

# Evaluate model performance on test data
rmse, r2, test_residuals = ev.evaluate_model_with_plots(final_model, X_test, y_test)
```

## 5. Analysis of Real Experimental Data

Load and preprocess real experimental data:

```python
import pandas as pd
import numpy as np

# Load real data
real_data_df = pd.read_excel(config.real_data_filename)
X_real, y_real = dp.preprocess_real_data(real_data_df)

# Evaluate model on real data
rmse, r2, real_residuals, real_data_df, y_real_pred = ev.evaluate_model_real_with_plots(
    final_model, X_real, y_real, real_data_df
)
```

Calculate p-values to identify significant interactions:

```python
# Calculate p-values and save results
ev.calculate_p_values(
    test_residuals, 
    real_residuals, 
    real_data_df, 
    X_real, 
    y_real, 
    y_real_pred, 
    config.output_predictions_path
)
```

## 6. Feature Importance Analysis

Visualize the importance of different features in the model:

```python
# Plot feature importance
pl.plot_feature_importance(final_model)
```

## 7. Running the Pipeline for Another Condition

You can run the same analysis pipeline for different conditions by updating the configuration:

```python
# Update configuration for a different condition (e.g., DIAB vs NG)
config.filename = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/NULLDIST_DIAB_vs_NG_allpathways.h5"
config.real_data_filename = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/HERE/DIAB DEGs_allpathways.xlsx"
config.output_model_path = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/Noise Model/noisemodel_DIAB.pkl"
config.output_predictions_path = "/Users/troks27/Desktop/DIAB_NG_TrokaChatML/TrokaChat/Noise Model/noise_model_predictions_DIAB.csv"

print("\nUpdated configuration for DIAB analysis:")
print("1. filename:", config.filename)
print("2. real_data_filename:", config.real_data_filename)
print("3. output_model_path:", config.output_model_path)
print("4. output_predictions_path:", config.output_predictions_path)
```

Then repeat the analysis steps for the new condition.

## 8. Complete Example: Using a Pre-Trained Model

If you've already trained a model, you can load it and proceed with evaluation:

```python
# Load a previously trained model
final_model = mt.load_model(config.output_model_path)

# List groups in the HDF5 file
groups = dp.list_h5_groups(config.filename)

# Split datasets into training and testing groups
train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)

# Preprocess and evaluate on the test data
batch_size = 100
X_test, y_test = dp.preprocess_test_data(test_groups, config.filename, batch_size, n_jobs=9)

# Evaluate on test data
rmse, r2, test_residuals = ev.evaluate_model_with_plots(final_model, X_test, y_test)

# Load and preprocess real data
real_data_df = pd.read_excel(config.real_data_filename)
X_real, y_real = dp.preprocess_real_data(real_data_df)

# Evaluate on real data
rmse, r2, real_residuals, real_data_df, y_real_pred = ev.evaluate_model_real_with_plots(
    final_model, X_real, y_real, real_data_df
)

# Calculate p-values and save results
ev.calculate_p_values(
    test_residuals, 
    real_residuals, 
    real_data_df, 
    X_real, 
    y_real, 
    y_real_pred, 
    config.output_predictions_path
)

# Plot feature importance
pl.plot_feature_importance(final_model)
```

## Conclusion

This vignette demonstrates how to use the TrokaChat Noise Model package to:
1. Train machine learning models on null distribution data
2. Apply the models to real experimental data
3. Identify statistically significant cell-cell communication interactions
4. Visualize and interpret the results

The significant interactions (p < 0.05) identified by this approach represent communication pathways that likely have biological relevance in your experimental conditions.