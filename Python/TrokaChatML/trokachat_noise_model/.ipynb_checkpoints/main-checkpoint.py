# main.py
from trokachat_noise_model import data_processing as dp
from trokachat_noise_model import model_training as mt
from trokachat_noise_model import evaluation as ev
from trokachat_noise_model import plotting as pl
from trokachat_noise_model import config
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np


# Import model if previously computed
final_model = mt.load_model(config.output_model_path)

# List groups in the HDF5 file
groups = dp.list_h5_groups(config.filename)

# Split datasets into training and testing groups
train_groups, test_groups = train_test_split(groups, test_size=0.2, random_state=42)

# Perform K-Fold CV with Bayesian optimization
final_model, best_params = mt.perform_kfold_cv(train_groups, config.filename, n_trials=50, n_jobs=None, batch_size=100)
print("Best hyperparameters:", best_params)

# Save model
mt.save_model(final_model, config.output_model_path)

# Preprocess and evaluate on the test data
batch_size = 100
X_test, y_test = dp.preprocess_test_data(test_groups, config.filename, batch_size, n_jobs=9)

# Evaluate performance and plot residuals on test data
rmse, r2, test_residuals = ev.evaluate_model_with_plots(final_model, X_test, y_test)

# Load and preprocess real data
real_data_df = pd.read_excel(config.real_data_filename)
X_real, y_real = dp.preprocess_real_data(real_data_df)

# Evaluate performance and plot residuals on real data
rmse, r2, real_residuals, real_data_df, y_real_pred = ev.evaluate_model_real_with_plots(final_model, X_real, y_real, real_data_df)

# Calculate adj. p-values, plot, and save results from real data
ev.calculate_p_values(test_residuals, real_residuals, real_data_df, X_real, y_real, y_real_pred, config.output_predictions_path)

# Plot model features
pl.plot_feature_importance(final_model)