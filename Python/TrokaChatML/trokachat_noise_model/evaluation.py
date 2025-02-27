# evaluation.py
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from . import plotting as pl

def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    r2 = r2_score(y_test, y_pred)
    print(f"RMSE on test data: {rmse:.4f}")
    print(f"R² on test data: {r2:.4f}")
    residuals = y_test - y_pred
    return rmse, r2, residuals, y_pred, y_test 

def evaluate_model_with_plots(model, X_test, y_test):
    rmse, r2, residuals, y_pred, y_test = evaluate_model(model, X_test, y_test)
    
    pl.plot_residuals_distribution(residuals)
    pl.plot_predicted_vs_actual_values(y_pred, y_test)
    
    return rmse, r2, residuals

def evaluate_model_real(model, X_real, y_real, real_data_df):
    # Predict on the real data
    y_real_pred = model.predict(X_real)
    rmse = np.sqrt(mean_squared_error(y_real, y_real_pred))
    r2 = r2_score(y_real, y_real_pred)
    print(f"RMSE on real data: {rmse:.4f}")
    print(f"R² on real data: {r2:.4f}")
    residuals = y_real - y_real_pred
    return rmse, r2, residuals, real_data_df, y_real, y_real_pred


def evaluate_model_real_with_plots(model, X_real, y_real, real_data_df):
    rmse, r2, residuals, real_data_df, y_real, y_real_pred, = evaluate_model_real(model, X_real, y_real, real_data_df)

    pl.plot_residuals_distribution(residuals)
    pl.plot_predicted_vs_actual_values(y_real_pred, y_real)
    
    return rmse, r2, residuals, real_data_df, y_real_pred

def calculate_p_values(test_residuals, real_residuals, real_data_df, X_real, y_real, y_real_pred, output_filepath):
    mean_null_residuals = np.mean(test_residuals)
    std_null_residuals = np.std(test_residuals)
    z_scores_real = (real_residuals - mean_null_residuals) / std_null_residuals
    p_values_real = norm.sf(z_scores_real)
    _, p_values_adjusted, _, _ = multipletests(p_values_real, method='fdr_bh')

    df_real = pd.DataFrame(X_real)
    df_real['Communication Score'] = real_data_df['Communication Score']
    df_real['Uniqueness Score'] = real_data_df['Uniqueness Score']
    df_real['logTF Communication Score'] = y_real
    df_real['Predicted Communication Score'] = y_real_pred
    df_real['Residual'] = real_residuals
    df_real['z_score'] = z_scores_real
    df_real['p_value'] = p_values_real
    df_real['adjusted_p_value'] = p_values_adjusted

    df_real.to_csv(output_filepath, index=False)
    print(f"Predictions with p-values saved to {output_filepath}")

    pl.plot_predicted_vs_actual(df_real, 'Predicted Communication Score', 'logTF Communication Score')
