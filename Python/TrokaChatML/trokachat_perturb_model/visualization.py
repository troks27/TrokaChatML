"""
Visualization Module

Functions for plotting and visualizing cell communication data and model results.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
import lightgbm as lgb
from matplotlib_venn import venn2


def plot_venn_diagram(control_df, exp_df, control_condition="NL", exp_condition="LS"):
    """
    Create a Venn diagram showing overlap between control and experimental conditions.
    
    Parameters:
    -----------
    control_df : pandas.DataFrame
        Control condition DataFrame
    exp_df : pandas.DataFrame
        Experimental condition DataFrame
    control_condition : str, default="NL"
        Label for control condition
    exp_condition : str, default="LS"
        Label for experimental condition
    """
    # Merge the DataFrames to find common and unique rows
    merged_df = exp_df.merge(
        control_df, 
        on=['Ligand', 'Receptor', 'Source', 'Target'], 
        suffixes=('_exp', '_control'), 
        how='outer',
        indicator=True
    )
    
    # Identify rows based on the adjusted_p_value conditions
    common_condition = (
        (merged_df['_merge'] == 'both') & 
        (
            ((merged_df['adjusted_p_value_exp'] < 0.05) & (merged_df['pct_filter_exp'] == 1)) | 
            ((merged_df['adjusted_p_value_control'] < 0.05) & (merged_df['pct_filter_control'] == 1))
        )
    )
    
    # Condition for rows that are only in the experimental DataFrame
    exp_only_condition = (
        (merged_df['_merge'] == 'left_only') & 
        (merged_df['adjusted_p_value_exp'] < 0.05) & 
        (merged_df['pct_filter_exp'] == 1)
    )

    # Condition for rows that are only in the control DataFrame
    control_only_condition = (
        (merged_df['_merge'] == 'right_only') & 
        (merged_df['adjusted_p_value_control'] < 0.05) & 
        (merged_df['pct_filter_control'] == 1)
    )
    
    # Count the number of rows that match each condition
    num_common = common_condition.sum()
    num_exp_only = exp_only_condition.sum()
    num_control_only = control_only_condition.sum()

    # Plot Venn diagram
    plt.figure(figsize=(8, 8))
    venn = venn2(subsets=(num_exp_only, num_control_only, num_common), 
                 set_labels=(exp_condition, control_condition))
    
    # Customize Venn diagram appearance
    venn.get_label_by_id('10').set_text(f'Only {exp_condition}\n{num_exp_only}')
    venn.get_label_by_id('01').set_text(f'Only {control_condition}\n{num_control_only}')
    venn.get_label_by_id('11').set_text(f'Both\n{num_common}')

    plt.title('Venn Diagram of Matching Rows')
    plt.show()


def plot_feature_importance(model, importance_type="gain"):
    """
    Plot feature importance for a LightGBM model.
    
    Parameters:
    -----------
    model : lightgbm.LGBMRegressor
        Trained LightGBM model
    importance_type : str, default="gain"
        Type of importance to plot ("gain" or "split")
    """
    lgb.plot_importance(model, importance_type=importance_type, 
                       figsize=(7, 6), 
                       title=f"LightGBM Feature Importance ({importance_type.capitalize()})")
    plt.show()


def plot_predicted_vs_actual(df_real, x_col, y_col, p_value_col='adjusted_p_value', threshold=0.05):
    """
    Plot predicted vs actual values with significant points highlighted.
    
    Parameters:
    -----------
    df_real : pandas.DataFrame
        DataFrame containing the data
    x_col : str
        Column name for x-axis (predicted values)
    y_col : str
        Column name for y-axis (actual values)
    p_value_col : str, default='adjusted_p_value'
        Column name for p-values
    threshold : float, default=0.05
        P-value threshold for significance
    """
    plt.figure(figsize=(5, 5))

    # Highlight significant points in red
    significant = df_real[p_value_col] < threshold
    plt.scatter(df_real.loc[significant, x_col], df_real.loc[significant, y_col], 
               color='red', alpha=0.6, label=f'Significant (p < {threshold})')

    # Plot non-significant points in blue
    plt.scatter(df_real.loc[~significant, x_col], df_real.loc[~significant, y_col], 
               color='blue', alpha=0.3, label='Non-significant')

    plt.plot([df_real[y_col].min(), df_real[y_col].max()],
             [df_real[y_col].min(), df_real[y_col].max()],
             'r--', linewidth=2)

    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Predicted vs Actual Values on Real Data")
    plt.legend()
    plt.show()


def plot_residuals_distribution(residuals, save_path=None):
    """
    Plot the distribution of residuals.
    
    Parameters:
    -----------
    residuals : array-like
        Residuals from model prediction
    save_path : str, optional
        Path to save the plot, if None doesn't save
    """
    # Set theme
    sns.set_theme(context="talk", style="white", font_scale=1.2)
    
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    sns.histplot(residuals, kde=True, color='gray', edgecolor='black', ax=ax)
    
    # Remove top and right spines
    sns.despine()

    # Increase font sizes
    ax.set_xlabel("Residual", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)

    # Clean layout
    plt.tight_layout()

    # Save as vector if desired
    if save_path:
        plt.savefig(save_path, format='pdf')
    
    plt.show()


def plot_predicted_vs_actual_values(y_pred, y_test, save_path=None):
    """
    Create a scatter plot of predicted vs actual values.
    
    Parameters:
    -----------
    y_pred : array-like
        Predicted values
    y_test : array-like
        Actual values
    save_path : str, optional
        Path to save the plot, if None doesn't save
    """
    sns.set_theme(context="talk", style="white", font_scale=1.2)
    
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    
    # Scatter plot
    ax.scatter(y_pred, y_test, alpha=0.5, color='black', s=10, edgecolors='none')
    
    # Perfect prediction line
    min_val = min(y_test.min(), y_pred.min())
    max_val = max(y_test.max(), y_pred.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=1)
    
    # Remove top and right spines
    sns.despine()

    # Labels
    ax.set_xlabel("Predicted", fontsize=12)
    ax.set_ylabel("Actual", fontsize=12)

    # Tight layout
    plt.tight_layout()

    # Save as vector if desired
    if save_path:
        plt.savefig(save_path, format='pdf')
    
    plt.show()


def plot_qq(residuals, save_path=None):
    """
    Create a Q-Q plot of the residuals to check normality.
    
    Parameters:
    -----------
    residuals : array-like
        Residuals from model prediction
    save_path : str, optional
        Path to save the plot, if None doesn't save
    """
    sns.set_theme(context="talk", style="white", font_scale=1.2)
    
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    
    # Create a Q-Q plot using scipy.stats.probplot
    stats.probplot(residuals, dist="norm", plot=ax)
    
    # Cleanup
    sns.despine()
    ax.set_title("")  # Remove if you want a cleaner look
    ax.set_xlabel("Theoretical Quantiles", fontsize=12)
    ax.set_ylabel("Ordered Residuals", fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, format="pdf")
    
    plt.show()


def plot_standardized_residuals_vs_predicted(y_pred, residuals, save_path=None):
    """
    Plot standardized residuals against predicted values.
    
    Parameters:
    -----------
    y_pred : array-like
        Predicted values
    residuals : array-like
        Residuals from model prediction
    save_path : str, optional
        Path to save the plot, if None doesn't save
    """
    sns.set_theme(context="talk", style="white", font_scale=1.2)
    
    # Compute standardized residuals
    std_residuals = (residuals - np.mean(residuals)) / np.std(residuals)
    
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    
    # Scatter plot
    ax.scatter(y_pred, std_residuals, alpha=0.5, edgecolors='none', color='black')
    
    # Add a horizontal line at 0 for reference
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1)
    
    sns.despine()
    ax.set_xlabel("Predicted", fontsize=12)
    ax.set_ylabel("Standardized Residuals", fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, format="pdf")
    
    plt.show()


def plot_elbow_residuals(residuals, rmse=None, save_path=None):
    """
    Create an elbow plot of sorted residuals.
    
    Parameters:
    -----------
    residuals : array-like
        Residuals from model prediction
    rmse : float, optional
        RMSE value to draw as horizontal line
    save_path : str, optional
        Path to save the plot, if None doesn't save
    """
    # Sort the data
    sorted_data = sorted(residuals)
    
    # Calculate the median
    median_value = np.median(sorted_data)
    
    # Create the elbow plot
    plt.figure(figsize=(10, 6))
    plt.plot(sorted_data, marker='o')
    plt.title('Elbow Plot of Residuals')
    plt.xlabel('Index')
    plt.ylabel('Residuals')
    
    # Add dashed lines at +/- RMSE if provided
    if rmse is not None:
        plt.axhline(y=rmse, color='r', linestyle='--', label=f'RMSE = {rmse}')
        plt.axhline(y=-rmse, color='r', linestyle='--')
    
    # Add a line for the median value
    plt.axhline(y=median_value, color='g', linestyle='-', label=f'Median = {median_value:.2f}')
    
    plt.legend()
    
    if save_path:
        plt.savefig(save_path, format="pdf")
    
    plt.show()
