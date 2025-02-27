# plotting.py
import matplotlib.pyplot as plt
import seaborn as sns
import lightgbm as lgb

def plot_feature_importance(model):
    lgb.plot_importance(model, importance_type="gain", figsize=(7, 6), title="LightGBM Feature Importance (Gain)")
    plt.show()
    
    lgb.plot_importance(model, importance_type="split", figsize=(7, 6), title="LightGBM Feature Importance (Split)")
    plt.show()

def plot_predicted_vs_actual(df_real, x_col, y_col, p_value_col='adjusted_p_value'):
    plt.figure(figsize=(10, 6))

    # Highlight significant points in red
    significant = df_real[p_value_col] < 0.05
    plt.scatter(df_real.loc[significant, x_col], df_real.loc[significant, y_col], color='red', alpha=0.6, label='Significant (p < 0.05)')

    # Plot non-significant points in blue
    plt.scatter(df_real.loc[~significant, x_col], df_real.loc[~significant, y_col], color='blue', alpha=0.3, label='Non-significant')

    plt.plot([df_real[y_col].min(), df_real[y_col].max()],
             [df_real[y_col].min(), df_real[y_col].max()],
             'r--', linewidth=2)

    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Predicted vs Actual Values on Real Data")
    plt.legend()
    plt.show()

def plot_residuals_distribution(residuals):
    plt.figure(figsize=(10, 6))
    sns.histplot(residuals, kde=True)
    plt.title("Residuals Distribution")
    plt.xlabel("Residual")
    plt.ylabel("Frequency")
    plt.show()

def plot_predicted_vs_actual_values(y_pred, y_test):
    plt.figure(figsize=(10, 6))
    plt.scatter(y_pred, y_test, alpha=0.3)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', linewidth=2)
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Predicted vs Actual Values")
    plt.show()
