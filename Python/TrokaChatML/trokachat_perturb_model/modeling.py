"""
Modeling Module

Functions for training, evaluating, and using machine learning models
for cell communication perturbation prediction.
"""

import pandas as pd
import numpy as np
import lightgbm as lgb
import pickle
import time
import warnings
import multiprocessing
from sklearn.metrics import root_mean_squared_error, r2_score, make_scorer
from sklearn.model_selection import KFold
from skopt import BayesSearchCV
from skopt.space import Real, Integer


def create_lgb_dataset(df):
    """
    Create a LightGBM dataset from a DataFrame.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing features and target
        
    Returns:
    --------
    tuple
        (lgb_data, feature_columns) - LightGBM dataset and feature column names
    """
    excluded_columns = [col for col in df.columns if col not in [
        'Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label', 
        'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor'
    ]]
    
    feature_columns = [col for col in df.columns if col not in excluded_columns]
    cat_indices = [i for i, col in enumerate(feature_columns) if df[col].dtype.name == 'category']
    
    X = df[feature_columns]
    for col in feature_columns:
        X[col] = X[col].astype('category')
    
    y = df['logTF Communication Score']
    lgb_data = lgb.Dataset(X, label=y, categorical_feature=feature_columns, free_raw_data=False)
    
    return lgb_data, feature_columns


def perform_kfold_cv(df, n_trials=50, n_jobs=None, random_state=42):
    """
    Perform Bayesian optimization with k-fold cross-validation to find the best hyperparameters.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing features and target
    n_trials : int, default=50
        Number of parameter settings that are sampled
    n_jobs : int, optional
        Number of jobs to run in parallel, if None uses (CPU cores - 1)
    random_state : int, default=42
        Random state for reproducibility
        
    Returns:
    --------
    tuple
        (final_model, best_params) - The best model and its parameters
    """
    warnings.filterwarnings("ignore")
    start_time = time.time()

    if n_jobs is None:
        num_cores = multiprocessing.cpu_count()
        n_jobs = num_cores - 1
    print(f"Using {n_jobs} CPU cores for parallel processing")

    lgb_train_data, feature_names = create_lgb_dataset(df)

    param_space = {
        'learning_rate': Real(0.01, 0.6, 'uniform'),
        'num_leaves': Integer(8, 1500),
        'colsample_bytree': Real(0.5, 1.0, 'uniform'),  # feature_fraction
        'subsample': Real(0.1, 1.0, 'uniform'),         # bagging_fraction
        'max_depth': Integer(3, 16),
        'min_child_samples': Integer(50, 500),          # min_data_in_leaf
        'min_split_gain': Real(0, 15),                  # min_gain_to_split
        'reg_lambda': Real(0.01, 100),                  # lambda_l2
    }

    bayes_search = BayesSearchCV(
        estimator=lgb.LGBMRegressor(
            objective='regression', 
            metric='rmse', 
            boosting_type='gbdt', 
            n_jobs=n_jobs, 
            verbose=-1,
            random_state=random_state
        ),
        search_spaces=param_space,
        n_iter=n_trials,
        cv=KFold(n_splits=5, shuffle=True, random_state=random_state),
        scoring=make_scorer(root_mean_squared_error, greater_is_better=False),
        n_jobs=n_jobs,
        return_train_score=True,
        random_state=random_state,
        refit=True
    )
    
    print("Starting Bayesian optimization with cross-validation...")
    bayes_search.fit(lgb_train_data.data, lgb_train_data.label, categorical_feature=feature_names)

    best_params = bayes_search.best_params_
    final_model = bayes_search.best_estimator_

    elapsed_time = time.time() - start_time
    print(f"Total Time elapsed: {elapsed_time:.2f} seconds")

    return final_model, best_params


def save_model(model, filepath):
    """
    Save a model to disk.
    
    Parameters:
    -----------
    model : object
        Model to save
    filepath : str
        Path where to save the model
    """
    with open(filepath, 'wb') as model_file:
        pickle.dump(model, model_file)
    print(f"Model saved as {filepath}")


def load_model(filepath):
    """
    Load a model from disk.
    
    Parameters:
    -----------
    filepath : str
        Path where the model is saved
        
    Returns:
    --------
    object
        Loaded model
    """
    with open(filepath, 'rb') as model_file:
        return pickle.load(model_file)


def preprocess_test_dataset(df):
    """
    Preprocess a dataset for testing with LightGBM.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame to preprocess
        
    Returns:
    --------
    tuple
        (lgb_data, feature_columns) - LightGBM dataset and feature column names
    """
    # Define columns to exclude from features
    excluded_columns = [col for col in df.columns if col not in [
        'Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label', 
        'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor'
    ]]

    # Select feature columns and identify categorical indices
    feature_columns = [col for col in df.columns if col not in excluded_columns]
    categorical_columns = [
        'Ligand', 'Receptor', 'Source', 'Target', 'Pathway Label', 
        'source_target', 'ligand_receptor', 'source_ligand', 'target_receptor'
    ]
    
    # Ensure categorical columns are of type 'category'
    for col in categorical_columns:
        df[col] = df[col].astype('category')
    
    # Prepare feature matrix and target vector
    X = df[feature_columns]
    y = df['logTF Communication Score']
    
    # Create LightGBM Dataset
    lgb_data = lgb.Dataset(X, label=y, categorical_feature=categorical_columns, free_raw_data=False)
    
    return lgb_data, feature_columns


def evaluate_model(model, X_test, y_test):
    """
    Evaluate a model's performance on test data.
    
    Parameters:
    -----------
    model : object
        Trained model to evaluate
    X_test : pandas.DataFrame
        Test features
    y_test : pandas.Series
        Test target values
        
    Returns:
    --------
    tuple
        (rmse, r2, residuals, y_pred, y_test) - Evaluation metrics and values
    """
    y_pred = model.predict(X_test)
    rmse = root_mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    print(f"RMSE on test data: {rmse:.4f}")
    print(f"R² on test data: {r2:.4f}")
    residuals = y_test - y_pred
    return rmse, r2, residuals, y_pred, y_test


def evaluate_model_with_plots(model, X_test, y_test):
    """
    Evaluate a model and create diagnostic plots.
    
    Parameters:
    -----------
    model : object
        Trained model to evaluate
    X_test : pandas.DataFrame
        Test features
    y_test : pandas.Series
        Test target values
        
    Returns:
    --------
    tuple
        (rmse, r2, residuals) - Evaluation metrics and residuals
    """
    from .visualization import (
        plot_residuals_distribution, 
        plot_predicted_vs_actual_values, 
        plot_qq, 
        plot_standardized_residuals_vs_predicted
    )
    
    rmse, r2, residuals, y_pred, y_test = evaluate_model(model, X_test, y_test)
    
    # Create diagnostic plots
    plot_residuals_distribution(residuals)
    plot_predicted_vs_actual_values(y_pred, y_test)
    plot_qq(residuals)
    plot_standardized_residuals_vs_predicted(y_pred, residuals)
    
    return rmse, r2, residuals


def train_final_model(df, best_params, random_state=42):
    """
    Train a final model with the best hyperparameters.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with training data
    best_params : dict
        Dictionary of best hyperparameters
    random_state : int, default=42
        Random state for reproducibility
        
    Returns:
    --------
    object
        Trained final model
    """
    lgb_data, feature_columns = preprocess_test_dataset(df)
    
    # Create the final model with the best hyperparameters
    final_model = lgb.LGBMRegressor(
        objective='regression',
        metric='rmse',
        boosting_type='gbdt',
        n_jobs=-1,  # Use all available cores
        verbose=-1,
        random_state=random_state,
        **best_params
    )
    
    # Train the model on the entire dataset
    final_model.fit(lgb_data.data, lgb_data.label, categorical_feature=feature_columns)
    
    return final_model


def predict_perturbations(model, df):
    """
    Predict perturbations using the trained model.
    
    Parameters:
    -----------
    model : object
        Trained model for predictions
    df : pandas.DataFrame
        DataFrame with features for prediction
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with predictions and residuals
    """
    # Preprocess the data
    lgb_data, feature_columns = preprocess_test_dataset(df)
    
    # Make predictions
    df['predicted_score'] = model.predict(lgb_data.data)
    
    # Calculate residuals
    df['residuals'] = df['logTF Communication Score'] - df['predicted_score']
    
    return df
