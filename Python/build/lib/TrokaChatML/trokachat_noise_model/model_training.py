# model_training.py
import lightgbm as lgb
from skopt import BayesSearchCV
from skopt.space import Real, Integer
from sklearn.metrics import mean_squared_error, make_scorer
from sklearn.model_selection import KFold
from joblib import Parallel, delayed
from tqdm import tqdm
import pandas as pd
import multiprocessing
import time
import warnings
from .data_processing import load_h5_batches, preprocess_datasets, preprocess_and_combine_datasets
import pickle

def create_lgb_dataset(combined_data):
    excluded_columns = ['logTF Communication Score', 'Communication Score', 'Uniqueness Score']
    feature_columns = [col for col in combined_data.columns if col not in excluded_columns]
    cat_indices = [i for i, col in enumerate(feature_columns) if combined_data[col].dtype.name == 'category']
    X = combined_data[feature_columns]
    for col in feature_columns:
        X[col] = X[col].astype('category')
    y = combined_data['logTF Communication Score']
    lgb_data = lgb.Dataset(X, label=y, categorical_feature=feature_columns, free_raw_data=False)
    return lgb_data, feature_columns

def perform_kfold_cv(train_groups, filename, n_trials=50, n_jobs=None, batch_size=100):
    warnings.filterwarnings("ignore")
    start_time = time.time()

    if n_jobs is None:
        num_cores = multiprocessing.cpu_count()
        n_jobs = num_cores - 1
    print(f"Using {n_jobs} CPU cores for parallel processing")
    
    def load_and_preprocess_wrapper(groups):
        data_dict = load_h5_batches(filename, groups)
        data_dict = preprocess_datasets(data_dict)
        return preprocess_and_combine_datasets(data_dict)

    print("Loading and preprocessing data in batches...")
    with Parallel(n_jobs=n_jobs) as parallel:
        combined_data_batches = parallel(delayed(load_and_preprocess_wrapper)(train_groups[start:start + batch_size]) for start in tqdm(range(0, len(train_groups), batch_size), desc="Batches"))

    combined_data = pd.concat(combined_data_batches, ignore_index=True)
    lgb_train_data, feature_names = create_lgb_dataset(combined_data)

    param_space = {
        'learning_rate': Real(0.01, 0.3, 'uniform'),
        'num_leaves': Integer(8, 1000),
        'colsample_bytree': Real(0.5, 1.0, 'uniform'),  # feature_fraction
        'subsample': Real(0.5, 1.0, 'uniform'),         # bagging_fraction
        'max_depth': Integer(3, 16),
        'min_child_samples': Integer(50, 500),          # min_data_in_leaf
        'min_split_gain': Real(0, 15),                  # min_gain_to_split
        'reg_lambda': Real(0.01, 100),                  # lambda_l2
    }

    bayes_search = BayesSearchCV(
        estimator=lgb.LGBMRegressor(objective='regression', metric='rmse', boosting_type='gbdt', n_jobs=n_jobs, verbose=-1),
        search_spaces=param_space,
        n_iter=n_trials,
        cv=KFold(n_splits=5),
        scoring=make_scorer(mean_squared_error, greater_is_better=False),
        n_jobs=n_jobs
    )
    
    print("Starting Bayesian optimization with cross-validation...")
    bayes_search.fit(lgb_train_data.data, lgb_train_data.label, categorical_feature=feature_names)

    best_params = bayes_search.best_params_
    final_model = bayes_search.best_estimator_

    elapsed_time = time.time() - start_time
    print(f"Total Time elapsed: {elapsed_time:.2f} seconds")

    return final_model, best_params

def save_model(model, filepath):
    with open(filepath, 'wb') as model_file:
        pickle.dump(model, model_file)
    print(f"Model saved as {filepath}")

def load_model(filepath):
    with open(filepath, 'rb') as model_file:
        return pickle.load(model_file)