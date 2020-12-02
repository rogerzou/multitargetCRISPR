# -*- coding: utf-8 -*-
""" Machine learning analysis of ChIP-seq after multi-targeting Cas9
"""
__author__ = "Roger Zou"
__license__ = "MIT"
__version__ = "0.9"
__maintainer__ = "Roger Zou"

from . import mtss as m
import numpy as np
from scipy import sparse
import pickle
from sklearn import svm
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.inspection import permutation_importance
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, roc_auc_score
from matplotlib import pyplot as pp


def getXy_all(data, yindex, epi=True, mm=2):
    fl = m.load_nparray(data)
    head = m.load_npheader(data).split(', ')
    return _getXy_helper(fl, head, epi, mm, yindex)


def getXy_2orLess(data, yindex, epi=True, mm=2):
    fl = m.load_nparray(data)
    head = m.load_npheader(data).split(', ')
    fl = fl[fl[:, head.index('mismatches')].astype(int) <= 2, :]
    return _getXy_helper(fl, head, epi, mm, yindex)


def getXy_noMM(data, yindex, epi=True, mm=2):
    fl = m.load_nparray(data)
    head = m.load_npheader(data).split(', ')
    fl = fl[fl[:, head.index('mismatches')] == '0', :]   # get non-mismatched columns only
    return _getXy_helper(fl, head, epi, mm, yindex)


def _getXy_helper(fl, head, epi, emm, yindex):
    y = fl[:, head.index(yindex)].astype(float)
    if emm == 2:     # one-hot encoding of exact mismatch
        obser_ind = head.index('observed target sequence')
        expec_ind = head.index('expected target sequence')
        onehot_mm = []
        for i in range(fl.shape[0]):
            obs_i = fl[i, obser_ind]
            exp_i = fl[i, expec_ind]
            onehot_mm.append([x + 1 if obs_i[j] != exp_i[j] else x for j, x in enumerate([0] * len(exp_i))])
        onehot_mm = np.asarray(onehot_mm)
        label_mm = ["mm_%i" % (20 - x) for x in range(20)]
    elif emm == 1:       # one-hot encoding of mismatch type
        index_mm = head.index('mm_type')
        int_mm = LabelEncoder().fit_transform(fl[:, index_mm])
        onehot_mm = OneHotEncoder(sparse=False).fit_transform(int_mm.reshape(-1, 1))
        label_mm = ["mm_%i" % x for x in range(onehot_mm.shape[1])]
    else:
        onehot_mm = []
        label_mm = []
    # one-hot encoding of chromatin state
    if epi:
        index_epista, index_epiend = head.index('h3k4me1'), head.index('rna')
        epigen = fl[:, index_epista:index_epiend].astype(float)
        label_epi = head[index_epista:index_epiend]
        if emm >= 1:
            X = np.column_stack((onehot_mm, epigen))
            labels = label_mm + list(label_epi)
        else:
            X = epigen
            labels = list(label_epi)
    else:
        if emm >= 1:
            X = onehot_mm
            labels = label_mm
        else:
            return ValueError("Empty data features for model!")
    return X, y, labels


def pca(X, num_components=6):
    pca = PCA(n_components=num_components)
    X = pca.fit_transform(X)
    return X


def SVMTrainDefault(X, y, modelfile, classifier=False):
    if classifier:
        clf = svm.SVC(gamma='auto', probability=True)
    else:
        clf = svm.SVR(gamma='auto')
    clf.fit(X, y)
    pickle.dump(clf, open(modelfile, 'wb'))
    return clf


def LassoTrainDefault(X, y, modelfile):
    X_sp = sparse.coo_matrix(X)
    sparse_lasso = Lasso(alpha=1, fit_intercept=False, max_iter=1000)
    sparse_lasso.fit(X_sp, y)
    pickle.dump(sparse_lasso, open(modelfile, 'wb'))
    return sparse_lasso


def LinearRegressionTrainDefault(X, y, modelfile):
    reg = LinearRegression().fit(X, y)
    pickle.dump(reg, open(modelfile, 'wb'))
    return reg


def NeuralNetworkTrainDefault(X, y, modelfile, classifier=False):
    # Define and fit base model
    params = {'solver': 'lbfgs', 'max_iter': 5000, 'random_state': 42,
              'hidden_layer_sizes': (8,), 'alpha': 0.1}
    mlp = MLPClassifier(**params) if classifier else MLPRegressor(**params)
    mlp.fit(X, y)
    evaluate(X, y, mlp, classifier)
    # Save and return base estimator
    pickle.dump(mlp, open(modelfile, 'wb'))
    return mlp


def NeuralNetworkTrainGridCV(X, y, modelfile, classifier=False):
    # Define base model
    params = {'random_state': 42, 'max_iter': 5000}
    mlp = MLPClassifier(**params) if classifier else MLPRegressor(**params)
    # Find best hyperparameters and then fit on all training data
    hyperparams = {
        'alpha': [0.00001, 0.0001, 0.001, 0.01, 0.1, 1],
        'hidden_layer_sizes': [(10, 4), (10, 10, 10), (40, 20, 5), (100, 50, 10), (40, 20, 10, 5)],
        'solver': ['lbfgs', 'adam']
    }
    mlp_grid = GridSearchCV(estimator=mlp, param_grid=hyperparams, cv=5, verbose=2, n_jobs=12)
    mlp_grid.fit(X, y)
    # Get the best estimator and calculate results (correlation coefficient)
    best_grid = mlp_grid.best_estimator_
    evaluate(X, y, best_grid, classifier)
    # Save and return best estimator
    pickle.dump(best_grid, open(modelfile, 'wb'))
    return best_grid


def RandomForestTrainDefault(X, y, modelfile, classifier=False):
    # Define and fit base model:
    params = {'random_state': 42}
    rf = RandomForestClassifier(**params) if classifier else RandomForestRegressor(**params)
    rf.fit(X, y)
    # Calculate results (correlation coefficient)
    evaluate(X, y, rf, classifier)
    # Save and return base estimator
    pickle.dump(rf, open(modelfile, 'wb'))
    return rf


def RandomForestTrainGridCV(X, y, modelfile, classifier=False):
    # Define base model:
    params = {'random_state': 42}
    rf = RandomForestClassifier(**params) if classifier else RandomForestRegressor(**params)
    # Define search space:
    hyperparams = {
        'n_estimators': [1400],
        'max_depth': [int(x) for x in np.linspace(10, 110, num=3)],
        'min_samples_split': [5, 10],
        'min_samples_leaf': [2, 4],
        'bootstrap': [True, False]
    }
    # Find best hyperparameters and then fit on all training data
    rf_random = GridSearchCV(estimator=rf, param_grid=hyperparams, cv=5, verbose=2, n_jobs=4)
    rf_random.fit(X, y)
    # Get the best estimator and calculate results (correlation coefficient)
    best_random = rf_random.best_estimator_
    evaluate(X, y, best_random, False)
    # Save and return best estimator
    pickle.dump(best_random, open(modelfile, 'wb'))
    return best_random


def ModelTest(X, y, modelfile, classifier=False):
    model = pickle.load(open(modelfile, 'rb'))
    print(modelfile)
    print(model)
    evaluate(X, y, model, classifier)
    np.savetxt(modelfile[0:-4] + "_out.csv",
               np.column_stack((y, model.predict(X))), fmt='%s', delimiter=',')


def evaluate(X_test, y_test, model, classifier):
    y_pred = model.predict(X_test)
    if classifier:
        y_prob = model.predict_proba(X_test)
        print(confusion_matrix(y_test, y_pred))
        AUC_ROC(y_test, y_prob)
    else:
        CORREL(y_test, y_pred)


def FeatureImportance(X, y, modelfile, labels, count=10):
    # load model and labels of variables
    model = pickle.load(open(modelfile, 'rb'))
    labels = np.asarray(labels, dtype=object)
    # calculate feature importance, sorted by increasing importance
    result = permutation_importance(model, X, y, n_repeats=10, random_state=42, n_jobs=8)
    sorted_idx = result.importances_mean.argsort()
    # generate box plot of permutation importance
    fig, ax = pp.subplots()
    importances_i = result.importances[sorted_idx].T[:, -count:]
    labels_i = labels[sorted_idx][-count:]
    ax.boxplot(importances_i, vert=False, labels=labels_i)
    ax.set_title("Permutation Importances (test set)")
    fig.tight_layout()
    fig.savefig(modelfile[0:-4] + "_fi.png")
    pp.close(fig)
    # save raw permutation importances to CSV file
    np.savetxt(modelfile[0:-4] + "_fi.csv", importances_i, fmt='%s', delimiter=',', header=",".join(labels_i))


def data_split(X, y, test_size=0.3):
    return train_test_split(X, y, test_size=test_size, shuffle=True, random_state=42)


def MAPE(y_test, y_pred):
    errors = abs(y_pred - y_test)
    mape = 100 * np.mean(errors / y_test)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))
    return accuracy


def CORREL(y_test, y_pred):
    corr = np.corrcoef(y_test, y_pred)
    print('Model Performance')
    print('Correlation coefficient = {:0.2f}%.'.format(corr[1, 0]))
    return


def AUC_ROC(y_test, y_prob):
    macro_roc_auc_ovo = roc_auc_score(y_test, y_prob, multi_class="ovo", average="macro")
    weighted_roc_auc_ovo = roc_auc_score(y_test, y_prob, multi_class="ovo", average="weighted")
    macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
    weighted_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="weighted")
    print("One-vs-One ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
          "(weighted by prevalence)"
          .format(macro_roc_auc_ovo, weighted_roc_auc_ovo))
    print("One-vs-Rest ROC AUC scores:\n{:.6f} (macro),\n{:.6f} "
          "(weighted by prevalence)"
          .format(macro_roc_auc_ovr, weighted_roc_auc_ovr))