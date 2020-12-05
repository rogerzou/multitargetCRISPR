
import numpy as np
from scipy import stats
from . import mtss as m
import sklearn
import pickle


def calc_correlation_with_y(xyfile, head):
    X, y = load_Xy_matrix(xyfile)
    Xpos = X[y > 0, :]
    Xneg = X[y <= 0, :]
    ttest = stats.ttest_ind(Xpos, Xneg)
    Xpos_avg = np.mean(Xpos, axis=0)
    Xneg_avg = np.mean(Xneg, axis=0)
    np.savetxt(xyfile + "ttest.csv", np.vstack((Xpos_avg, Xneg_avg, ttest.statistic, ttest.pvalue)),
               fmt='%s', delimiter=',', header=head)


def save_Xy_matrix(y_file, X_files, outfile):
    y_data = m.load_nparray(y_file)[:, 2:]
    X_data = []
    datanone = y_data == 'None'
    for x_file in X_files:
        x = m.load_nparray(x_file)[:, 2:]
        X_data.append(x)
        datanone += (x == 'None')
    y_data = y_data[~datanone.any(axis=1)].astype(float)
    [num_samp, num_time] = y_data.shape
    num_input = len(X_data)
    for k in range(num_input):
        X_data[k] = X_data[k][~datanone.any(axis=1)].astype(float)
    X = np.zeros((num_samp * (num_time - 4), num_input * 2))
    Xalt = []
    y = np.zeros(num_samp * (num_time - 4))
    for j in range(num_time - 4):
        Xalt_i = np.zeros((num_samp, num_input * 2))
        for i in range(num_samp):
            x_vals = []
            for k in range(num_input):
                x_vals.append(X_data[k][i, j])
                x_vals.append(stats.linregress(np.arange(0, 4), X_data[k][i, j:j + 4])[0])
            y_slope = stats.linregress(np.arange(0, 4), y_data[i, j:j + 4])[0]
            X[i * (num_time - 4) + j, :] = x_vals
            y[i * (num_time - 4) + j] = y_slope
            Xalt_i[i, :] = x_vals
        Xalt.append(Xalt_i)
    np.savetxt(outfile + ".csv", np.hstack((y.reshape(-1, 1), X)), fmt='%s', delimiter=',')
    pickle.dump((X, y, Xalt), open(outfile + ".pickle", 'wb'))
    return X, y, Xalt


def remove_outliers(X, y, outfile):
    avg_y, std_y = np.average(y), np.std(y)
    y[y > avg_y + 3 * std_y] = avg_y + 3 * std_y
    y[y < avg_y - 3 * std_y] = avg_y - 3 * std_y
    for i in range(X.shape[1]):
        x_i = X[:, i]
        avg_i, std_i = np.average(x_i), np.std(x_i)
        x_i[x_i > avg_i + 3 * std_i] = avg_i + 3 * std_i
        x_i[x_i < avg_i - 3 * std_i] = avg_i - 3 * std_i
        X[:, i] = x_i
    np.savetxt(outfile + ".csv", np.hstack((y.reshape(-1, 1), X)), fmt='%s', delimiter=',')
    return X, y


def load_Xy_matrix(infile):
    data = np.loadtxt(infile, delimiter=',').astype(float)
    return data[:, 1:], data[:, 0]


def modify_matrix(X, y, classifier=False, normalize=False):
    X = sklearn.preprocessing.normalize(X) if normalize else X
    y = (y > 0).astype(int) if classifier else y
    return X, y


def from_cutsite_predict_slope(modelfile, initial0, X, outpath):
    model = pickle.load(open(modelfile, 'rb'))
    num_iters = len(X)
    num_sites, num_input = X[0].shape
    outnpy = np.zeros((num_sites, num_iters + 1))
    outnpy[:, 0] = initial0
    for i in range(num_iters):
        outnpy[:, i+1] = outnpy[:, i] + model.predict(X[i])
    np.savetxt(outpath + ".csv", np.asarray(outnpy), fmt='%s', delimiter=',')
