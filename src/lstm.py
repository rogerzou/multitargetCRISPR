
import numpy as np
from scipy import stats
from . import mtss as m
from . import ml as ml


def series_to_supervised(f_h2ax, f_insulation, f_4Cseq):
    h2ax = m.load_nparray(f_h2ax)[1:, 2:]
    insu = m.load_nparray(f_insulation)[1:, 2:]
    fCsq = m.load_nparray(f_4Cseq)[1:, 2:]
    head = np.asarray(m.load_npheader(f_h2ax).split(',')[2:]).astype(float)
    datanone = np.logical_or(h2ax == 'None', insu == 'None')
    h2ax = h2ax[~datanone.any(axis=1)].astype(float)
    insu = insu[~datanone.any(axis=1)].astype(float)
    fCsq = fCsq[~datanone.any(axis=1)].astype(float)
    fCsq = np.log(fCsq)
    [num_samp, num_time] = h2ax.shape
    X = np.zeros((num_samp * (num_time - 4), 4))
    y = np.zeros(num_samp * (num_time - 4))
    for i in range(num_samp):
        for j in range(num_time - 4):
            insu_slope = stats.linregress(np.arange(0, 4), insu[i, j:j + 4])[0]
            fCsq_slope = stats.linregress(np.arange(0, 4), fCsq[i, j:j + 4])[0]
            h2ax_slope = stats.linregress(np.arange(0, 4), h2ax[i, j:j + 4])[0]
            X[i * (num_time - 4) + j, :] = [insu[i, j], fCsq[i, j], insu_slope, fCsq_slope]
            y[i * (num_time - 4) + j] = h2ax_slope
    X_train, X_test, y_train, y_test = ml.data_split(X, y)
    ml.NeuralNetworkTrainGridCV(X_train, y_train, 'lstm_model.sav')
    ml.ModelTest(X_test, y_test, 'lstm_model.sav')
