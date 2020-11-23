from math import sqrt
from numpy import concatenate
from matplotlib import pyplot
from pandas import read_csv
from pandas import DataFrame
from pandas import concat
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import mean_squared_error
from keras.models import Sequential
from keras.layers import Dense, Dropout, LSTM
import tensorflow as tf
import numpy as np
from . import mtss as m


def series_to_supervised(f_h2ax, f_insulation):
    h2ax = m.load_nparray(f_h2ax)[1:, 2:]
    insu = m.load_nparray(f_insulation)[1:, 2:]
    head = np.asarray(m.load_npheader(f_h2ax).split(',')[2:]).astype(float)
    datanone = np.logical_or(h2ax == 'None', insu == 'None')
    h2ax = h2ax[~datanone.any(axis=1)].astype(float)
    insu = insu[~datanone.any(axis=1)].astype(float)
    # head = np.repeat(np.transpose(head.reshape(-1, 1)), repeats=num_samp, axis=0)
    X = np.expand_dims(insu, axis=2)
    y = h2ax
    [num_samp, num_time, num_feat] = X.shape
    model = Sequential()
    model.add(LSTM(units=50, return_sequences=True, input_shape=(num_time, num_feat)))
    model.add(Dropout(0.2))
    model.add(LSTM(units=50, return_sequences=True))
    model.add(Dropout(0.2))
    model.add(Dense(1))
    # compile model
    model.compile(loss='mse', optimizer='adam')
    # fit model
    model.fit(X, y, epochs=100, batch_size=32)
    model.save('lstm_model.h5')
    yhat = np.squeeze(model.predict(X))
    np.savetxt("lstm_yhat.csv", yhat, delimiter=',')
    np.savetxt("lstm_y.csv", y, delimiter=',')
    print(model)



