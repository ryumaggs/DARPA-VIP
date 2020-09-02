import torch
import torch.nn as nn
import torchvision.datasets as datasets
import torchvision.transforms as transforms
from torch.autograd import Variable
from torch.utils import data
from sklearn.model_selection import KFold
import scipy.io as scio
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.multiclass import OneVsRestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.metrics import accuracy_score
from sklearn.svm import LinearSVC
from skmultilearn.problem_transform import ClassifierChain
from sklearn.naive_bayes import GaussianNB
from skmultilearn.problem_transform import LabelPowerset
from sklearn.neural_network import MLPClassifier

musicnet_data = scio.loadmat('musicnet_low_dim.mat')
low_dim = musicnet_data['low_dim']


all_scales = {}
num_scale = low_dim[0][0].shape[1]
for i in range(num_scale):
    scale_i = []
    for j in range(low_dim.shape[1]):
        p = low_dim[0][j][0][i]
        if p.size != 0:
            p = p.reshape(p.shape[0])
        else:
            p = all_scales[i-1][j]
        scale_i.append(p)
    all_scales[i] = np.array(scale_i)
for i in range(7):
    scale_i_sizes = {}
#     size_scale_i = all_scales[i][0].shape[0]
    for j in range(low_dim.shape[1]):
        size_j = all_scales[i][j].shape[0]
        if not scale_i_sizes.get(size_j):
            scale_i_sizes[size_j] = 1
        else:
            scale_i_sizes[size_j] += 1
    #print(f"scale {i} num_sizes: {len(scale_i_sizes)}, distribution: {scale_i_sizes}")


max_size_for_scale = {}
for i in range(len(all_scales)):
    size_scale_i = all_scales[i][0].shape[0]
    for j in range(low_dim.shape[1]):
        size_j = all_scales[i][j].shape[0]
        size_scale_i = max(size_j, size_scale_i)
    max_size_for_scale[i] = size_scale_i
#print(max_size_for_scale)

for i in range(len(all_scales)):
    embedded_scale_i = []
    for j in range(low_dim.shape[1]):
        embedded_vect = all_scales[i][j]
        if all_scales[i][j].shape[0] < max_size_for_scale[i]:
            embedded_vect = np.zeros((max_size_for_scale[i],))
            embedded_vect[:all_scales[i][j].shape[0]] = all_scales[i][j]
        embedded_scale_i.append(embedded_vect)
    all_scales[i] = np.asarray(embedded_scale_i)   
np.save("musicnet_normalized_by_scale", all_scales, allow_pickle=True)


load_data = np.load('musicnet_normalized_by_scale.npy',allow_pickle=True)
Labels = np.load('arr_0.npy',allow_pickle=True)

# training network with rough scale representation (scale 0)
correct_train = 0
correct_test = 0
total_stdev_test = 0
total_stdev_train = 0
data_splits = 10

scale_i_data = np.matrix(load_data.item()[0])
kf = KFold(n_splits=data_splits, shuffle=True)
kf.get_n_splits(scale_i_data, Labels)
(train_acc, test_acc) = (0,0)
for train_index, test_index in kf.split(scale_i_data,Labels):
    X_train, X_test = scale_i_data[train_index], scale_i_data[test_index]
    y_train, y_test = Labels[train_index], Labels[test_index]
    
    #Initialize network
    classifier = MLPClassifier(solver='adam', alpha=1e-5,hidden_layer_sizes=(120,), random_state=1, max_iter=30000)
    classifier.fit(X_train, y_train)

    #predicting training data
    #predictions_train = classifier.predict(X_train)
    '''
    #accuracy check
    correct_split_train = 0
    stdev_train = 0
    pred_array_length_train = predictions_train.toarray().shape[0]
    for i in range(pred_array_length_train):
        correct_split_train = correct_split_train + accuracy_score(y_train[i],predictions_train.toarray()[i])
        stdev_train = stdev_train + np.std(np.array(y_train[i]) - np.array(predictions_train.toarray()[i]))

    total_stdev_train = total_stdev_train + stdev_train/pred_array_length_train
    avg_correct_train = 100*correct_split_train/pred_array_length_train
    correct_train = correct_train + avg_correct_train'''

    #predicting test data
    predictions_test = classifier.predict(X_test)
    print(predictions_test.shape)

    #accuracy check
    correct_split_test = 0
    stdev_test = 0
    pred_array_length = predictions_test.shape[0]
    for i in range(pred_array_length):
        accurate = accuracy_score(y_test[i],predictions_test[i])
        stdev = np.std(np.array(y_test[i]) - np.array(predictions_test[i]))
        correct_split_test = correct_split_test + accurate
        stdev_test = stdev_test + stdev
        print("Acc  ", accurate)
        print("Std  ", stdev)
        print("done")
    correct_test = correct_test + 100*correct_split_test/pred_array_length
    total_stdev_test = total_stdev_test + stdev_test/pred_array_length

'''print("Accuracy for train data is: ", correct_train/data_splits, "%")
print("Standard deviation for train data is: ", total_stdev_train/data_splits)'''
'''
print("Accuracy for test data is: ", correct_test/data_splits, "%")
print("Standard deviation for test data is: ", total_stdev_test/data_splits)'''