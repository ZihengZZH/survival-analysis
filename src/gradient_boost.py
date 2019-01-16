import os
import numpy as np
import pandas as pd
import pickle
import datetime
import matplotlib.pyplot as plt
from smart_open import smart_open
from multiprocessing import cpu_count, Pool
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split

from src.utility import load_data_clinical
from src.utility import load_data_RNASeq


NO_TREES = 100
MAX_DEPTH = 3
MAX_FEATURES = 'log2'
LEARNING_RATE = 0.1
LOSS = 'deviance'
MODEL_PATH = './models/models_gradient_boost/'
MODEL_LIST_PATH = './models/models_gradient_boost/model_list.txt'


def save_model(gbrt, gbrt_name):
    # para gbrt: GB model to be saved
    # para gbrt_name: name of the GB model
    model_name = gbrt_name
    if os.path.isdir(MODEL_PATH + model_name):
        print("\nmodel %s already existed" % model_name)
    else:
        os.mkdir(MODEL_PATH + model_name)
        with smart_open(MODEL_PATH + model_name + '/model.sav', 'wb') as save_path:
            pickle.dump(gbrt, save_path)
        # readme file
        readme_notes = np.array(["This %s model is trained on %s" % (model_name, str(datetime.datetime.now()))])
        np.savetxt(MODEL_PATH + model_name + '/readme.txt', readme_notes, fmt="%s")
        # add model names to the list for further load
        with smart_open(MODEL_PATH + 'model_list.txt', 'a+') as f:
            f.write(MODEL_PATH + model_name + '\n')


def load_model(model_no):
    # para model_no:
    with smart_open(MODEL_LIST_PATH, 'rb', encoding='utf-8') as model_list:
        for line_no, line in enumerate(model_list):
            if line_no == model_no - 1:
                model_path = str(line).replace('\n', '')
                with smart_open(model_path + '/model.sav', 'rb') as f:
                    gbrt = pickle.load(f)
                break
    return gbrt


def run_gradient_boost(load=False, model_no=1):
    # para load: whether or not to load pre-trained model
    # para model_no: if load, which model to load
    data_RNASeq_labels = load_data_RNASeq()
    data_RNASeq_labels = data_RNASeq_labels.drop(columns=['gene'])

    data_labels = data_RNASeq_labels['label']
    data_RNASeq = data_RNASeq_labels.drop(columns=['label'])

    # train/test split 
    print("\nsplitting the training/test dataset ...")
    X_train, X_test, y_train, y_test = train_test_split(data_RNASeq, data_labels)

    if load:
        gbrt = load_model(model_no)
    else:
        print("\ntraining a Gradient Boosting Tree classifier ...")
        gbrt = GradientBoostingClassifier(n_estimators=NO_TREES, random_state=0, max_features=MAX_FEATURES, max_depth=MAX_DEPTH, learning_rate=LEARNING_RATE)
        gbrt_name = "n_estimators=%s,max_features=%s,max_depth=%s,learning_rate=%s" % (str(NO_TREES), MAX_FEATURES, str(MAX_DEPTH), str(LEARNING_RATE))

        gbrt.fit(X_train, y_train)

        print("\ntraining DONE.\n\nsaving the GB classifier ...")
        save_model(gbrt, gbrt_name)


    print("\ntesting the Gradient Boosting Tree classifier ...\n")
    print("Accuracy on training set: %.3f" % gbrt.score(X_train, y_train))
    print("Accuracy on test set: %.3f" % gbrt.score(X_test, y_test))
