import os
import numpy as np
import pandas as pd
import xgboost as xgb
import pickle
import datetime
import matplotlib.pyplot as plt
from smart_open import smart_open
from multiprocessing import cpu_count, Pool
from xgboost import plot_importance
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import confusion_matrix, mean_squared_error
from sklearn.metrics import average_precision_score, precision_recall_curve
from sklearn.utils.fixes import signature

from src.utility import load_data_clinical
from src.utility import load_data_RNASeq


NO_TREES = 200
MAX_DEPTH = 3
MIN_CHILD_WEIGHT = 5
G = 0.5
SUBSAMPLE = 0.6
NO_JOBS = cpu_count() * 3
MODEL_PATH = './models/models_xgboost/'
MODEL_LIST_PATH = './models/models_xgboost/model_list.txt'
FEATURE_IMP_PATH = './results/feature_importance_xgbt.txt'
PRECISION_RECALL_CURVE_PATH = './results/curve_xgboost.png'


def save_model(xgbt, xgbt_name):
    # para xgbt: XGBoost model to be saved
    # para xgbt_name: name of the XGBoost model
    model_name = xgbt_name
    if os.path.isdir(MODEL_PATH + model_name):
        print("\nmodel %s already existed" % model_name)
    else:
        os.mkdir(MODEL_PATH + model_name)
        with smart_open(MODEL_PATH + model_name + '/model.sav', 'wb') as save_path:
            pickle.dump(xgbt, save_path)
        # readme file
        readme_notes = np.array(["This %s model is trained on %s" % (model_name, str(datetime.datetime.now()))])
        np.savetxt(MODEL_PATH + model_name + '/readme.txt', readme_notes, fmt="%s")
        # add model names to the list for further load
        with smart_open(MODEL_PATH + 'model_list.txt', 'a+') as f:
            f.write(MODEL_PATH + model_name + '\n')


def load_model(model_no):
    # para model_no: which XGBoost model to load
    with smart_open(MODEL_LIST_PATH, 'rb', encoding='utf-8') as model_list:
        for line_no, line in enumerate(model_list):
            if line_no == model_no - 1:
                model_path = str(line).replace('\n', '')
                with smart_open(model_path + '/model.sav', 'rb') as f:
                    xgbt = pickle.load(f)
                break
    return xgbt


def save_importances(save2file):
    # para save2file: list(str)
    with smart_open(FEATURE_IMP_PATH, 'w+', encoding='utf-8') as f:
        for i in range(len(save2file)):
            f.write(save2file[i])
            f.write("\n")


'''ONLY SHOW TOP 50 FEATURES'''
def show_important_feature(xgbt, data, save=True, img=False):
    # para xgbt: XGBoost model to draw important features
    # para data: RNA_Seq data w/ index
    # para save: whether or not to save importances to file
    # para img: whether or not to show the image
    feature_names = list(data.columns.values)
    importances = xgbt.feature_importances_
    print("\nFeature ranking:")
    indices = np.argsort(importances)[::-1]
    save2file = []
    feature_names_img, feature_importance_img = [], []
    for f in range(50):
        print("%d. feature %d %s (%f)" % (f+1, indices[f], feature_names[indices[f]], importances[indices[f]]))
        if save:
            save2file.append("%d. feature %d %s (%f)" % (f+1, indices[f], feature_names[indices[f]], importances[indices[f]]))
        if img:
            feature_names_img.append(feature_names[indices[f]])
            feature_importance_img.append(importances[indices[f]])
    
    if save:
        save_importances(save2file)
    if img:
        plt.barh(range(50), feature_importance_img, align='center')
        plt.yticks(np.arange(50), feature_names_img)
        plt.xlabel("Feature Importance")
        plt.ylabel("Feature")
        plt.show()


def run_xgboost_classifier(load=False, model_no=1):
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
        print("\nload pre-trained no.%d model" % model_no)
        xgb_model = load_model(model_no)
    else:
        print("\ntraining a XGBoost classifier ...")
        xgb_model = xgb.XGBClassifier(min_child_weight=MIN_CHILD_WEIGHT, gamma=G, subsample=SUBSAMPLE, n_estimators=NO_TREES, max_depth=MAX_DEPTH)
        xgbt_name = "min_child_weight=%s,gamma=%s,subample=%s,n_estimators=%s,max_depth=%s" % (str(MIN_CHILD_WEIGHT), str(G), str(SUBSAMPLE), str(NO_TREES), str(MAX_DEPTH))

        xgb_model.fit(X_train, y_train)

        print("\ntraining DONE. \n\nsaving the XGBoost classifer ...")
        save_model(xgb_model, xgbt_name)
        
    print("\ntesting the XGBoost classifier ...")
    y_pred = xgb_model.predict(X_test)
    print(mean_squared_error(y_test, y_pred))
    print(confusion_matrix(y_test, y_pred))

    # show & save the top 50 important features
    show_important_feature(xgb_model, data_RNASeq, img=False)
    # draw the precision recall curve for the classifier
    draw_precision_recall_curve(y_test, y_pred)


def draw_precision_recall_curve(y_test, y_score):
    fig = plt.figure()
    average_precision = average_precision_score(y_test, y_score)
    precision, recall, _ = precision_recall_curve(y_test, y_score)
    step_kwargs = ({'step': 'post'} if 'step' in signature(plt.fill_between).parameters else {})
    plt.step(recall, precision, color='b', alpha=0.2, where='post')
    plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision))
    fig.savefig(PRECISION_RECALL_CURVE_PATH, dpi=300)


'''ONLY EXECUTE ONCE'''
def tune_hyperparameters():
    # loda the data
    data_RNASeq_labels = load_data_RNASeq()
    data_RNASeq_labels = data_RNASeq_labels.drop(columns=['gene'])

    data_labels = data_RNASeq_labels['label']
    data_RNASeq = data_RNASeq_labels.drop(columns=['label'])

    # train/test split 
    print("\nsplitting the training/test dataset ...")
    X_train, X_test, y_train, y_test = train_test_split(data_RNASeq, data_labels)

    params = {
        "min_child_weight": [1, 5 ,10],
        "gamma": [0.5, 1, 2],
        "subsample": [0.6, 0.8, 1.0],
        "max_depth": [3, 5],
        "n_estimators": [50, 200]
    }

    print("\nrunning the Grid Search for XGBoost classifier ...")
    clf = GridSearchCV(xgb.XGBClassifier(), params, cv=3, n_jobs=NO_JOBS, verbose=10)

    clf.fit(X_train, y_train)
    print(clf.best_score_)
    print(clf.best_estimator_)

    # save tuned hyperparameters
    with smart_open("./results/best_params_xgboost.txt", 'w', encoding='utf-8') as f:
        f.write(str(clf.best_params_) + str(clf.best_score_))
    print("\nbest hyperparameters for XGBOOST has been written to file.")

