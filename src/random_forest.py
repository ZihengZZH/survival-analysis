import os
import numpy as np
import pandas as pd
import pickle
import datetime
import matplotlib.pyplot as plt
from smart_open import smart_open
from multiprocessing import cpu_count, Pool
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.utils.fixes import signature
from sklearn.model_selection import GridSearchCV

from src.utility import load_data_clinical
from src.utility import load_data_RNASeq


NO_TREES = 500
MAX_FEATURES = 0.1
CRITERION = 'entropy'
NO_JOBS = cpu_count() * 3
MODEL_PATH = './models/models_random_forest/'
MODEL_LIST_PATH = './models/models_random_forest/model_list.txt'
FEATURE_IMP_PATH = './results/feature_importance_rf.txt'
ROC_CURVE_PATH = './results/roc_curve_random_forest.png'


def save_model(forest, forest_name):
    # para forest: RF model to be saved
    # para forest_name: name of the RF model
    model_name = forest_name
    if os.path.isdir(MODEL_PATH + model_name):
        print("\nmodel %s already existed" % model_name)
    else:
        os.mkdir(MODEL_PATH + model_name)
        with smart_open(MODEL_PATH + model_name + '/model.sav', 'wb') as save_path:
            pickle.dump(forest, save_path)
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
                    forest = pickle.load(f)
                break
    return forest


def save_importances(save2file):
    # para save2file: list(str)
    with smart_open(FEATURE_IMP_PATH, 'w+', encoding='utf-8') as f:
        for i in range(len(save2file)):
            f.write(save2file[i])
            f.write("\n")


def show_important_feature(forest, data, save=True, img=False):
    # para forest: RF model to draw important features
    # para data: RNA_Seq data w/ index
    # para save: whether or not to save importances to file
    # para img: whether or not to show the image
    n_features = data.shape[1]
    feature_names = list(data.columns.values)
    importances = forest.feature_importances_
    print("\nFeature ranking:")
    indices = np.argsort(importances)[::-1]
    save2file = []
    feature_names_img, feature_importance_img = [], []
    # only output top 50 important features
    for f in range(50):
        print("%d. feature %d %s (%f)" % (f+1, indices[f], feature_names[indices[f]], importances[indices[f]]))
        if save:
            save2file.append("%d. feature\t %d\t %s\t (%f)" % (f+1, indices[f], feature_names[indices[f]], importances[indices[f]]))
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


def run_random_forest(load=False, model_no=1):
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
        forest = load_model(model_no)
    else:
        print("\ntraining a Random Forest classifier ...")
        forest = RandomForestClassifier(n_estimators=NO_TREES, random_state=0, max_features=MAX_FEATURES, criterion=CRITERION, n_jobs=NO_JOBS)
        forest_name = "n_estimators=%s,max_features=%s,criterion=%s,n_jobs=%s" % (NO_TREES, MAX_FEATURES, CRITERION, NO_JOBS)

        forest.fit(X_train, y_train)

        print("\ntraining DONE.\n\nsaving the RF classifier ...")
        save_model(forest, forest_name)

    print("\ntesting the Random Forest classifier ...\n")
    y_pred = forest.predict(X_test)
    print("Accuracy on training set: %.3f" % forest.score(X_train, y_train))
    print("Accuracy on test set: %.3f" % forest.score(X_test, y_test))

    show_important_feature(forest, data_RNASeq, img=False)
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
    fig.savefig(ROC_CURVE_PATH, dpi=300)


'''ONLY EXECUTE ONCE'''
def tune_hyperparameters():
    # load the data
    data_RNASeq_labels = load_data_RNASeq()
    data_RNASeq_labels = data_RNASeq_labels.drop(columns=['gene'])

    data_labels = data_RNASeq_labels['label']
    data_RNASeq = data_RNASeq_labels.drop(columns=['label'])

    # train/test split 
    print("\nsplitting the training/test dataset ...")
    X_train, X_test, y_train, y_test = train_test_split(data_RNASeq, data_labels)

    parameters = {
        "n_estimators": [100, 200, 300, 400, 500],
        "max_features": [0.1, 0.2, 0.25, 0.3, 0.4, 0.5],
        "criterion": ["entropy"]
    }

    print("\nrunning the Grid Search for Random Forest classifier ...")
    clf = GridSearchCV(RandomForestClassifier(), parameters, cv=3, n_jobs=NO_JOBS, verbose=10)

    clf.fit(X_train, y_train)
    print(clf.score(X_train, y_train))
    print(clf.best_params_)

    # np.savetxt("./results/all_params_random_forest.txt", clf.cv_results_)
    with smart_open("./results/best_params_random_forest.txt", 'w', encoding='utf-8') as f:
        f.write(str(clf.best_params_) + str(clf.best_score_))
    print("\nbest hyperparameters for RF has been written to file.")