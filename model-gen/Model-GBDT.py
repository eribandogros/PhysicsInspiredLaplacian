import numpy as np
import argparse
import pandas as pd
from sklearn import svm
from sklearn import tree
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import metrics 
import scipy.stats as stats
import pickle
import os
import argparse

parser = argparse.ArgumentParser(description='GBDT model')
parser.add_argument('--dataset', type=str, default='../dataset')
args = parser.parse_args()
dataset=args.dataset
   
# Read data set in *.csv to data frame in Pandas
def read_dataset(feature_file, label_file):
    df_X = np.load(feature_file, allow_pickle=True)
    df_y = pd.read_csv(label_file, header=None, index_col=False)
    X = df_X # convert features in dataframe to numpy array
    y = df_y.values # convert label in dataframe to numpy array
    return X, y
    
ml_method="GradientBoostingClassifier"

def build_model():
    X_train, y_train = read_dataset('%s/train.npy'%(dataset), '%s/label_train.csv'%(dataset))
    y_train=np.ravel(y_train)

    X_test, y_test = read_dataset('%s/test.npy'%(dataset), '%s/label_test.csv'%(dataset))
    y_test=np.ravel(y_test)
    
    # Model parameters based on training set size
    num=np.shape(X_train)[0]
    if(num < 1000):
        i=10000; j=7; k=3; m=2
    elif(num >= 1000 and num < 10000):
        i=10000; j=8; k=4; m=3
    elif(num >= 10000):
        i=70000; j=5; k=9; m=9

    clf=globals()["%s"%ml_method](n_estimators=i, max_depth=j, min_samples_split=k, learning_rate=0.001, subsample=0.1*m, max_features='sqrt')
    
    clf.fit(X_train, y_train)
    res = clf.predict(X_test)

    save_model = False
    if save_model:
        model_name = 'model-GBDT.sav'
        pickle.dump(clf, open(model_name, 'wb'))

    res_positive_score = clf.predict_proba(X_test)

    mcc = metrics.matthews_corrcoef(y_test, res)
    accuracy = metrics.accuracy_score(y_test, res)
    recall = metrics.recall_score(y_test, res)
    precision = metrics.precision_score(y_test, res)
    F1_score = metrics.f1_score(y_test, res)
    roc_auc_score = metrics.roc_auc_score(y_test, res_positive_score[:,1])

    from sklearn.metrics import confusion_matrix
    conf_matrix = confusion_matrix(y_test, res)
    TP = conf_matrix[1][1]
    TN = conf_matrix[0][0]
    FP = conf_matrix[0][1]
    FN = conf_matrix[1][0]
    sensitivity = TP / float(TP + FN)
    specificity = TN / float(TN + FP)
    
    if not os.path.exists('results'):
        os.makedirs('results')
    
    f1 = open('results/prediction-GBDT.txt', 'w')
    print("mcc=%f accuracy=%f precision=%f roc_auc_score=%f sensitivity=%f specificity=%f"%(mcc,accuracy, precision, roc_auc_score, sensitivity, specificity), file=f1)
    f1.close()

build_model()
