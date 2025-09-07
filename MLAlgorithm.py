# ==========================================================
# Python ML Script for Gene-Based Classification (Multiple Classifiers)
# ==========================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    f1_score,
    roc_curve,
    confusion_matrix
)
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.inspection import permutation_importance

# ------------------------------
# 1. Load datasets
# ------------------------------
train_df = pd.read_csv("/Users/yangchengxuan/Downloads/Oncology Dataset/combined_dataset_train.csv")
test_df = pd.read_csv("/Users/yangchengxuan/Downloads/Oncology Dataset/combined_dataset_valid.csv")

# Drop first column (sample ID)
X_train = train_df.iloc[:, 1:-1]
y_train = train_df.iloc[:, -1]

X_test = test_df.iloc[:, 1:-1]
y_test = test_df.iloc[:, -1]

# print(np.unique(y_train, return_counts=True))
# print(np.unique(y_test, return_counts=True))
# print(train_df["Label"])


# Encode labels if needed
if y_train.dtype == 'object':
    le = LabelEncoder()
    y_train = le.fit_transform(y_train)
    y_test = le.transform(y_test)

# ------------------------------
# 2. Create folder to save figures
# ------------------------------
os.makedirs("figures", exist_ok=True)

# ------------------------------
# 3. Define classifiers
# ------------------------------
classifiers = {
    "RandomForest": RandomForestClassifier(n_estimators=500, random_state=42),
    # "GradientBoosting": GradientBoostingClassifier(n_estimators=500, random_state=42),
    "LogisticRegression": LogisticRegression(max_iter=1000, random_state=42),
    "SVM_linear": SVC(kernel='linear', probability=True, random_state=42),
    # "SVM_rbf": SVC(kernel='rbf', probability=True, random_state=42),
    # "KNN": KNeighborsClassifier(),
    "DecisionTree": DecisionTreeClassifier(random_state=42)
}

# ------------------------------
# 4. Train, evaluate, and visualize
# ------------------------------
for name, clf in classifiers.items():
    print(f"\n--- {name} ---")
    
    # Train
    clf.fit(X_train, y_train)
    
    # Predictions
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test) if hasattr(clf, "predict_proba") else None
    
    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Multi-class ROC AUC calculation
    if y_proba is not None:
        # Check if it's multi-class
        n_classes = len(np.unique(y_test))
        if n_classes > 2:
            # Multi-class using 'ovr' (one-vs-rest) strategy
            roc_auc = roc_auc_score(y_test, y_proba, multi_class='ovr', average='weighted')
        else:
            # Binary classification
            roc_auc = roc_auc_score(y_test, y_proba[:, 1])
    else:
        roc_auc = "N/A"
    
    print(f"Accuracy: {accuracy:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"ROC-AUC: {roc_auc}")
    
    # Confusion Matrix
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(5,4))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title(f"{name} Confusion Matrix")
    plt.savefig(f"figures/{name}_confusion_matrix.png", dpi=300)
    plt.close()
    
    # ROC Curve (Multi-class support)
    if y_proba is not None:
        n_classes = len(np.unique(y_test))
        plt.figure(figsize=(8,6))
        
        if n_classes > 2:
            # Multi-class: plot ROC curve for each class
            from sklearn.preprocessing import label_binarize
            from sklearn.metrics import roc_curve, auc
            
            # Binarize labels
            y_test_bin = label_binarize(y_test, classes=np.unique(y_test))
            
            # Calculate ROC curve for each class
            colors = ['blue', 'red', 'green', 'orange', 'purple']
            class_names = ['Normal/Non-tumor', 'HBV-HCC', ' HCC']
            
            for i in range(n_classes):
                fpr, tpr, _ = roc_curve(y_test_bin[:, i], y_proba[:, i])
                roc_auc_class = auc(fpr, tpr)
                plt.plot(fpr, tpr, color=colors[i], lw=2,
                        label=f'{class_names[i]} (AUC = {roc_auc_class:.2f})')
            
            plt.plot([0, 1], [0, 1], 'k--', lw=2, label='Random Classifier')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'{name} Multi-class ROC Curves')
            plt.legend(loc="lower right")
        else:
            # Binary classification
            fpr, tpr, _ = roc_curve(y_test, y_proba[:, 1])
            plt.plot(fpr, tpr, label=f'ROC Curve (AUC={roc_auc:.2f})')
            plt.plot([0,1],[0,1],'k--')
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"{name} ROC Curve")
            plt.legend()
        
        plt.savefig(f"figures/{name}_roc_curve.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # Feature importance
    if hasattr(clf, "feature_importances_"):  # Tree-based
        feat_imp = pd.DataFrame({
            "Feature": X_train.columns,
            "Importance": clf.feature_importances_
        }).sort_values(by="Importance", ascending=False)
    elif hasattr(clf, "coef_"):  # Linear models
        if len(clf.coef_.shape) > 1:  # multi-class
            coef_abs = np.mean(np.abs(clf.coef_), axis=0)
        else:
            coef_abs = np.abs(clf.coef_)
        feat_imp = pd.DataFrame({
            "Feature": X_train.columns,
            "Importance": coef_abs
        }).sort_values(by="Importance", ascending=False)
    else:  # Non-linear SVM, KNN
        result = permutation_importance(clf, X_test, y_test, n_repeats=3, random_state=42, n_jobs=-1)
        feat_imp = pd.DataFrame({
            "Feature": X_train.columns,
            "Importance": result.importances_mean
        }).sort_values(by="Importance", ascending=False)
    
    print("Top 20 predictive features:")
    print(feat_imp.head(20))
    
    plt.figure(figsize=(8,6))
    sns.barplot(x='Importance', y='Feature', data=feat_imp.head(20), palette='viridis')
    plt.title(f"{name} Top 20 Feature Importances")
    plt.savefig(f"figures/{name}_feature_importance.png", dpi=300)
    plt.close()

# ------------------------------
# 5. Correlation analysis
# ------------------------------
feature_corr = X_train.copy()
feature_corr['Label'] = y_train
correlations = feature_corr.corr()['Label'].drop('Label').sort_values(key=abs, ascending=False)
print("\nTop correlated features with Label:")
print(correlations.head(10))
