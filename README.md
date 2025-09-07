# Oncology Dataset Analysis Pipeline

A comprehensive machine learning pipeline for analyzing hepatocellular carcinoma (HCC) gene expression data from multiple GEO datasets.

## 📊 Datasets Used

This project utilizes three Gene Expression Omnibus (GEO) datasets:

- **GSE14520**: Hepatocellular carcinoma expression profiles
- **GSE25097**: HCC tumor vs. non-tumor liver tissue comparison (268 HCC tumor, 243 adjacent non-tumor, 40 cirrhotic, and 6 healthy liver samples)
- **GSE121248**: Additional HCC expression data

## 🔧 Prerequisites

### R Dependencies
```r
install.packages(c("GEOquery"))
```

### Python Dependencies
```bash
pip install pandas numpy matplotlib seaborn scikit-learn
```

## 📁 Project Structure

```
Oncology Dataset/
├── ProcessingData_1dataset.R          # R script for GSE25097 processing
├── merge_datasets.py                  # Python script for dataset merging
├── MLAlgorithm.py                     # Machine learning pipeline
├── GSE25097 Series Matrix.txt         # Raw GSE25097 data
├── merged_features_labels.csv         # Preprocessed GSE14520 & GSE121248 data
└── figures/                           # Generated plots and visualizations
```

## 🚀 Usage Instructions

### Step 1: Data Preprocessing (GSE25097)

**⚠️ Important: Update file paths in the script before running**

1. Open `ProcessingData_1dataset.R`
2. Update the input file path on line 11:
   ```r
   input_file <- "/YOUR/PATH/TO/GSE25097 Series Matrix.txt"
   ```
3. Run the R script:
   ```bash
   Rscript ProcessingData_1dataset.R
   ```

**What this script does:**
- Reads GSE25097 Series Matrix file
- Maps probe IDs to gene symbols using GPL10687 platform
- Transposes data (samples as rows, genes as columns)
- Creates binary labels: `tumor=1`, `others (healthy/cirrhotic/non_tumor)=0`
- Splits data into training (80%) and validation (20%) sets
- Outputs: `GSE25097.csv`, `GSE25097_train.csv`, `GSE25097_valid.csv`

### Step 2: Dataset Merging

**⚠️ Important: Ensure input files are in the correct directory**

1. Place the following files in your working directory:
   - `merged_features_labels.csv` (preprocessed GSE14520 & GSE121248)
   - `GSE25097.csv` (output from Step 1)

2. Run the merging script:
   ```bash
   python merge_datasets.py
   ```

**What this script does:**
- Merges two datasets keeping only intersection features (9,923 common genes)
- Adjusts label encoding:
  - `Label 0`: Normal/Non-tumor samples from both datasets
  - `Label 1`: HBV-HCC samples from merged_features_labels.csv
  - `Label 2`: Tumor samples from GSE25097.csv (originally labeled as 1)
- Performs stratified 80:20 train/validation split
- Outputs: `combined_dataset.csv`, `combined_dataset_train.csv`, `combined_dataset_valid.csv`

### Step 3: Machine Learning Analysis

**⚠️ Important: Update file paths in MLAlgorithm.py**

1. Open `MLAlgorithm.py`
2. Update the dataset paths on lines 28-29:
   ```python
   train_df = pd.read_csv("/YOUR/PATH/TO/combined_dataset_train.csv")
   test_df = pd.read_csv("/YOUR/PATH/TO/combined_dataset_valid.csv")
   ```

3. Run the machine learning pipeline:
   ```bash
   python MLAlgorithm.py
   ```

**What this script does:**
- Trains multiple classifiers: Random Forest, Logistic Regression, SVM, Decision Tree
- Evaluates performance using accuracy, F1-score, and multi-class ROC-AUC
- Generates visualizations:
  - Confusion matrices
  - Multi-class ROC curves
  - Feature importance plots
- Identifies top predictive genes for each classifier

## 📈 Data Processing Pipeline

### 1. Raw Data Processing
- **GSE25097**: Processed from Series Matrix format
  - Probe-to-gene mapping using GPL10687 platform
  - Sample classification based on tissue type
  - Quality control and data normalization

### 2. Multi-Dataset Integration
- **Feature Alignment**: Only genes present in both datasets are retained
- **Label Harmonization**: Consistent labeling scheme across datasets
- **Batch Effect Consideration**: Datasets processed separately then merged

### 3. Label Encoding Strategy
```
Original Labels → Final Labels
GSE25097: tumor → 2 (HCC)
GSE25097: healthy/cirrhotic/non_tumor → 0 (Normal/Non-tumor)
GSE14520/GSE121248: HBV-HCC → 1 (HBV-HCC)
GSE14520/GSE121248: Normal → 0 (Normal/Non-tumor)
```



## 📝 Citation

If you use this pipeline in your research, please cite the original GEO datasets:
- GSE14520: [Original paper citation]
- GSE25097: [Original paper citation]  
- GSE121248: [Original paper citation]

## 🤖 Development Tools

This project was developed with assistance from:
- **Cursor**: AI-powered code editor for enhanced development workflow
- **ChatGPT**: AI assistant for code generation, debugging, and documentation

The combination of these AI tools significantly accelerated the development process, enabling rapid prototyping, code optimization, and comprehensive documentation generation.

## 📄 License

This project is provided for research and educational purposes. Please ensure compliance with GEO data usage policies and cite original data sources appropriately.

---

**Note**: This pipeline is designed for research purposes. For clinical applications, additional validation and regulatory compliance may be required.
