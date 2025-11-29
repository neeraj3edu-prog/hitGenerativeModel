# **HGM - Hierarchical Generative Model**

![IIT Kanpur Logo](resources/iitklogo.png)

## **Overview**

HGM is a Hierarchical Generative Model designed for molecular generation and optimization using SMILES sequences. This project integrates deep learning techniques for data processing, training, beam search, and sampling, enabling efficient molecular design for drug discovery and materials science applications.

The repository includes:
- **Command-line pipeline** for batch processing
- **Streamlit web application** for interactive experiment management
- **Boltz Protein-Ligand Binding Analysis** for binding affinity prediction
- **Retro Synthesis Generator** for synthesis pathway generation
- **Bioactivity Prediction Pipeline** for property prediction

---

## **Quick Start**

### **Prerequisites**
- Python 3.10 or newer
- CUDA compatible GPU (optional but recommended for faster training)
- Google Cloud SDK (for deployment)

### **Local Setup**

1. **Clone the repository**:
   ```bash
   git clone https://github.com/neeraj3edu-prog/hitGenerativeModel.git
   cd hitGenerativeModel
   ```

2. **Create and activate environment**:
   ```bash
   conda env create -f environment.yml
   conda activate oam_hgm
   ```
   
   Or using pip:
   ```bash
   pip install -r requirements.txt
   ```

3. **Download pretrained model** (optional):
   - Download from [this link](https://drive.google.com/file/d/1hyMgwQnU9V7u5cKER9dSS_0pwJneWluj/view?usp=drive_link)
   - Place in `/pretrained` directory

4. **Run the application**:
   ```bash
   streamlit run app.py
   ```
   Access at `http://localhost:8501`

---

## **Project Structure**

```
oam_iitk/
├── configs/              # Configuration files for model and pipeline settings
├── input_data/           # Raw input files for training
├── output_data/          # Generated molecules per experiment
├── pretrained/           # Pretrained model files (.h5)
├── funcs/                # Utility functions for data processing
├── processes/            # Process abstractions for pipeline routines
├── memory/               # Generated interim files, models per experiment
├── resources/            # UI resources and static files
├── app.py                # Streamlit application entry point
├── app_link.py           # Connects Streamlit app with backend pipeline
├── main.py               # CLI entry point for pipeline execution
├── bpp/                  # Bioactivity prediction pipeline module
│   ├── models/           # Trained ML models for property prediction
│   ├── dataset_config.py # Configuration for prediction datasets
│   ├── feature_engineering.py # Molecular feature extraction
│   ├── predict.py        # Property prediction functions
│   ├── train_models.py   # ML model training scripts
│   ├── training_data/    # Training datasets
│   └── transformers.py   # Data transformation utilities
├── deploy_now.sh         # Quick deployment script
└── README.md             # This file
```

---

## **Features**

### **1. Molecular Generation**
- **Hierarchical Generative Model** using LSTM/GRU layers
- **Beam Search** for targeted molecule generation
- **Temperature-controlled Sampling** for diversity
- **SMILES validation** and canonicalization

### **2. Bioactivity Prediction**
- Predict EC50 values for GLP-1R, GCGR, and GIP receptors
- Binding affinity (Kd) prediction
- Ensemble ML models with Morgan fingerprints
- Target transformation for skewed distributions

### **3. Boltz Protein-Ligand Binding Analysis**
- Analyze protein-ligand binding affinity
- Confidence scores (PTM, iPTM, pLDDT)
- Affinity scores (IC50, Kd, ΔG)
- **Interactive 3D structure visualization** with NGL Viewer
- Multiple representation modes (cartoon, ball & stick, surface)

### **4. Retro Synthesis Generator**
- Generate synthesis pathways for molecules
- Configurable beam size for pathway diversity
- Integration with external retro synthesis API

### **5. Experiment Management**
- Create and configure experiments via web UI
- Track experiment status and progress
- Browse and analyze completed experiments
- Download results in CSV format
- Email notifications for completion

### **6. Interactive Visualizations**
- Training history charts
- Molecular structure visualization (2D)
- 3D protein-ligand structure viewer
- Tanimoto similarity heatmaps
- Property distribution plots

---

## **Usage**

### **Command Line Interface**

1. Modify `config.ini` with your parameters
2. Run the pipeline:
   ```bash
   conda activate oam_hgm
   python main.py
   ```
3. Results available in `/output_data/{experiment_name}`

### **Web Application**

#### **Create Experiment**
1. **Setup Tab**: Upload SMILES files, name experiment
2. **Processing Tab**: Configure train/validation split, sequence lengths
3. **Model Tab**: Select pretrained model, adjust architecture
4. **Beam Search Tab**: Configure beam width
5. **Sampling Tab**: Set temperature and sample count
6. **Communication Tab**: Configure email notifications
7. **Pipeline Execution Tab**: Review and run

#### **Browse Experiments**
- View all experiments and their status
- Analyze training history with interactive charts
- Explore generated molecules with 2D/3D visualization
- **Retro Synthesis**: Generate synthesis pathways
- **Boltz Analysis**: Analyze protein-ligand binding
- Download results and molecular properties

---

## **API Integrations**

### **Retro Synthesis API**
- **Method**: POST
- **Payload**: `{"smiles": "...", "beam_size": 3}`
- Integration with external retro synthesis service

### **Boltz Binding API**
- **Method**: POST
- **Payload**: `{"ligand_smiles": "..."}`
- **Timeout**: 15 minutes (900 seconds)
- Protein-ligand binding affinity prediction

### **3D Structure Files**
- **Format**: CIF files
- **Viewer**: NGL Viewer (embedded)
- Interactive 3D molecular structure visualization

---

## **Backend Architecture**

### **Core Components**

1. **Hierarchical Generative Model (HGM)**
   - Multi-level deep learning architecture for SMILES generation
   - LSTM/GRU layers with attention mechanisms
   - Teacher forcing with categorical cross-entropy loss
   - Sampling and beam search for generation

2. **Data Processing Pipeline**
   - SMILES canonicalization and validation
   - Data augmentation with multiple valid SMILES
   - Tokenization and vocabulary building
   - Train/validation splits with sequence constraints

3. **Bioactivity Prediction Pipeline (BPP)**
   - Ensemble ML models for property prediction
   - Morgan fingerprints and molecular descriptors
   - Target transformation for distribution handling

4. **Experiment Management System**
   - Isolated experiment environments
   - Status tracking through JSON files
   - Configuration inheritance
   - Resumability for interrupted experiments

### **Process Flow**

1. **Configuration Initialization** → Parse inputs, create experiment directory
2. **Data Preprocessing** → Filter, augment, tokenize SMILES data
3. **Model Training** → Initialize architecture, train with callbacks
4. **Molecule Generation** → Beam search and sampling from checkpoints
5. **Post-Processing** → Property prediction, filtering, visualization

### **Performance Optimizations**

- **Multiprocessing** for parallel sampling
- **GPU Acceleration** for training and inference
- **Memory-mapped files** for large datasets
- **Caching** for fingerprints and descriptors

---

## **Configuration Options**

### **Processing**
- `split`: Train/validation split ratio (0.5-1.0)
- `min_len`: Minimum SMILES sequence length
- `max_len`: Maximum SMILES sequence length
- `augmentation`: Data augmentation factor

### **Model**
- `pretrained_model`: Path to pretrained .h5 file
- `epochs`: Number of training epochs
- `lr`: Learning rate
- `neurons`: Neural network layer sizes
- `dropouts`: Dropout rates for regularization
- `batch_size`: Training batch size

### **Beam Search**
- `width`: Beam width for molecule generation
- `from_epoch`: Starting epoch for beam search

### **Sampling**
- `temp`: Temperature for sampling (higher = more diversity)
- `n_sample`: Number of molecules to sample
- `last_n_epochs`: Number of recent epochs to sample from

---

## **Security Features**

1. **Authentication System**
   - Password-protected access
   - Secure credential storage

2. **Input Validation**
   - Sanitized user inputs
   - Error handling for invalid inputs

3. **Resource Limits**
   - Controlled computational resource usage
   - Timeouts for long-running operations

---

## **Email Notifications**

To enable email notifications:

1. Create `.env` file in project root:
   ```
   EMAIL_USER=your_email@gmail.com
   EMAIL_PASSWORD=your_app_password
   ```
2. Enter email in Communication tab
3. Receive notifications on experiment completion

---

## **Deployment**

For detailed deployment instructions, see [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md).

**Quick Deploy to GCP:**
```bash
./deploy_now.sh
```

**Server Information:**
- Deployed on Google Cloud Platform
- See [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) for server setup details

---

## **Troubleshooting**

### **Common Issues**

- **Missing pretrained model**: Download to `/pretrained` directory
- **CUDA errors**: Check GPU drivers and CUDA installation
- **Memory errors**: Reduce batch size or model complexity
- **Invalid SMILES**: Ensure input file contains valid SMILES strings
- **Indentation errors**: Check Python version compatibility (3.10 recommended)

### **Logs & Debugging**

- Detailed logs in `logs/` directory
- Real-time progress in Streamlit app
- Check `memory/{experiment_name}/status.json` for status

---

## **Advanced Usage**

### **Custom Token Vocabulary**
Edit `configs/fixed_params.py` and update `INDICES_TOKEN` and `TOKEN_INDICES` dictionaries.

### **Custom Neural Network Architecture**
Modify `funcs/helpers_training.py` and adjust the `SeqModel` class.

---

## **Contributing**

Contributions are welcome! Please submit a Pull Request.

---

## **License**

This project is licensed under the MIT License.

---

## **Contact**

For queries or support:
- neeraj3edu@gmail.com
- akshaykakar@gmail.com
- guptavasudelhi@gmail.com

---

## **Acknowledgments**

- IIT Kanpur for project support
- Google Cloud Platform for hosting infrastructure
- RDKit and NGL Viewer communities
