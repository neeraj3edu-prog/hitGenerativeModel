# Deployment Guide: HGM Streamlit Application

Complete guide for deploying the HGM Streamlit application to Google Cloud Platform or any cloud provider.

---

## **Overview**

This guide covers deploying the Hierarchical Generative Model (HGM) Streamlit application to a cloud server. The application can be deployed on:
- Google Cloud Platform (GCP)
- AWS EC2
- Azure Virtual Machines
- Any Linux server with Python support

---

## **Prerequisites**

- Cloud provider account (GCP, AWS, Azure, etc.)
- SSH access to your server
- Python 3.10 or newer on the server
- Domain name (optional, for production)

---

## **Server Requirements**

### **Minimum Specifications**
- **CPU**: 2 cores
- **RAM**: 4 GB
- **Storage**: 20 GB
- **OS**: Ubuntu 22.04 LTS or newer

### **Recommended Specifications**
- **CPU**: 4+ cores
- **RAM**: 8+ GB
- **Storage**: 50+ GB
- **GPU**: Optional (for model training)

---

## **Quick Deployment**

For quick updates to an existing deployment:

```bash
# From your local machine
./deploy_now.sh
```

This script will:
1. Transfer `app.py` to the server
2. Restart the Streamlit application
3. Verify the application is running

---

## **Full Deployment from Scratch**

### **Step 1: Set Up Your Server**

#### **Create a Server Instance**

**For GCP:**
```bash
gcloud compute instances create your-server-name \
  --zone=your-zone \
  --machine-type=e2-standard-2 \
  --image-family=ubuntu-2204-lts \
  --image-project=ubuntu-os-cloud \
  --boot-disk-size=50GB
```

**For AWS:**
```bash
# Use AWS Console or CLI to create an EC2 instance
# Choose Ubuntu 22.04 LTS AMI
# Select t3.medium or larger instance type
```

#### **Connect to Your Server**

**For GCP:**
```bash
gcloud compute ssh your-server-name --zone=your-zone
```

**For AWS/Other:**
```bash
ssh -i your-key.pem ubuntu@your-server-ip
```

---

### **Step 2: Transfer Application Files**

#### **Option A: Using SCP**

```bash
# From your local machine
scp -r /path/to/hitGenerativeModel user@server-ip:~/
```

#### **Option B: Using Git**

```bash
# On the server
git clone https://github.com/neeraj3edu-prog/hitGenerativeModel.git
cd hitGenerativeModel
```

---

### **Step 3: Install Dependencies**

On the server:

```bash
# Update system packages
sudo apt update
sudo apt upgrade -y

# Install build tools
sudo apt install -y build-essential python3-dev

# Install Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
~/miniconda3/bin/conda init bash
source ~/.bashrc
```

---

### **Step 4: Set Up Python Environment**

```bash
# Navigate to project directory
cd ~/hitGenerativeModel

# Create conda environment with Python 3.10
conda create -n hgm_env python=3.10 -y
conda activate hgm_env

# Install packages via conda (pre-compiled binaries)
conda install -c conda-forge streamlit pandas numpy matplotlib plotly requests python-dotenv pillow rdkit scikit-learn -y

# Install remaining packages via pip
pip install joblib

# If you have requirements.txt
pip install -r requirements.txt
```

---

### **Step 5: Configure Streamlit**

Create Streamlit configuration:

```bash
mkdir -p ~/.streamlit

cat > ~/.streamlit/config.toml << 'EOF'
[server]
port = 8501
address = "0.0.0.0"
headless = true
enableCORS = false
enableXsrfProtection = false

[browser]
serverAddress = "0.0.0.0"
gatherUsageStats = false
EOF
```

---

### **Step 6: Configure Firewall**

#### **For GCP:**

```bash
# From your local machine
gcloud compute firewall-rules create allow-streamlit \
  --direction=INGRESS \
  --priority=1000 \
  --network=default \
  --action=ALLOW \
  --rules=tcp:8501 \
  --source-ranges=0.0.0.0/0
```

#### **For AWS:**

- Go to EC2 Console → Security Groups
- Add inbound rule: TCP port 8501 from 0.0.0.0/0

#### **For Ubuntu UFW:**

```bash
sudo ufw allow 8501/tcp
sudo ufw enable
```

---

### **Step 7: Set Up Secrets (Optional)**

If your application uses secrets:

```bash
mkdir -p ~/.streamlit
nano ~/.streamlit/secrets.toml
```

Add your secrets in TOML format:

```toml
[passwords]
admin = "your-secure-password"

[api]
key = "your-api-key"
```

---

### **Step 8: Launch the Application**

#### **Option A: Run in Foreground (Testing)**

```bash
cd ~/hitGenerativeModel
conda activate hgm_env
streamlit run app.py
```

Press `Ctrl+C` to stop.

#### **Option B: Run in Background with nohup**

```bash
cd ~/hitGenerativeModel
conda activate hgm_env
nohup streamlit run app.py > streamlit.log 2>&1 &

# Check if running
ps aux | grep streamlit

# View logs
tail -f streamlit.log

# Stop the application
pkill -f "streamlit run app.py"
```

#### **Option C: Run with Screen (Recommended for Development)**

```bash
# Start a new screen session
screen -S streamlit

# Inside screen, run the app
cd ~/hitGenerativeModel
conda activate hgm_env
streamlit run app.py

# Detach from screen: Press Ctrl+A then D

# Reattach to screen later
screen -r streamlit

# List all screen sessions
screen -ls

# Kill the screen session
screen -X -S streamlit quit
```

---

### **Step 9: Access the Application**

Get your server's external IP address and access:

```
http://YOUR_SERVER_IP:8501
```

---

## **Production Deployment with systemd**

For production-ready deployment that automatically starts on server reboot:

### **1. Create systemd Service File**

```bash
sudo nano /etc/systemd/system/streamlit-hgm.service
```

Add the following content (replace `YOUR_USERNAME` with your actual username):

```ini
[Unit]
Description=Streamlit HGM Application
After=network.target

[Service]
Type=simple
User=YOUR_USERNAME
WorkingDirectory=/home/YOUR_USERNAME/hitGenerativeModel
Environment="PATH=/home/YOUR_USERNAME/miniconda3/envs/hgm_env/bin"
ExecStart=/home/YOUR_USERNAME/miniconda3/envs/hgm_env/bin/streamlit run app.py
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

### **2. Enable and Start the Service**

```bash
# Reload systemd
sudo systemctl daemon-reload

# Enable service to start on boot
sudo systemctl enable streamlit-hgm

# Start the service
sudo systemctl start streamlit-hgm

# Check service status
sudo systemctl status streamlit-hgm
```

### **3. Manage the Service**

```bash
# Stop the service
sudo systemctl stop streamlit-hgm

# Restart the service
sudo systemctl restart streamlit-hgm

# View logs
sudo journalctl -u streamlit-hgm -f

# Disable auto-start on boot
sudo systemctl disable streamlit-hgm
```

---

## **Updating the Application**

When you make changes to the application:

### **1. Transfer Updated Files**

```bash
# From your local machine
scp app.py user@server-ip:~/hitGenerativeModel/
```

Or use the quick deploy script:
```bash
./deploy_now.sh
```

### **2. Restart the Application**

If using systemd:
```bash
sudo systemctl restart streamlit-hgm
```

If using screen/nohup:
```bash
pkill -f "streamlit run app.py"
cd ~/hitGenerativeModel
conda activate hgm_env
nohup streamlit run app.py > streamlit.log 2>&1 &
```

---

## **Security Best Practices**

### **1. Restrict Firewall Access**

For production, limit access to specific IP addresses:

```bash
# Update firewall rule to allow only your IP
# Replace YOUR_IP with your actual IP address
```

### **2. Use HTTPS with Nginx Reverse Proxy**

```bash
sudo apt update
sudo apt install nginx certbot python3-certbot-nginx

# Configure Nginx as reverse proxy
sudo nano /etc/nginx/sites-available/streamlit
```

Add:

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

```bash
sudo ln -s /etc/nginx/sites-available/streamlit /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx

# Get SSL certificate
sudo certbot --nginx -d your-domain.com
```

### **3. Enable Authentication**

The application has built-in password authentication. Ensure `.streamlit/secrets.toml` is properly configured with strong passwords.

### **4. Regular Updates**

```bash
# Keep system updated
sudo apt update && sudo apt upgrade -y

# Update Python packages
conda activate hgm_env
conda update --all
```

---

## **Monitoring**

### **Set Up Basic Monitoring**

Create a monitoring script:

```bash
nano ~/monitor_streamlit.sh
```

Add:

```bash
#!/bin/bash
if ! pgrep -f "streamlit run app.py" > /dev/null; then
    echo "Streamlit is not running. Starting..."
    cd ~/hitGenerativeModel
    source ~/miniconda3/bin/activate hgm_env
    nohup streamlit run app.py > streamlit.log 2>&1 &
fi
```

```bash
chmod +x ~/monitor_streamlit.sh

# Add to crontab to check every 5 minutes
crontab -e
```

Add:
```
*/5 * * * * /home/YOUR_USERNAME/monitor_streamlit.sh
```

---

## **Troubleshooting**

### **Check if Streamlit is Running**

```bash
ps aux | grep streamlit
netstat -tuln | grep 8501
```

### **View Application Logs**

If using systemd:
```bash
sudo journalctl -u streamlit-hgm -f
```

If using nohup:
```bash
tail -f ~/hitGenerativeModel/streamlit.log
```

### **Common Issues**

#### **Cannot connect to the application**
- Check if firewall rules are configured
- Verify Streamlit is running: `ps aux | grep streamlit`
- Check if port 8501 is listening: `netstat -tuln | grep 8501`

#### **Missing Python packages**
```bash
conda activate hgm_env
pip install missing-package-name
```

#### **Build/compilation errors**
```bash
sudo apt update
sudo apt install -y build-essential python3-dev
```

#### **Python version compatibility**
Use Python 3.10 for best compatibility:
```bash
conda create -n hgm_env python=3.10 -y
```

---

## **Backup and Recovery**

### **Backup Important Files**

```bash
# Backup configuration and data
tar -czf hgm-backup-$(date +%Y%m%d).tar.gz \
  ~/hitGenerativeModel \
  ~/.streamlit/secrets.toml \
  ~/.streamlit/config.toml
```

### **Restore from Backup**

```bash
tar -xzf hgm-backup-YYYYMMDD.tar.gz -C ~/
```

---

## **Performance Optimization**

### **1. Use GPU for Training**

If your server has a GPU:

```bash
# Install CUDA toolkit
# Install TensorFlow with GPU support
conda install tensorflow-gpu
```

### **2. Optimize Memory Usage**

Edit Streamlit config:

```toml
[server]
maxUploadSize = 200
maxMessageSize = 200
```

### **3. Enable Caching**

The application uses Streamlit's caching. Ensure cache directory has sufficient space:

```bash
df -h ~/.streamlit/cache
```

---

## **Quick Reference Commands**

```bash
# Connect to server
ssh user@server-ip

# Activate environment
conda activate hgm_env

# Start app (screen)
screen -S streamlit
cd ~/hitGenerativeModel
streamlit run app.py

# Detach from screen
Ctrl+A then D

# Reattach to screen
screen -r streamlit

# View logs (systemd)
sudo journalctl -u streamlit-hgm -f

# Restart service
sudo systemctl restart streamlit-hgm

# Kill streamlit process
pkill -f "streamlit run app.py"

# Check if running
ps aux | grep streamlit
```

---

## **Additional Resources**

- [Streamlit Documentation](https://docs.streamlit.io/)
- [Google Cloud Platform Documentation](https://cloud.google.com/docs)
- [AWS EC2 Documentation](https://docs.aws.amazon.com/ec2/)
- [Nginx Documentation](https://nginx.org/en/docs/)

---

## **Support**

For deployment issues or questions:
- Check the main [README.md](README.md) for project documentation
- Review application logs for error messages
- Ensure all prerequisites are met

---

**Last Updated**: November 29, 2025
