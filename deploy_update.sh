#!/bin/bash

################################################################################
# Quick Deployment Script for OAM IITK Application Updates
# 
# This script:
# 1. Transfers updated files to the remote server
# 2. Restarts the Streamlit application
#
# Usage: 
#   chmod +x deploy_update.sh
#   ./deploy_update.sh
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SERVER_IP="136.116.43.177"
SERVER_USER="neeraj"
SERVER_PATH="~/oam_iitk"
LOCAL_PATH="/Users/neeraj/Study/IIT_K/Project/GIT/oam_iitk"

# GCP Configuration (if using gcloud)
GCP_SERVER_NAME="rsgpt-server-cpu"
GCP_ZONE="us-central1-a"
GCP_PROJECT="ml-project-477222"

################################################################################
# Helper Functions
################################################################################

print_header() {
    echo -e "\n${GREEN}========================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${GREEN}========================================${NC}\n"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

################################################################################
# Step 1: Check Connection Method
################################################################################

check_connection() {
    print_header "Step 1: Checking Connection Method"
    
    # Check if gcloud is available
    if command -v gcloud &> /dev/null; then
        print_info "gcloud CLI found. Using GCP connection."
        CONNECTION_METHOD="gcloud"
    elif command -v ssh &> /dev/null; then
        print_info "Using direct SSH connection."
        CONNECTION_METHOD="ssh"
    else
        print_error "Neither gcloud nor ssh found. Please install one of them."
        exit 1
    fi
}

################################################################################
# Step 2: Transfer Files
################################################################################

transfer_files() {
    print_header "Step 2: Transferring Updated Files"
    
    if [ "$CONNECTION_METHOD" = "gcloud" ]; then
        print_info "Transferring files via gcloud scp..."
        
        # Transfer only the app.py file (faster for updates)
        gcloud compute scp "$LOCAL_PATH/app.py" \
            "$GCP_SERVER_NAME:$SERVER_PATH/" \
            --zone "$GCP_ZONE" \
            --project "$GCP_PROJECT"
        
        print_success "Files transferred successfully!"
        
    elif [ "$CONNECTION_METHOD" = "ssh" ]; then
        print_info "Transferring files via rsync..."
        
        # Use rsync for efficient transfer
        rsync -avz --progress \
            --exclude '.git' \
            --exclude '__pycache__' \
            --exclude '*.pyc' \
            --exclude '.DS_Store' \
            --exclude 'memory/*' \
            "$LOCAL_PATH/app.py" \
            "$SERVER_USER@$SERVER_IP:$SERVER_PATH/"
        
        print_success "Files transferred successfully!"
    fi
}

################################################################################
# Step 3: Restart Application
################################################################################

restart_application() {
    print_header "Step 3: Restarting Streamlit Application"
    
    if [ "$CONNECTION_METHOD" = "gcloud" ]; then
        print_info "Connecting to server and restarting application..."
        
        gcloud compute ssh --zone "$GCP_ZONE" "$GCP_SERVER_NAME" --project "$GCP_PROJECT" << 'EOF'
# Stop existing Streamlit process
pkill -f "streamlit run app.py" || true
sleep 2

# Navigate to project directory
cd ~/oam_iitk

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate oam_env

# Start Streamlit in background
nohup streamlit run app.py > streamlit.log 2>&1 &

echo "Streamlit application restarted!"
echo "Check logs with: tail -f ~/oam_iitk/streamlit.log"
EOF
        
    elif [ "$CONNECTION_METHOD" = "ssh" ]; then
        print_info "Connecting to server and restarting application..."
        
        ssh "$SERVER_USER@$SERVER_IP" << 'EOF'
# Stop existing Streamlit process
pkill -f "streamlit run app.py" || true
sleep 2

# Navigate to project directory
cd ~/oam_iitk

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate oam_env

# Start Streamlit in background
nohup streamlit run app.py > streamlit.log 2>&1 &

echo "Streamlit application restarted!"
echo "Check logs with: tail -f ~/oam_iitk/streamlit.log"
EOF
    fi
    
    print_success "Application restarted successfully!"
}

################################################################################
# Step 4: Verify Deployment
################################################################################

verify_deployment() {
    print_header "Step 4: Verifying Deployment"
    
    print_info "Waiting 5 seconds for application to start..."
    sleep 5
    
    print_info "Checking if application is running..."
    
    if [ "$CONNECTION_METHOD" = "gcloud" ]; then
        gcloud compute ssh --zone "$GCP_ZONE" "$GCP_SERVER_NAME" --project "$GCP_PROJECT" \
            --command "ps aux | grep 'streamlit run app.py' | grep -v grep"
    elif [ "$CONNECTION_METHOD" = "ssh" ]; then
        ssh "$SERVER_USER@$SERVER_IP" "ps aux | grep 'streamlit run app.py' | grep -v grep"
    fi
    
    if [ $? -eq 0 ]; then
        print_success "Application is running!"
    else
        print_warning "Application may not be running. Check logs on the server."
    fi
}

################################################################################
# Step 5: Display Access Information
################################################################################

display_info() {
    print_header "Deployment Complete!"
    
    echo -e "${GREEN}✓ Files transferred${NC}"
    echo -e "${GREEN}✓ Application restarted${NC}"
    echo ""
    echo -e "${BLUE}Access your application at:${NC}"
    echo -e "  ${YELLOW}http://$SERVER_IP:8501${NC}"
    echo ""
    echo -e "${BLUE}To view logs:${NC}"
    if [ "$CONNECTION_METHOD" = "gcloud" ]; then
        echo -e "  gcloud compute ssh --zone \"$GCP_ZONE\" \"$GCP_SERVER_NAME\" --project \"$GCP_PROJECT\" -- 'tail -f ~/oam_iitk/streamlit.log'"
    elif [ "$CONNECTION_METHOD" = "ssh" ]; then
        echo -e "  ssh $SERVER_USER@$SERVER_IP 'tail -f ~/oam_iitk/streamlit.log'"
    fi
    echo ""
    echo -e "${BLUE}To stop the application:${NC}"
    if [ "$CONNECTION_METHOD" = "gcloud" ]; then
        echo -e "  gcloud compute ssh --zone \"$GCP_ZONE\" \"$GCP_SERVER_NAME\" --project \"$GCP_PROJECT\" -- 'pkill -f \"streamlit run app.py\"'"
    elif [ "$CONNECTION_METHOD" = "ssh" ]; then
        echo -e "  ssh $SERVER_USER@$SERVER_IP 'pkill -f \"streamlit run app.py\"'"
    fi
    echo ""
}

################################################################################
# Main Execution
################################################################################

main() {
    print_header "OAM IITK Deployment Script"
    
    echo -e "${BLUE}Server IP:${NC} $SERVER_IP"
    echo -e "${BLUE}Local Path:${NC} $LOCAL_PATH"
    echo -e "${BLUE}Remote Path:${NC} $SERVER_PATH"
    echo ""
    
    read -p "Continue with deployment? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Deployment cancelled."
        exit 0
    fi
    
    check_connection
    transfer_files
    restart_application
    verify_deployment
    display_info
    
    print_success "Deployment completed successfully!"
}

# Run main function
main
