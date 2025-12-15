#!/bin/bash

################################################################################
# Quick Deployment Script for OAM IITK Application
# Uses gcloud to deploy to rsgpt-server-cpu
################################################################################

set -e  # Exit on error

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Configuration
GCLOUD="/Users/neeraj/google-cloud-sdk/bin/gcloud"
SERVER_NAME="rsgpt-server-cpu"
ZONE="us-central1-a"
PROJECT="ml-project-477222"
LOCAL_FILE="/Users/neeraj/Study/IIT_K/Project/GIT/oam_iitk/app.py"
REMOTE_PATH="~/oam_iitk/"

echo -e "${BLUE}🚀 Deploying OAM IITK Application${NC}\n"

# Step 1: Transfer file
echo -e "${BLUE}📦 Step 1: Transferring app.py to server...${NC}"
$GCLOUD compute scp "$LOCAL_FILE" "$SERVER_NAME:$REMOTE_PATH" \
  --zone "$ZONE" \
  --project "$PROJECT"

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ File transferred successfully!${NC}\n"
else
    echo -e "${YELLOW}❌ File transfer failed!${NC}"
    exit 1
fi

# Step 2: Restart application
echo -e "${BLUE}🔄 Step 2: Restarting Streamlit application...${NC}"
$GCLOUD compute ssh --zone "$ZONE" "$SERVER_NAME" --project "$PROJECT" << 'EOF'
# Stop existing process
echo "Stopping existing Streamlit process..."
pkill -f "streamlit run app.py" || true
sleep 2

# Navigate to project directory
cd ~/oam_iitk

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate oam_env

# Start Streamlit in background
echo "Starting Streamlit application..."
nohup streamlit run app.py > streamlit.log 2>&1 &

sleep 3

# Verify it's running
if ps aux | grep -v grep | grep "streamlit run app.py" > /dev/null; then
    echo "✅ Streamlit is running!"
    echo "📊 Process info:"
    ps aux | grep -v grep | grep "streamlit run app.py"
else
    echo "⚠️  Warning: Streamlit process not found. Check logs."
fi
EOF

if [ $? -eq 0 ]; then
    echo -e "\n${GREEN}✅ Application restarted successfully!${NC}\n"
else
    echo -e "\n${YELLOW}⚠️  Application restart may have issues. Check logs.${NC}\n"
fi

# Step 3: Display access information
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}🎉 Deployment Complete!${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "\n${BLUE}🌐 Access your application at:${NC}"
echo -e "   ${YELLOW}http://34.68.97.108:8501${NC}\n"

echo -e "${BLUE}📋 Useful commands:${NC}"
echo -e "   View logs:"
echo -e "   ${YELLOW}$GCLOUD compute ssh --zone \"$ZONE\" \"$SERVER_NAME\" --project \"$PROJECT\" -- 'tail -f ~/oam_iitk/streamlit.log'${NC}\n"

echo -e "   Check status:"
echo -e "   ${YELLOW}$GCLOUD compute ssh --zone \"$ZONE\" \"$SERVER_NAME\" --project \"$PROJECT\" -- 'ps aux | grep streamlit'${NC}\n"

echo -e "${GREEN}✨ New Features Deployed:${NC}"
echo -e "   • Boltz Protein-Ligand Binding Analysis"
echo -e "   • 15-minute timeout for complex analyses"
echo -e "   • Confidence & affinity scores display"
echo -e "   • 3D structure file downloads\n"
