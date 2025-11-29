#!/bin/bash

# HGM Deployment Script for Google Cloud Platform
# This script automates the deployment process

set -e  # Exit on error

echo "🚀 HGM Google Cloud Deployment Script"
echo "======================================"
echo ""

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    echo "❌ Error: gcloud CLI is not installed."
    echo "Please install it from: https://cloud.google.com/sdk/docs/install"
    exit 1
fi

# Get project ID
echo "📋 Step 1: Project Configuration"
echo "Current gcloud project: $(gcloud config get-value project 2>/dev/null)"
read -p "Enter your Google Cloud Project ID (or press Enter to use current): " PROJECT_ID

if [ -z "$PROJECT_ID" ]; then
    PROJECT_ID=$(gcloud config get-value project 2>/dev/null)
fi

if [ -z "$PROJECT_ID" ]; then
    echo "❌ Error: No project ID specified."
    exit 1
fi

echo "Using project: $PROJECT_ID"
gcloud config set project $PROJECT_ID

# Enable required APIs
echo ""
echo "🔧 Step 2: Enabling Required APIs..."
gcloud services enable run.googleapis.com --quiet
gcloud services enable containerregistry.googleapis.com --quiet
gcloud services enable cloudbuild.googleapis.com --quiet
echo "✅ APIs enabled"

# Choose deployment method
echo ""
echo "📦 Step 3: Choose Deployment Method"
echo "1) Cloud Run (Recommended - Serverless, Auto-scaling)"
echo "2) Compute Engine (VM with full control)"
echo "3) Build Docker image only (for testing)"
read -p "Enter your choice (1-3): " DEPLOY_METHOD

case $DEPLOY_METHOD in
    1)
        echo ""
        echo "🚢 Deploying to Cloud Run..."
        
        # Get region
        read -p "Enter region (default: us-central1): " REGION
        REGION=${REGION:-us-central1}
        
        # Get memory
        read -p "Enter memory allocation (default: 4Gi): " MEMORY
        MEMORY=${MEMORY:-4Gi}
        
        # Get CPU
        read -p "Enter CPU allocation (default: 2): " CPU
        CPU=${CPU:-2}
        
        # Deploy
        echo "Deploying with: Region=$REGION, Memory=$MEMORY, CPU=$CPU"
        gcloud run deploy hgm-app \
            --source . \
            --region $REGION \
            --platform managed \
            --allow-unauthenticated \
            --memory $MEMORY \
            --cpu $CPU \
            --timeout 3600 \
            --max-instances 10 \
            --quiet
        
        echo ""
        echo "✅ Deployment complete!"
        echo "Your application URL:"
        gcloud run services describe hgm-app --region $REGION --format='value(status.url)'
        ;;
        
    2)
        echo ""
        echo "💻 Creating Compute Engine VM..."
        
        # Get zone
        read -p "Enter zone (default: us-central1-a): " ZONE
        ZONE=${ZONE:-us-central1-a}
        
        # Get machine type
        read -p "Enter machine type (default: n1-standard-4): " MACHINE_TYPE
        MACHINE_TYPE=${MACHINE_TYPE:-n1-standard-4}
        
        # Ask about GPU
        read -p "Do you need GPU support? (y/n, default: n): " GPU_SUPPORT
        
        if [ "$GPU_SUPPORT" = "y" ]; then
            gcloud compute instances create hgm-vm \
                --zone=$ZONE \
                --machine-type=$MACHINE_TYPE \
                --accelerator=type=nvidia-tesla-t4,count=1 \
                --image-family=pytorch-latest-gpu \
                --image-project=deeplearning-platform-release \
                --boot-disk-size=100GB \
                --metadata="install-nvidia-driver=True" \
                --quiet
        else
            gcloud compute instances create hgm-vm \
                --zone=$ZONE \
                --machine-type=$MACHINE_TYPE \
                --image-family=ubuntu-2204-lts \
                --image-project=ubuntu-os-cloud \
                --boot-disk-size=100GB \
                --quiet
        fi
        
        # Create firewall rule
        gcloud compute firewall-rules create allow-streamlit \
            --allow tcp:8501 \
            --source-ranges 0.0.0.0/0 \
            --description "Allow Streamlit traffic" \
            --quiet 2>/dev/null || echo "Firewall rule already exists"
        
        echo ""
        echo "✅ VM created successfully!"
        echo "External IP:"
        gcloud compute instances describe hgm-vm --zone=$ZONE --format='get(networkInterfaces[0].accessConfigs[0].natIP)'
        echo ""
        echo "To SSH into your VM:"
        echo "  gcloud compute ssh hgm-vm --zone=$ZONE"
        echo ""
        echo "Then run these commands on the VM:"
        echo "  git clone https://github.com/vasudevgupta31/oam_iitk.git"
        echo "  cd oam_iitk"
        echo "  pip3 install -r requirements.txt"
        echo "  streamlit run app.py --server.port 8501 --server.address 0.0.0.0"
        ;;
        
    3)
        echo ""
        echo "🐳 Building Docker image..."
        docker build -t gcr.io/$PROJECT_ID/hgm-app:latest .
        
        read -p "Do you want to push to Container Registry? (y/n): " PUSH_IMAGE
        if [ "$PUSH_IMAGE" = "y" ]; then
            docker push gcr.io/$PROJECT_ID/hgm-app:latest
            echo "✅ Image pushed to gcr.io/$PROJECT_ID/hgm-app:latest"
        fi
        
        read -p "Do you want to test locally? (y/n): " TEST_LOCAL
        if [ "$TEST_LOCAL" = "y" ]; then
            echo "Starting container on http://localhost:8080"
            docker run -p 8080:8080 gcr.io/$PROJECT_ID/hgm-app:latest
        fi
        ;;
        
    *)
        echo "❌ Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "🎉 Deployment process complete!"
echo ""
echo "📚 Next Steps:"
echo "1. Set up environment variables for email notifications"
echo "2. Configure Cloud Storage for persistent data"
echo "3. Set up monitoring and logging"
echo "4. Review the DEPLOYMENT_GUIDE.md for more details"
echo ""
echo "For detailed instructions, see: DEPLOYMENT_GUIDE.md"
