#!/bin/bash

# Quick deployment using gcloud
echo "🚀 Deploying updated app.py to server..."

# Step 1: Transfer file
echo "📦 Transferring app.py..."
gcloud compute scp app.py rsgpt-server-cpu:~/oam_iitk/ \
  --zone "us-central1-a" \
  --project "ml-project-477222"

if [ $? -eq 0 ]; then
    echo "✅ File transferred successfully!"
    
    # Step 2: Restart application
    echo "🔄 Restarting application..."
    gcloud compute ssh --zone "us-central1-a" "rsgpt-server-cpu" --project "ml-project-477222" << 'EOF'
# Stop existing process
pkill -f "streamlit run app.py" || true
sleep 2

# Navigate and restart
cd ~/oam_iitk
source ~/miniconda3/etc/profile.d/conda.sh
conda activate oam_env
nohup streamlit run app.py > streamlit.log 2>&1 &

echo "✅ Application restarted!"
echo "📊 Access at: http://136.116.43.177:8501"
EOF

    echo ""
    echo "🎉 Deployment complete!"
    echo "🌐 Access your application at: http://136.116.43.177:8501"
else
    echo "❌ File transfer failed. Please check your gcloud configuration."
    exit 1
fi
