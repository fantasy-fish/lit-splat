#!/bin/bash
set -e  # Exit on error

# Build and deploy to Cloud Run
gcloud run deploy gaussian-splat-viewer \
    --source . \
    --platform managed \
    --region us-central1 \
    --allow-unauthenticated \
    --set-env-vars BUCKET_NAME=splat-upload \
    --memory 2Gi \
    --timeout 3600 \
    --min-instances 1 \
    --max-instances 10 \
    --port 8080 \
    --cpu 2 \
    --concurrency 80 \
    --execution-environment gen2

# Wait for deployment to complete
echo "Waiting for deployment to stabilize..."
sleep 10

# Get the service URL
SERVICE_URL=$(gcloud run services describe gaussian-splat-viewer \
    --platform managed \
    --region us-central1 \
    --format 'value(status.url)')

echo "Service deployed at: $SERVICE_URL"

# Update config.js with the service URL
echo "Updating config.js with service URL..."
sed -i.bak "s|YOUR_CLOUD_RUN_URL|$SERVICE_URL|g" static/config.js

# Verify CORS configuration
echo "Verifying CORS configuration..."
if gsutil cors get gs://splat-upload > /dev/null 2>&1; then
    echo "CORS already configured"
else
    echo "Setting up CORS..."
    gsutil cors set cors.json gs://splat-upload
fi

# Print deployment info
echo "Deployment complete!"
echo "Website URL: $SERVICE_URL"