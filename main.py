from flask import Flask, request, send_from_directory, jsonify
from google.cloud import storage
import os
import uuid
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__, static_url_path='')
bucket_name = os.environ.get('BUCKET_NAME', 'splat-upload')

try:
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
except Exception as e:
    logger.error(f"Failed to initialize storage client: {e}")

@app.route('/')
def root():
    try:
        return send_from_directory('static', 'index.html')
    except Exception as e:
        logger.error(f"Error serving index.html: {e}")
        return jsonify({"error": "Internal server error"}), 500

@app.route('/healthz')
def health_check():
    return jsonify({"status": "healthy"}), 200

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return {'error': 'No file provided'}, 400
    
    file = request.files['file']
    if not file.filename.endswith(('.ply', '.lsplat')):
        return {'error': 'Invalid file format'}, 400

    unique_filename = f"{uuid.uuid4()}-{file.filename}"
    blob = bucket.blob(f"scenes/{unique_filename}")
    
    blob.upload_from_string(
        file.read(),
        content_type='application/octet-stream'
    )

    return {
        'url': f"https://storage.googleapis.com/{bucket_name}/scenes/{unique_filename}"
    }

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))