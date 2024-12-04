FROM nginx:alpine

# Install Python and dependencies
RUN apk add --no-cache python3 py3-pip python3-dev build-base

# Create and activate virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Set working directory
WORKDIR /app

# Install Python dependencies in virtual environment
COPY requirements.txt .
RUN . /opt/venv/bin/activate && pip install -r requirements.txt

# Copy application code
COPY main.py .

# Copy static files to Nginx directory
COPY static/ /usr/share/nginx/html/

# Copy Nginx configuration
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Install supervisor
RUN apk add --no-cache supervisor

# Copy supervisor configuration
COPY supervisord.conf /etc/supervisord.conf

EXPOSE 8080

# Start services using supervisor
CMD ["/usr/bin/supervisord", "-c", "/etc/supervisord.conf"]