FROM python:3.11-slim

WORKDIR /app

# Minimal system dependencies only
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy files
COPY requirements.txt .
COPY rabies_proteomics_project.py .
COPY rabies_visualization_extensions.py . 

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Create output directories
RUN mkdir -p /app/rabies_output /app/visualizations

ENV PYTHONUNBUFFERED=1

EXPOSE 8501 8050

CMD ["python", "rabies_proteomics_project.py"]
