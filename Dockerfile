FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Build essentials for compiling Python packages
    build-essential \
    gcc \
    g++ \
    gfortran \
    \
    # Scientific computing dependencies
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    \
    # BioPython optional dependencies
    ncbi-blast+ \
    \
    # Multiple sequence alignment
    muscle \
    \
    # Graphics/visualization support
    libpng-dev \
    libfreetype6-dev \
    libz-dev \
    \
    # For py3Dmol and 3D visualization
    libglu1-mesa \
    libgl1-mesa-glx \
    \
    # For network analysis and graph layout
    graphviz \
    libgraphviz-dev \
    \
    # For downloading and data handling
    curl \
    wget \
    \
    # Development utilities
    git \
    vim \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY requirements.txt .
COPY rabies_proteomics_project.py .
COPY rabies_visualization_extensions.py .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Create output directories
RUN mkdir -p /app/rabies_output /app/visualizations

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Expose port for Streamlit and Dash apps
EXPOSE 8501 8050

# Default command
CMD ["python", "rabies_proteomics_project.py"]
