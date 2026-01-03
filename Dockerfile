# pyNCA Docker Image
# Multi-stage build for optimized image size

# ---- Builder stage ----
FROM python:3.11-slim as builder

WORKDIR /app

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install dependencies
COPY pyproject.toml .
COPY pynca/ pynca/
COPY README.md .
COPY LICENSE .

# Install package with all dependencies
RUN pip install --no-cache-dir --user .[all]

# ---- Runtime stage ----
FROM python:3.11-slim as runtime

WORKDIR /app

# Copy installed packages from builder
COPY --from=builder /root/.local /root/.local

# Make sure scripts in .local are in PATH
ENV PATH=/root/.local/bin:$PATH

# Copy the package
COPY pynca/ pynca/
COPY pyproject.toml .
COPY README.md .
COPY LICENSE .
COPY examples/ examples/

# Set Python environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Create a non-root user (optional, for security)
# RUN useradd --create-home appuser
# USER appuser

# Default command
CMD ["python", "-c", "import pynca; print(f'pyNCA version: {pynca.__version__}')"]

# Labels
LABEL maintainer="pyNCA Authors"
LABEL description="Python Non-Compartmental Analysis Package"
LABEL version="0.1.0"
