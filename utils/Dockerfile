FROM 'google/cloud-sdk:slim'

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    man-db \
    pkg-config \
    python3-venv \
    software-properties-common \
    unzip \
    wget \
    tabix \
    && \
    # clean up apt cache
    rm -rf /var/lib/apt/lists/*

# Install Java 8 for hail
RUN mkdir -p /etc/apt/keyrings
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | tee /etc/apt/keyrings/adoptium.asc && echo "deb [signed-by=/etc/apt/keyrings/adoptium.asc] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list && apt-get update && apt install -y --no-install-recommends temurin-8-jdk

# Install python packages
RUN apt-get update \
    && apt-get dist-upgrade -y \
    && apt-get install -y --no-install-recommends\
    python3 python3-venv python3-pip libpq-dev

# Create a virtual environment and activate it
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip to latest version
RUN python3 -m pip install --upgrade pip

COPY generate_output_vcf_pk.py generate_output_vcf_pk.py

# Install hail and other python libraries
ENV HAIL_VERSION="0.2.122"
RUN python3 -m pip install \
    hail==${HAIL_VERSION} \
    gnomad

RUN export SPARK_HOME=$(find_spark_home.py) && \
    curl https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar \
    >$SPARK_HOME/jars/gcs-connector-hadoop2-2.0.1.jar && \
    mkdir -p $SPARK_HOME/conf && \
    touch $SPARK_HOME/conf/spark-defaults.conf && \
    sed -i $SPARK_HOME/conf/spark-defaults.conf \
    -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.enable.*:spark.hadoop.google.cloud.auth.service.account.enable true:' \
    -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile.*:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile /gsa-key/key.json:'

ENV PYSPARK_SUBMIT_ARGS="--driver-memory 8g --executor-memory 8g pyspark-shell"
