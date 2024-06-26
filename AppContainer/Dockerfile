FROM 409124030547.dkr.ecr.us-east-1.amazonaws.com/foundry-plasmidseq-base:2024.04.17.12.34

# Install the function's dependencies
COPY requirements.txt .
RUN pip3 install -r requirements.txt --target "${LAMBDA_TASK_ROOT}" -U --no-cache-dir


ENV PYTHONPATH "/var/lang/lib:/lib64:/usr/lib64:/var/runtime:/var/runtime/lib:/var/task:/var/task/lib:/opt/lib"

# Copy function code
COPY /app ${LAMBDA_TASK_ROOT}

# Set the CMD to your handler
CMD ["app.handler"]