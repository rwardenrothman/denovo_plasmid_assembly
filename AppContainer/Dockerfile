FROM 409124030547.dkr.ecr.us-east-1.amazonaws.com/foundry-plasmidseq-base:2023.09.29.13.39

# Install the function's dependencies
ARG SPACE_UN
ARG SPACE_TOKEN
RUN pip3 install geneiousdb -i https://$SPACE_UN:$SPACE_TOKEN@pypi.pkg.jetbrains.space/grobio-foundry/p/wf-dev/grobio-foundry-pypi/simple
RUN pip3 install grobiosimplegui -i https://$SPACE_UN:$SPACE_TOKEN@pypi.pkg.jetbrains.space/grobio-foundry/p/wf-dev/grobio-foundry-pypi/simple
RUN pip3 install aiohttp pandas pydantic requests
RUN pip3 install idtscitoolsapi -i https://$SPACE_UN:$SPACE_TOKEN@pypi.pkg.jetbrains.space/grobio-foundry/p/wf-dev/grobio-foundry-pypi/simple
COPY requirements.txt .
RUN pip3 install -r requirements.txt --target "${LAMBDA_TASK_ROOT}" -U --no-cache-dir


ENV PYTHONPATH "/var/lang/lib:/lib64:/usr/lib64:/var/runtime:/var/runtime/lib:/var/task:/var/task/lib:/opt/lib"

# Copy function code
COPY /app ${LAMBDA_TASK_ROOT}

# Set the CMD to your handler
CMD ["app.handler"]