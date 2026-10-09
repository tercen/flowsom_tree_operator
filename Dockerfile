# FlowSOM 2.4.0 (Bioconductor 3.15), ggplot2, jsonlite and the tercen R client come with the
# runtime; nothing else is installed.
FROM tercen/runtime-flowsuite:3.15-1

COPY main.R operator.json /operator/
COPY R /operator/R
WORKDIR /operator

ENV TERCEN_SERVICE_URI=https://tercen.com
ENV OPENBLAS_NUM_THREADS=1

ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
