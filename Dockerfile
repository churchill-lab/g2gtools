FROM python:3.12
LABEL author="Matthew Vincent <matt.vincent@jax.org>"

RUN apt-get update && apt-get install -y \
    vim

COPY . /src/g2gtools

WORKDIR /src/g2gtools

RUN pip install -U pip; cd /src/g2gtools; pip install .; cd

WORKDIR /app

CMD ["g2gtools"]
