FROM python:3.10
LABEL author="Matthew Vincent <matt.vincent@jax.org>"
LABEL version="1.0.0"

COPY . /src/g2gtools

WORKDIR /src/g2gtools

RUN pip install --no-cache-dir .

WORKDIR /app

CMD ["g2gtools"]



