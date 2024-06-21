FROM python:3.10

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 5000

ENV POSTGRES_PASSWORD=${POSTGRES_PASSWORD} \
    OPENAI_KEY=${OPENAI_KEY}

CMD ["python3", "-m", "flask", "--app", "server.py",  "run", "--host=0.0.0.0"]
