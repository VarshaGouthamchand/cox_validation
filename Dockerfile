FROM python:3.6

WORKDIR /app

RUN touch input.txt
RUN touch output.txt

COPY requirements.txt /app
COPY validate.py /app

RUN pip install -r requirements.txt

CMD ["python", "validate.py"]
