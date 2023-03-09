FROM python:3.6

WORKDIR /app
COPY requirements.txt /app
COPY get_lp.py /app

RUN pip install -r requirements.txt

RUN touch input.txt
RUN touch output.txt

CMD ["python", "get_lp.py"]