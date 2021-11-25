FROM python:3

WORKDIR /usr/src/sbi

COPY requirements.sbi ./

RUN pip3 install numpy>=1.17
RUN pip3 install torch==1.8.2+cpu torchvision==0.9.2+cpu torchaudio==0.8.2 -f https://download.pytorch.org/whl/lts/1.8/torch_lts.html
RUN pip3 install --no-cache-dir -r requirements.sbi
