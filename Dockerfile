FROM node:lts-slim

RUN apt-get update
RUN apt-get -y install python3 python3-dev
RUN apt-get -y install python3-pip
RUN apt-get -y install python-is-python3
RUN apt-get -y install cdo

WORKDIR /usr/src/app
COPY package.json ./
RUN npm install

COPY requirements.txt ./
RUN pip install -r requirements.txt --break-system-packages

RUN apt-get -y install gfortran 


COPY . .

EXPOSE 4000
CMD ["npm", "start"]
