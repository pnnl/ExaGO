FROM node:16.13-alpine3.13

WORKDIR /usr/src/app

COPY . .
RUN npm install --legacy-peer-deps

CMD ["npm","start"]
