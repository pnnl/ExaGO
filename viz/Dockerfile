FROM node:16.13-alpine3.13

RUN apk add --no-cache python3 make g++

WORKDIR /usr/src/app

COPY . .

RUN npm install --legacy-peer-deps

EXPOSE 8080

CMD ["npm", "start", "--", "--host", "0.0.0.0", "--disable-host-check"]
