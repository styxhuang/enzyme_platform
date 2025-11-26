FROM node:20-alpine AS base
WORKDIR /app

# Copy package files and install dependencies (including devDependencies for build)
COPY package*.json ./
RUN npm install

# Copy build script and dist template
COPY scripts ./scripts
COPY dist ./dist

# Run build script to copy molstar from node_modules
RUN npm run build

# Install only production dependencies for runtime
RUN npm ci --production

COPY server.js ./server.js

ENV PORT=5500
EXPOSE 5500
CMD ["node", "server.js"]