import path from 'path';
import { fileURLToPath } from 'url';
import express from 'express';
import compression from 'compression';
import morgan from 'morgan';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const app = express();
const distDir = path.join(__dirname, 'dist');
const port = Number(process.env.PORT) || 5500;

app.disable('x-powered-by');
app.use(compression());
app.use(morgan('dev'));
app.use(express.static(distDir, { extensions: ['html'] }));

app.get('*', (req, res, next) => {
  const filePath = path.join(distDir, 'index.html');
  res.sendFile(filePath, err => {
    if (err) next(err);
  });
});

app.listen(port, () => {
console.log(`Frontend server listening on http://localhost:${port}`);
});

