#!/usr/bin/env node

/**
 * Build script to copy molstar assets from node_modules to dist/vendor
 */

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const sourceDir = path.join(__dirname, '../node_modules/molstar/build/viewer');
const targetDir = path.join(__dirname, '../dist/vendor/molstar');

// Create target directory if it doesn't exist
if (!fs.existsSync(targetDir)) {
  fs.mkdirSync(targetDir, { recursive: true });
}

// Files to copy
const filesToCopy = [
  { src: 'molstar.js', dest: 'molstar.js' },
  { src: 'molstar.css', dest: 'molstar.css' }
];

console.log('Building frontend assets...');
console.log(`Source: ${sourceDir}`);
console.log(`Target: ${targetDir}`);

filesToCopy.forEach(({ src, dest }) => {
  const sourcePath = path.join(sourceDir, src);
  const targetPath = path.join(targetDir, dest);
  
  if (fs.existsSync(sourcePath)) {
    fs.copyFileSync(sourcePath, targetPath);
    console.log(`✓ Copied ${src} -> ${dest}`);
  } else {
    console.warn(`⚠ File not found: ${sourcePath}`);
  }
});

console.log('Build complete!');

