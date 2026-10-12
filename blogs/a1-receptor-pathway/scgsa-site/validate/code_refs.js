// Checks that every "ours" code reference in static/pipeline.js still points
// at a line that exists in static/scgsa.js (the side panel finds snippets by
// these strings). Run after editing either file: node validate/code_refs.js
const fs = require("fs"), path = require("path");
const dir = path.join(__dirname, "..", "static");
const src = fs.readFileSync(path.join(dir, "scgsa.js"), "utf8");
const pj = fs.readFileSync(path.join(dir, "pipeline.js"), "utf8");
const finds = [...pj.matchAll(/find: "((?:[^"\\]|\\.)*)"/g)].map(m => m[1].replace(/\\"/g, '"'));
const missing = finds.filter(f => !src.includes(f));
finds.forEach(f => console.log(src.includes(f) ? "ok     " : "MISSING", f));
console.log(missing.length ? `FAIL: ${missing.length} reference(s) not found` : `PASS: all ${finds.length} references found`);
process.exit(missing.length ? 1 : 0);
