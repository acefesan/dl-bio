"""Add a content hash to every static/*.css|js reference in the site's HTML
(static/style.css -> static/style.css?v=1a2b3c4d), so a browser holding an old
copy fetches the new one as soon as the file changes. Rerun after editing
anything in static/:  python3 stamp_assets.py
"""
import hashlib, re
from pathlib import Path

site = Path(__file__).resolve().parent
ref = re.compile(r'((?:\.\./)?static/([\w.-]+\.(?:css|js)))(?:\?v=[0-9a-f]+)?')
hashes = {p.name: hashlib.sha1(p.read_bytes()).hexdigest()[:8] for p in (site / "static").iterdir() if p.suffix in (".css", ".js")}
changed = 0
for html in sorted(site.rglob("*.html")):
    s = html.read_text()
    t = ref.sub(lambda m: f"{m.group(1)}?v={hashes[m.group(2)]}" if m.group(2) in hashes else m.group(0), s)
    if t != s:
        html.write_text(t)
        changed += 1
print(f"stamped {changed} pages; {len(hashes)} assets")
