"""Which KEGG pathways contain the adenosine receptors, and what else is in them?

Why KEGG: the scGSA paper (Wang & Thakar 2024) used the 2023 KEGG pathway
collection as its real-data gene sets, so this is the collection whose
behaviour the benchmark speaks to.

Source: KEGG REST (https://rest.kegg.jp), human (hsa). Pulled on the run date
below, so membership is whatever KEGG currently lists (can differ slightly from
the 2023 snapshot the paper used).
  /list/pathway/hsa   pathway ids and names
  /link/hsa/pathway   pathway -> gene (NCBI gene ids)
  /list/hsa           NCBI gene id -> symbol
Output: results/kegg/kegg_human.json
"""
import json, re, time, urllib.request, datetime
from pathlib import Path

def get(path):
    for i in range(4):
        try:
            return urllib.request.urlopen("https://rest.kegg.jp/" + path, timeout=60).read().decode()
        except Exception as e:
            time.sleep(2 * (i + 1))
    raise RuntimeError(path)

names = {l.split("\t")[0].replace("path:", ""): l.split("\t")[1].replace(" - Homo sapiens (human)", "")
         for l in get("list/pathway/hsa").strip().split("\n")}
sym = {}
for l in get("list/hsa").strip().split("\n"):
    gid, rest = l.split("\t")[0], l.split("\t")[3] if len(l.split("\t")) > 3 else l.split("\t")[1]
    sym[gid] = rest.split(";")[0].split(",")[0].strip()
members = {}
for l in get("link/hsa/pathway").strip().split("\n"):
    p, g = l.split("\t"); members.setdefault(p.replace("path:", ""), set()).add(sym.get(g, g))
out = {"fetched": datetime.date.today().isoformat(),
       "pathways": {p: {"name": names.get(p, p), "genes": sorted(g)} for p, g in members.items()}}
Path("results/kegg/kegg_human.json").write_text(json.dumps(out))
adora = ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"]
print(len(out["pathways"]), "pathways with genes;", "ADORA membership:")
for a in adora:
    print(a, [p for p, d in out["pathways"].items() if a in d["genes"]])
