#!/usr/bin/env python3
"""Download and verify the official ViMOP database parts with aria2c."""
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(16 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

MANIFEST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/data/work/vimop/db_manifest.yaml")
DEST = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("/data/databases/vimop/v1.0.6")
STAGE = DEST.parent / ".staging-v1.0.6"

text = MANIFEST.read_text()
dbs = {}
for match in re.finditer(r"^  (virus|contaminants|centrifuge):\n(.*?)(?=^  (?:virus|contaminants|centrifuge):|\Z)", text, re.M | re.S):
    name, block = match.group(1), match.group(2)
    directory = re.search(r"^    checksum_directory: (\S+)", block, re.M).group(1)
    zipped = re.search(r"^    checksum_zipped: (\S+)", block, re.M).group(1)
    files = re.findall(r"^    - name: (\S+)\n      url: (\S+)\n      checksum: (\S+)", block, re.M)
    dbs[name] = {"directory": directory, "zipped": zipped, "files": files}
if set(dbs) != {"virus", "contaminants", "centrifuge"}:
    raise SystemExit(f"manifest parse failed: {sorted(dbs)}")

STAGE.mkdir(parents=True, exist_ok=True)
for db, spec in dbs.items():
    if (DEST / db).exists():
        print(f"SKIP {db}: verified destination already exists at {DEST / db}", flush=True)
        continue
    stage = STAGE / db
    stage.mkdir(parents=True, exist_ok=True)
    urlfile = stage / "aria2.urls"
    with urlfile.open("w") as fh:
        for fname, url, _ in spec["files"]:
            fh.write(f"{url}\n  out={fname}\n")
    subprocess.run([
        "/data/miniconda3/bin/aria2c", "--continue=true", "--auto-file-renaming=false",
        "--allow-overwrite=false", "--file-allocation=none", "-x16", "-s16", "-j16",
        "--input-file", str(urlfile), "--dir", str(stage),
    ], check=True)
    for fname, _, expected in spec["files"]:
        path = stage / fname
        got = sha256_file(path)
        if got != expected:
            raise SystemExit(f"{db}/{fname}: checksum mismatch {got} != {expected}")
    merged = stage / "merged.tar.xz"
    with merged.open("wb") as out:
        for fname, _, _ in spec["files"]:
            with (stage / fname).open("rb") as part:
                shutil.copyfileobj(part, out, length=16 * 1024 * 1024)
    got = sha256_file(merged)
    if got != spec["zipped"]:
        raise SystemExit(f"{db}: merged checksum mismatch {got} != {spec['zipped']}")
    extract = stage / "extract"
    extract.mkdir(exist_ok=True)
    subprocess.run(["tar", "-xf", str(merged), "-C", str(extract)], check=True)
    candidates = [p for p in extract.iterdir() if p.is_dir()]
    if len(candidates) != 1:
        raise SystemExit(f"{db}: expected one extracted directory, found {candidates}")
    actual_dir = candidates[0]
    entries = []
    for p in sorted(actual_dir.rglob("*")):
        if p.is_file():
            entries.append(sha256_file(p))
    # `find ... -exec sha256sum | sort | awk '{print $1}'` sorts complete
    # sha256sum lines, therefore effectively sorts the hexadecimal hashes.
    entries.sort()
    # Nextflow's `sha256sum` output has a trailing newline before the final
    # hash; preserve it exactly when reproducing checksumDir().
    got = hashlib.sha256(("\n".join(entries) + "\n").encode()).hexdigest()
    if got != spec["directory"]:
        raise SystemExit(f"{db}: directory checksum mismatch {got} != {spec['directory']}")
    final = DEST / db
    if final.exists():
        print(f"VERIFIED {db}: existing destination retained at {final}", flush=True)
    else:
        shutil.copytree(actual_dir, final)
        print(f"VERIFIED {db}: {len(spec['files'])} parts -> {final}", flush=True)

print(f"ALL VERIFIED: {DEST}")
