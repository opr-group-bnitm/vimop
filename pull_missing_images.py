#!/usr/bin/env python3
"""Download ViMOP images through mirror.gcr.io and import them into Docker's containerd.

Docker Hub pulls are extremely slow on this host; aria2 verifies every registry blob
by its manifest digest before an OCI import.
"""
import hashlib, json, os, pathlib, shutil, subprocess, sys, tarfile, urllib.request

ROOT = pathlib.Path('/data/work/vimop/image_cache')
MIRROR = 'https://mirror.gcr.io'
IMAGES = [
    ('oprgroup/medaka', '1.0.3'),
    ('oprgroup/report', '1.0.0'),
    ('oprgroup/structural_variants', '1.0.0'),
]

def get_json(url):
    req=urllib.request.Request(url,headers={'Accept':'application/vnd.docker.distribution.manifest.v2+json','User-Agent':'docker/29'})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

def sha256_file(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024), b''): h.update(b)
    return h.hexdigest()

def run(cmd):
    print('+', ' '.join(map(str,cmd)), flush=True)
    subprocess.run(cmd, check=True)

for repo, tag in IMAGES:
    name = repo.split('/')[-1]
    out = ROOT/name/tag
    out.mkdir(parents=True, exist_ok=True)
    manifest = get_json(f'{MIRROR}/v2/{repo}/manifests/{tag}')
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2))
    descs = [manifest['config']] + manifest['layers']
    for d in descs:
        digest=d['digest'].split(':',1)[1]
        blob=out/'blobs'/digest
        blob.parent.mkdir(parents=True, exist_ok=True)
        if blob.exists() and blob.stat().st_size == d['size'] and sha256_file(blob)==digest:
            print('verified existing', digest, flush=True); continue
        if blob.exists(): blob.unlink()
        url=f'{MIRROR}/v2/{repo}/blobs/{d["digest"]}'
        run(['aria2c','-x16','-s16','-k1M','--file-allocation=none',
             f'--checksum=sha-256={digest}','--max-tries=20','--retry-wait=10',
             '--timeout=120','--connect-timeout=30','--summary-interval=30',
             '-d',str(blob.parent),'-o',digest,url])
    # Make an OCI image layout with original Docker manifest/config/layers.
    layout=out/'oci'; shutil.rmtree(layout, ignore_errors=True)
    (layout/'blobs'/'sha256').mkdir(parents=True)
    for d in descs:
        digest=d['digest'].split(':',1)[1]
        shutil.copy2(out/'blobs'/digest, layout/'blobs'/'sha256'/digest)
    oci_manifest=dict(manifest)
    oci_manifest['mediaType']='application/vnd.oci.image.manifest.v1+json'
    mb=json.dumps(oci_manifest,separators=(',',':')).encode()
    md=hashlib.sha256(mb).hexdigest()
    (layout/'blobs'/'sha256'/md).write_bytes(mb)
    (layout/'oci-layout').write_text('{"imageLayoutVersion":"1.0.0"}')
    (layout/'index.json').write_text(json.dumps({'schemaVersion':2,'manifests':[{
        'mediaType':'application/vnd.oci.image.manifest.v1+json','digest':'sha256:'+md,
        'size':len(mb),'annotations':{'org.opencontainers.image.ref.name':f'{name}:{tag}'}}]}))
    tar=out/f'{name}-{tag}.oci.tar'
    if not tar.exists():
        run(['tar','-cf',str(tar),'-C',str(layout),'.'])
    # ctr needs root on this host; docker (containerd store, docker-group access)
    # loads the OCI tar fine. This daemon stores the loaded image under the
    # bare OCI ref name ("report:1.0.0") and cannot resolve it by that name,
    # so always tag by numeric image ID from `docker images --format`.
    def ref_map():
        r = subprocess.run(['docker', 'images', '--format', '{{.Repository}}:{{.Tag}} {{.ID}}'],
                           capture_output=True, text=True, check=True)
        m = {}
        for line in r.stdout.splitlines():
            ref, _, iid = line.rpartition(' ')
            m[ref] = iid
        return m
    rm = ref_map()
    if f'{repo}:{tag}' in rm:
        print('already tagged', f'{repo}:{tag}', flush=True)
    elif f'{name}:{tag}' in rm:
        run(['docker', 'tag', rm[f'{name}:{tag}'], f'{repo}:{tag}'])
    else:
        run(['docker', 'load', '-q', '-i', str(tar)])
        rm = ref_map()
        if f'{name}:{tag}' not in rm:
            raise SystemExit(f'load did not produce {name}:{tag}; refs={list(rm)}')
        run(['docker', 'tag', rm[f'{name}:{tag}'], f'{repo}:{tag}'])
    print('IMPORTED', f'{repo}:{tag}', flush=True)
