"""Download and verify the data used by the figures.

usage: python download_data.py [--dest data] [--source zenodo|nextcloud] [--url URL]
       python download_data.py --checksum archive.zip     (maintainers: print the SHA-256 to paste below)

The data are a single ZIP archive. Its SHA-256 is checked after download, and the integrity
of every file inside it (ZIP CRC-32) before extraction. Interrupted downloads resume.
"""
import argparse
import hashlib
import sys
import urllib.request
import zipfile
from pathlib import Path

ARCHIVE = 'mecp2_paper_data.zip'
SHA256 = "0e5c2ecafa4e38e247b5bf907e4d13f3177b19c59ec250e04ee65127d0a6d774" 
NEXTCLOUD_SHARE = 'https://neuro.files.academiccloud.de/s/g29FfykzfFse2fL'
ZENODO_RECORD = None  # TODO: record id once published, e.g. '1234567'
EXPECTED = ['Spontaneous', 'Stimulation_Contrast', 'Stimulation_Azimuth', 'MECP2_VR_Behavior']


def url_for(source: str) -> str:
    if source == 'zenodo':
        if ZENODO_RECORD is None:
            sys.exit('No Zenodo record yet: use --source nextcloud')
        return f'https://zenodo.org/records/{ZENODO_RECORD}/files/{ARCHIVE}?download=1'
    return NEXTCLOUD_SHARE.rstrip('/') + '/download'


def download(url: str, target: Path) -> None:
    part = target.with_suffix(target.suffix + '.part')
    done = part.stat().st_size if part.exists() else 0
    request = urllib.request.Request(url, headers={'Range': f'bytes={done}-'} if done else {})
    with urllib.request.urlopen(request) as response:
        if response.status != 206:  # server ignored the range: start over
            done = 0
        total = done + int(response.headers.get('Content-Length', 0))
        with open(part, 'ab' if done else 'wb') as f:
            while chunk := response.read(1 << 20):
                f.write(chunk)
                done += len(chunk)
                print(f'\r{done / 1e9:.2f} / {total / 1e9:.2f} GB', end='', file=sys.stderr)
    print(file=sys.stderr)
    part.rename(target)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dest', type=Path, default=Path('data'), help='folder to extract into')
    parser.add_argument('--source', choices=['zenodo', 'nextcloud'],
                        default='zenodo' if ZENODO_RECORD else 'nextcloud')
    parser.add_argument('--url', help='download from this URL instead')
    parser.add_argument('--keep-zip', action='store_true', help='keep the archive after extraction')
    parser.add_argument('--checksum', type=Path, metavar='ZIP', help='print the SHA-256 of an archive and exit')
    args = parser.parse_args()

    if args.checksum:
        print(sha256(args.checksum))
        return

    args.dest.mkdir(parents=True, exist_ok=True)
    archive = args.dest / ARCHIVE
    if not archive.exists():
        url = args.url or url_for(args.source)
        print(f'Downloading {url}', file=sys.stderr)
        download(url, archive)

    if SHA256 is None:
        print('Warning: no reference checksum set, skipping SHA-256 verification', file=sys.stderr)
    elif (digest := sha256(archive)) != SHA256:
        archive.rename(archive.with_suffix('.corrupt'))
        sys.exit(f'Checksum mismatch (got {digest}): the download is corrupt, run again to retry')

    with zipfile.ZipFile(archive) as z:
        if (bad := z.testzip()) is not None:
            sys.exit(f'Corrupt file in the archive: {bad}')
        z.extractall(args.dest)
    if not args.keep_zip:
        archive.unlink()

    datadir = next((p.parent for p in sorted(args.dest.rglob(EXPECTED[0])) if p.is_dir()), None)
    missing = [d for d in EXPECTED if datadir is None or not (datadir / d).is_dir()]
    if missing:
        sys.exit(f'Extracted, but missing folders: {", ".join(missing)}')
    print(f'Data ready. Run: python run_all.py {datadir}')


if __name__ == '__main__':
    main()
