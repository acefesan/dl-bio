#!/bin/bash
# Upload a PDF to reMarkable via SSH over Tailscale
# Usage: ./upload_to_remarkable.sh <file.pdf> ["Display Name"]
#
# The reMarkable must be online on Tailscale first:
#   tailscale status | grep remark
#
# If this is the first time, get the SSH password from the tablet:
#   Settings > General > Help > Copyrights and licenses (bottom)
# Then set up key auth: ssh-copy-id root@100.125.228.82

set -euo pipefail

REMARKABLE_HOST="${REMARKABLE_HOST:-root@100.125.228.82}"
XOCHITL_DIR="/home/root/.local/share/remarkable/xochitl"

file="$1"
name="${2:-$(basename "$file" .pdf)}"
uuid=$(uuidgen | tr '[:upper:]' '[:lower:]')
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

cp "$file" "${tmpdir}/${uuid}.pdf"

cat > "${tmpdir}/${uuid}.metadata" <<EOF
{
    "deleted": false,
    "lastModified": "$(date +%s)000",
    "metadatamodified": false,
    "modified": false,
    "parent": "",
    "pinned": false,
    "synced": false,
    "type": "DocumentType",
    "version": 1,
    "visibleName": "${name}"
}
EOF

cat > "${tmpdir}/${uuid}.content" <<EOF
{
    "extraMetadata": {},
    "fileType": "pdf",
    "fontName": "",
    "lastOpenedPage": 0,
    "lineHeight": -1,
    "margins": 100,
    "pageCount": 1,
    "textScale": 1,
    "transform": {
        "m11": 1, "m12": 1, "m13": 1,
        "m21": 1, "m22": 1, "m23": 1,
        "m31": 1, "m32": 1, "m33": 1
    }
}
EOF

mkdir -p "${tmpdir}/${uuid}.cache" "${tmpdir}/${uuid}.highlights" "${tmpdir}/${uuid}.thumbnails"

echo "Uploading '${name}' to reMarkable..."
scp -r ${tmpdir}/* "${REMARKABLE_HOST}:${XOCHITL_DIR}/"
ssh "${REMARKABLE_HOST}" "systemctl restart xochitl"

echo "Done! Uploaded '${name}' as ${uuid}"
echo "The tablet screen will briefly go blank while xochitl restarts."
