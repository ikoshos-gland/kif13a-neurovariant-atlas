"""Single probe (no loop) of the documented batch presigned-URL endpoint for ONE small file.

Context: GET /repo/v1/entity/{id}/file returned 403 although /permissions reports
canDownload=true, restrictionLevel=OPEN and no unmet access requirement. This one call
distinguishes "the token lacks the download scope" (403 here too) from an endpoint-specific
refusal. Result goes to raw/synapse/batch_probe.json. No token output.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import syn_common as sc  # noqa: E402

OUT = "/arf/scratch/avural/kif13a_eqtl/raw/synapse"
SYN = "syn61930118"  # specific_and_coloc.tsv, 59,535 bytes


def main() -> int:
    res = {"syn": SYN}
    try:
        fh = sc.file_handle(SYN)
        res["filehandle_ok"] = bool(fh)
        body = {"requestedFiles": [{"fileHandleId": fh["id"], "associateObjectId": SYN,
                                    "associateObjectType": "FileEntity"}],
                "includePreSignedURLs": True, "includeFileHandles": False}
        r = sc.call("POST", f"{sc.FILE}/fileHandle/batch", body, tries=1)
        item = (r.get("requestedFiles") or [{}])[0]
        res["batch_status"] = 200
        res["failureCode"] = item.get("failureCode")
        res["presigned_url_returned"] = bool(item.get("preSignedURL"))
        if item.get("preSignedURL"):
            import urllib.request
            dest = f"{OUT}/specific_and_coloc.tsv"
            with urllib.request.urlopen(item["preSignedURL"], timeout=120) as u:
                data = u.read()
            open(dest, "wb").write(data)
            import hashlib
            res["downloaded_bytes"] = len(data)
            res["md5_ok"] = hashlib.md5(data).hexdigest() == fh.get("contentMd5")
    except sc.AuthError as e:
        res["batch_status"] = e.status
        res["endpoint"] = e.endpoint
    except urllib.error.HTTPError as e:
        res["batch_status"] = e.code
    json.dump(res, open(f"{OUT}/batch_probe.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
