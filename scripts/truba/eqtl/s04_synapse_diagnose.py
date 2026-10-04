"""One-shot diagnosis of the Synapse 403 on GET /entity/{id}/file (no download retry).

Reads only permission / restriction metadata (view scope) for the root folder and one file and
writes raw/synapse/diagnosis.json. Prints booleans and requirement ids, never the token.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import syn_common as sc  # noqa: E402

OUT = "/arf/scratch/avural/kif13a_eqtl/raw/synapse"
IDS = ["syn61929919", "syn61930218", "syn61930118", "syn61930384"]


def safe(method, url, body=None):
    try:
        return {"ok": True, "data": sc.call(method, url, body, tries=1)}
    except sc.AuthError as e:
        return {"ok": False, "status": e.status, "endpoint": e.endpoint}
    except urllib.error.HTTPError as e:
        return {"ok": False, "status": e.code, "endpoint": method + " " + url.split("?")[0]}


def main() -> int:
    diag = {}
    diag["userProfile"] = safe("GET", f"{sc.REPO}/userProfile")
    if diag["userProfile"].get("ok"):
        up = diag["userProfile"]["data"]
        diag["userProfile"] = {"ok": True}
        uid = up.get("ownerId")
        cert = safe("GET", f"{sc.REPO}/user/{uid}/certifiedUserPassingRecord")
        diag["certified"] = {"ok": cert.get("ok"), "passed": (cert.get("data") or {}).get("passed"),
                             "status": cert.get("status")}
        ver = safe("GET", f"{sc.REPO}/user/{uid}/bundle?mask=63")
        if ver.get("ok"):
            b = ver["data"]
            diag["user_bundle"] = {"isCertified": b.get("isCertified"), "isVerified": b.get("isVerified"),
                                   "isACTMember": b.get("isACTMember"),
                                   "ORCID_present": bool(b.get("ORCID"))}
        else:
            diag["user_bundle"] = ver
    for i in IDS:
        d = {}
        p = safe("GET", f"{sc.REPO}/entity/{i}/permissions")
        if p.get("ok"):
            keep = ["canView", "canDownload", "canEdit", "canCertifiedUserEdit", "isCertifiedUser",
                    "isCertificationRequired", "canPublicRead", "isEntityOpenData", "ownerPrincipalId",
                    "canUpload", "canAddChild"]
            d["permissions"] = {k: p["data"].get(k) for k in keep}
        else:
            d["permissions"] = p
        r = safe("POST", f"{sc.REPO}/restrictionInformation",
                 {"objectId": i, "restrictableObjectType": "ENTITY"})
        d["restrictionInformation"] = r.get("data") if r.get("ok") else r
        a = safe("GET", f"{sc.REPO}/entity/{i}/accessRequirement")
        if a.get("ok"):
            d["accessRequirements"] = [{"id": x.get("id"), "concreteType": x.get("concreteType"),
                                        "accessType": x.get("accessType"),
                                        "termsOfUse_head": (x.get("termsOfUse") or "")[:300]}
                                       for x in a["data"].get("results", [])]
        else:
            d["accessRequirements"] = a
        diag[i] = d
    json.dump(diag, open(f"{OUT}/diagnosis.json", "w"), indent=1)
    print(json.dumps(diag, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
