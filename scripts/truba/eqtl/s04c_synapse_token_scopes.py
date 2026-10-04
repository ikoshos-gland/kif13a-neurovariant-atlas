"""One diagnostic (no download): which scopes do the user's Synapse access tokens carry?

GET /auth/v1/personalAccessToken lists token RECORDS (name, scopes, state, dates), never the token
strings. Prints only name, scopes, state; plus userName and certification state from the profile.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import syn_common as sc  # noqa: E402

AUTH = "https://repo-prod.prod.sagebase.org/auth/v1"
OUT = "/arf/scratch/avural/kif13a_eqtl/raw/synapse"


def main() -> int:
    res = {}
    try:
        recs, tok = [], None
        for _ in range(10):
            url = f"{AUTH}/personalAccessToken" + (f"?nextPageToken={tok}" if tok else "")
            r = sc.call("GET", url, tries=1)
            recs.extend(r.get("results", []))
            tok = r.get("nextPageToken")
            if not tok:
                break
        res["tokens"] = [{"name": x.get("name"), "scopes": x.get("scopes"), "state": x.get("state")}
                         for x in recs]
    except sc.AuthError as e:
        res["tokens"] = {"status": e.status, "endpoint": e.endpoint}
    except urllib.error.HTTPError as e:
        res["tokens"] = {"status": e.code}
    try:
        up = sc.call("GET", f"{sc.REPO}/userProfile", tries=1)
        uid = up.get("ownerId")
        b = sc.call("GET", f"{sc.REPO}/user/{uid}/bundle?mask=63", tries=1)
        res["profile"] = {"userName": up.get("userName"), "isCertified": b.get("isCertified"),
                          "isVerified": b.get("isVerified")}
    except (sc.AuthError, urllib.error.HTTPError) as e:
        res["profile"] = {"error": getattr(e, "status", getattr(e, "code", None))}
    json.dump(res, open(f"{OUT}/token_scopes.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
