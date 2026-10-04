r"""Ders 1: Bir genin kisitlilik sayilarini gnomAD'den cekmek.

Calistirmak icin:
    python dersler\ders_01_kisitlilik.py

Ogrendigimiz sey: bir veri tabanina soru sormak ve cevabi okumak.
Bu projede yaptigimiz her seyin altinda bu var.
"""

import requests

# gnomAD'nin tek bir adresi var, butun sorular buraya gidiyor.
API = "https://gnomad.broadinstitute.org/api"

# Soruyu GraphQL diliyle yaziyoruz.
# Duz ceviri: "Bana bu gen sembolunun, GRCh38 surumunde,
#              MANE transkriptinin kisitlilik sayilarini ver."
SORU = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    gene_id
    symbol
    mane_select_transcript { ensembl_id }
    gnomad_constraint {
      pli
      oe_lof
      oe_lof_upper
      mis_z
      oe_mis
      obs_mis
      exp_mis
      obs_lof
      exp_lof
    }
  }
}
"""


def kisitlilik(sembol):
    """Bir gen sembolu al, gnomAD'ye sor, cevabi sozluk olarak dondur."""
    cevap = requests.post(API, json={"query": SORU % sembol}, timeout=60)
    cevap.raise_for_status()          # HTTP hatasi varsa burada patlasin
    veri = cevap.json()

    # GraphQL hatalari HTTP 200 ile de gelebilir, o yuzden ayrica bakiyoruz.
    if "errors" in veri:
        raise RuntimeError(veri["errors"])

    gen = veri["data"]["gene"]
    k = gen["gnomad_constraint"]
    return {
        "sembol": gen["symbol"],
        "transkript": gen["mane_select_transcript"]["ensembl_id"],
        "pLI": k["pli"],
        "LOEUF": k["oe_lof_upper"],
        "oe_lof": k["oe_lof"],
        "missense_Z": k["mis_z"],
        "oe_mis": k["oe_mis"],
        "gozlenen_mis": k["obs_mis"],
        "beklenen_mis": k["exp_mis"],
    }


if __name__ == "__main__":
    for sembol in ["KIF13A", "KIF1A", "TCF4"]:
        s = kisitlilik(sembol)
        print(f"\n{s['sembol']}  ({s['transkript']})")
        print(f"  pLI          {s['pLI']:.4f}    tek kopya kaybina dayaniksiz mi (1'e yakin = evet)")
        print(f"  LOEUF        {s['LOEUF']:.3f}     dusuk = daha kisitli")
        print(f"  oe_lof       {s['oe_lof']:.3f}     beklenen gen susturucu varyantlarin kaci gorulmus")
        print(f"  missense Z   {s['missense_Z']:.2f}      yuksek = missense'e kisitli")
        print(f"  oe_mis       {s['oe_mis']:.3f}     {s['gozlenen_mis']:.0f} gozlenen / {s['beklenen_mis']:.1f} beklenen")
