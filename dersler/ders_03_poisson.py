r"""Ders 3: Gozlenen, beklenen ve Poisson. De novo analizinin cekirdegi.

Calistirmak icin:
    python dersler\ders_03_poisson.py

Dort asama:
  1. "beklenen" sayi nereden geliyor
  2. Poisson ile "bu kadar cok gormek sansa mi" sorusunu sormak
  3. p degeri hesaplarken tasmayan yontem (sf, 1-cdf degil)
  4. SIFIR gozlemin ne kadar bilgi verdigi: ust sinir hesabi
"""

from scipy.stats import poisson

# Kaplanis 2020, 31.058 gelisimsel bozukluk uclusu.
# "beklenen" degerler projede olculdu; burada onlari kullanip yorumu ogreniyoruz.
N_UCLU = 31_058

VERI = [
    # gen,        rol,                  gozlenen LoF, beklenen LoF
    ("ANKRD11",   "pozitif kontrol",    111,  0.308),
    ("CHD8",      "pozitif kontrol",     19,  0.497),
    ("KIF1A",     "paralog, missense",    0,  0.298),
    ("KIF13B",    "paralog, null",        2,  0.337),
    ("KIF13A",    "HEDEF",                0,  0.367),
]


def p_zenginlesme(gozlenen, beklenen):
    """P(en az bu kadar gormek | sadece sans).

    sf(k-1, lam) = P(X >= k). '1 - cdf' yerine bunu kullan:
    1 - cdf ucta 1.0'a yuvarlanip sahte bir p = 0 uretir.
    """
    if gozlenen == 0:
        return 1.0
    return poisson.sf(gozlenen - 1, beklenen)


def ust_sinir_lambda(gozlenen=0, guven=0.95):
    """Bu kadar az gozlemle uyumlu EN BUYUK gercek beklenti nedir.

    Sifir icin analitik: P(0 | lam) = exp(-lam) = 1 - guven  ->  lam = -ln(1 - guven)
    """
    from math import log
    if gozlenen == 0:
        return -log(1 - guven)
    # genel hal: chi-kare ile
    from scipy.stats import chi2
    return chi2.ppf(guven, 2 * (gozlenen + 1)) / 2


if __name__ == "__main__":
    print("1) 'beklenen' ne demek: her gende sansa bagli da novo mutasyon olur.")
    print("   beklenen = 2 x uclu_sayisi x (o gen icin generasyon basina mutasyon olasiligi)")
    print("   Ters cevirip o olasiligi gorelim:\n")
    for gen, _, _, bek in VERI:
        oran = bek / (2 * N_UCLU)
        print(f"   {gen:9s} beklenen {bek:.3f}  ->  gen basina olasilik {oran:.2e}")

    print(f"\n2) Poisson testi ({N_UCLU:,} uclu):\n")
    print(f"   {'gen':9s} {'rol':20s} {'gozl':>5s} {'bekl':>7s} {'kat':>8s} {'p':>12s}")
    print("   " + "-" * 68)
    for gen, rol, goz, bek in VERI:
        p = p_zenginlesme(goz, bek)
        kat = goz / bek if bek else float("nan")
        print(f"   {gen:9s} {rol:20s} {goz:5d} {bek:7.3f} {kat:8.1f} {p:12.2e}")

    print("\n3) TUZAK: '1 - cdf' ucta tasar ve sahte p = 0 verir.")
    goz, bek = 111, 0.308
    print(f"   ANKRD11 icin  1 - poisson.cdf({goz-1}, {bek}) = {1 - poisson.cdf(goz-1, bek)}")
    print(f"   ANKRD11 icin  poisson.sf({goz-1}, {bek})      = {poisson.sf(goz-1, bek):.2e}")
    print("   Ikincisi dogru. Ilki '0' diyor ve bu bir sayi degil, bir yuvarlama hatasi.")

    print("\n4) SIFIR gozlem ne kadar bilgi veriyor?")
    for etiket, bek in [("Kaplanis tek basina", 0.367), ("butun kataloglar birlesik", 1.16)]:
        lam_max = ust_sinir_lambda(0)
        print(f"   {etiket:28s} beklenen {bek:.3f}  ->  "
              f"gercek beklenti en fazla {lam_max:.2f}  ->  zenginlesme en fazla {lam_max/bek:.1f} kat")
    print("\n   Yani sifir gormek 'etki yok' demiyor, 'etki su kadardan buyuk olamaz' diyor.")
