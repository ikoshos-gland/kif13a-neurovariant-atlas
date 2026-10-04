# Nakavt ekseni (C ekseni): iki kopya KIF13A kaybı hastalık yapar mı

**Ölçüm tarihi: 3 Ekim 2026.** Bütün hesaplar TRUBA'da SLURM işleri olarak yapıldı (debug bölümü; iş
6439580 indirme, 6439641 hesap). Script'ler: `scripts/truba/nakavt/`. Küçük tablolar:
`data/interim/nakavt/`. İndirilen her dosyanın URL'si, UTC zamanı, bayt sayısı ve sha256'sı:
`data/interim/nakavt/SOURCES_MANIFEST.jsonl`. Doğrulayamadığım her şey açıkça **DOĞRULANMADI** diye
işaretli. Canlı getirilmemiş hiçbir PMID, DOI ya da erişim numarası bu belgede yok.

**Sorulan soru (hipotez H_C):** KIF13A için şimdiye kadar yalnız tek kopya kaybı (haployetmezlik)
soruldu ve büyük ölçüde elendi (`docs/denovo_ekseni.md`). Akraba gen KIF1A'da ise LoF ile hastalık
**çekiniktir**: DDG2P'de G2P00963 "biallelic_autosomal, absent gene product, loss of function"
(KIF1A-related neuropathy, hereditary sensory, definitive), G2P00959 ise monoallelik "altered gene
product structure" (NESCAV). İkisi de 28 Eylül 2026 DDG2P dosyasında bu belge için yeniden okundu.
H_C: iki kopya KIF13A kaybı ağır bir nörogelişimsel ya da ak madde fenotipi yapar.
Çürütücü gözlem: fenotipsiz yetişkin KIF13A nakavtı. Destekleyici gözlem: gelişimsel bozukluk
kohortlarında çekinik KIF13A yükü ya da yetişkin kohortlarda beklenenden az nakavt.

---

## 0. Hüküm

**BİLGİ VERMİYOR.** Dört kaynakta da (gnomAD ile RGC-ME, UK Biobank üzerinden örtüşür) KIF13A için ölçülen sinyal, toleranslı bir genden ve
bilinen bir çekinik gelişimsel bozukluk geninden beklenenler arasında ayırt edilemiyor. Bunun nedeni
ölçüm hatası değil, güç: aynı tasarım, biallelik LoF ile hastalık yaptığı kesin olan **KIF1A'yı da
ayırt edemiyor** (olabilirlik oranı 1,00 ile 1,03 arası).

En önemli üç gözlem:

1. **RGC-ME'de (983.578 kişi) iki homozigot KIF13A "pLOF" taşıyıcısı var** (Ek Tablo 6: `sum_homAA` 2,
   bileşik heterozigot 0). Bu, KIF13A için yayımlanmış ilk insan nakavt kaydıdır. Ama tablo varyantı
   vermiyor; yalnız kümülatif allel frekansını veriyor (7,73e-5). Bu frekans, gnomAD'deki
   **c.4581+1G>C** (son intronun verici bölgesi) ile 1,1 kat tutuyor; bu varyant NMD'den kaçar ve
   ekzon 38 (105 nt, 35 kodon) atlanırsa çerçeve içi bir iç delesyon üretir. Yani büyük olasılıkla bu iki
   kişi **gerçek bir KIF13A nakavtı değildir** (çıkarım, DOĞRULANMADI). Fenotipleri yayımlanmamış.
2. **Genes & Health'te (44.028 otozigot yetişkin) KIF13A nakavtı yok** (22 heterozigot, 0 homozigot,
   0 bileşik heterozigot). Kohortun kendi verisiyle kalibre edilmiş beklenen 0,17 (gerçek LoF için
   0,10), yani sıfır hiçbir şey söylemiyor (P(0) = 0,84).
3. **Chundru 2024'te (29.745 gelişimsel bozukluk triosu) KIF13A'da hiçbir biallelik genotip yok**;
   yazarların nötr beklentisi LoF/LoF için 0,0063. H_C doğru olsa ve KIF13A bilinen çekinik DD
   genlerinin ortalaması kadar zenginleşse bile beklenen 0,052 olurdu (P(0) = 0,95). Bu kohortun H_C'yi
   %95 güçle görmesi için KIF13A'nın 474 kat zenginleşmesi gerekirdi.

Kohortları birleştiren olabilirlik oranı (H_C benzeri / toleranslı benzeri, aşırı dağılıma dayanıklı
ampirik dağılımdan): RGC homozigotları LoF sayılırsa **0,45**, NMD kaçışı sayılırsa **1,57**. İkisi de
0,3 ile 3 arasında, yani yön bile sağlam değil. **popEVE** bağımsız olarak aynı tabloyu çiziyor:
KIF13A aday listesinin hiçbir üyesi gelişimsel bozukluk düzeyinde "ağır" sayılmıyor (katman A'da 52'de
1), KIF1A'da hastalık yapan **birebir aynı 68 değişim** KIF13A'da sıfır ağır skor alıyor.

Desteklenen en güçlü cümle: *"Bugünkü kamuya açık nakavt ve çekinik DD kaynakları KIF13A için çekinik
LoF hipotezini ne destekliyor ne çürütüyor; RGC-ME'deki iki homozigot taşıyıcı büyük olasılıkla NMD'den
kaçan bir son intron varyantı taşıyor ve fenotipleri bilinmiyor."*

---

## 1. Kaynaklar (hepsi 3 Ekim 2026'da canlı doğrulandı)

Kimlikler Europe PMC `search?query=EXT_ID:<PMID> AND SRC:MED` ile getirildi; PMID, PMCID, DOI, yıl,
dergi ve başlık dört makalede de görev metniyle tuttu.

| Kaynak | PMID / PMCID / DOI | Kapsam | Kullanılan dosya | Durum |
|---|---|---|---|---|
| Kim ... van Heel 2026, *Nat Genet*, Genes & Health | 41896352 / PMC13083238 / 10.1038/s41588-026-02553-7 | 44.028 Britanyalı Pakistanlı ve Bangladeşli ekzom, 2.991 nakavt gen | Ek Tablolar 1-26 (`MOESM4`), Ek Bilgi PDF (`MOESM1`) | açık, indirildi |
| Sun ve ark. 2024, *Nature*, RGC-ME | 38768635 / PMC11254753 / 10.1038/s41586-024-07556-0 | 983.578 ekzom, 4.848 nakavt gen | Ek Tablo 1, 2, 6 ve Ek Bilgi DOCX | açık, indirildi |
| Chundru ... Martin 2024, *Nat Genet* | 39313616 / PMC11525179 / 10.1038/s41588-024-01910-8 | 29.745 akrabasız DD triosu (DDD + GeneDx), 22 köken alt grubu | Ek Tablo 1-6, Ek Veri Dosyası (gen düzeyi, 16.655 satır) | açık, indirildi |
| Orenbuch ... Marks 2025, *Nat Genet*, popEVE | 41286104 / PMC12695638 / 10.1038/s41588-025-02400-1 | 31.058 ağır DD hastası, 442 gen, 123 yeni aday | Ek Tablo 1-11; `data.evemodel.org/popeve/v1.1/mutations/<NP>.csv` | açık, indirildi |
| DDG2P (Gene2Phenotype) | | 28 Eylül 2026 sürümü | `ftp.ebi.ac.uk/.../2026_09_28/DDG2P_2026-09-28.csv.gz` | açık, indirildi |
| gnomAD v4.1 | | popülasyon bazlı AC/AN/homozigot, MANE transkript sorgusu; gen kısıtlaması | GraphQL `gnomad.broadinstitute.org/api`; `gnomad.v4.1.constraint_metrics.tsv` | açık |
| Ensembl REST | | MANE ekzon yapısı (NMD sınırı) | `lookup/id/<ENST>?expand=1;utr=1` | açık |

Ek dosyaların indirme kalıbı:
`https://static-content.springer.com/esm/art%3A10.1038%2F<doi-son-ek>/MediaObjects/<dosya>`.

**Kapalı kalan iki kapı (kural gereği durdum):**

- **RGC-ME tarayıcısı** (`rgc-research.regeneron.com/me`) her ziyaretçiyi lisans sayfasına
  yönlendiriyor ve kullanım için "ACCEPT & CONTINUE" tıklaması istiyor; lisans ayrıca kazıyıcıları
  yasaklıyor. Sözleşme kabulü kullanıcı onayı gerektirdiği için tıklamadım. **KIF13A homozigotlarının
  hangi varyantta olduğu ancak burada görülebilir.** Kullanıcının bakması gereken iki satır:
  `6-17771113-C-G` (c.4581+1G>C) ve `6-17764129-C-T` (c.5399G>A, p.Trp1800Ter); ayrıca homozigot
  sayısı 1 veya fazla olan her KIF13A pLOF satırı.
- **Genes & Health özet veri kovası** (`gs://genesandhealth_publicdatasets/`) anonim istekte
  "Request violates VPC Service Controls" (HTTP 403) döndü. Bireysel veri zaten başvuruyla
  açılan güvenli ortamda.

---

## 2. Yöntem

### 2.1 Gerçek LoF tanımı (projeyle aynı)

`scripts/lof_axis.py` ile birebir aynı kural: güvenilir taşıyıcı (AC 1 veya fazla, AC0 filtresi yok),
LOFTEE HC, LOFTEE filtresi yok, ve 50 nt kuralına göre NMD'ye duyarlı (KIF13A'da CDS 4532 = kodon
1511 öncesi). Sağlama: bu çekimde KIF13A için **158** gerçek LoF varyantı çıktı (proje: 158),
KIF1A 102 (102), KIF13B 283 (283). p.Trp1800Ter (LC, END_TRUNC, son ekzon) otomatik dışarıda.

gnomAD v4.1 frekansları (`S["gnomad_q"]`, kümülatif allel frekansı q):

| Gen | sınıf | varyant | AC | q (birleşik) | q (SAS) | q (RGC kökenine ağırlıklı) |
|---|---|---|---|---|---|---|
| KIF13A | bütün pLoF | 204 | 765 | 4,78e-4 | 4,38e-4 | 4,69e-4 |
| KIF13A | LOFTEE HC (kaçış dahil) | 162 | 440 | 2,77e-4 | 3,72e-4 | 2,80e-4 |
| KIF13A | **gerçek LoF** | 158 | 308 | **1,95e-4** | **2,16e-4** | **1,92e-4** |
| KIF1A | gerçek LoF | 102 | 1.173 | 7,71e-4 | 4,33e-4 | 8,08e-4 |
| KIF13B | gerçek LoF | 283 | 853 | 5,95e-4 | 5,51e-4 | 5,96e-4 |

KIF1A'nın yüksek görünen q'su tek bir varyanttan gelir: c.4257-2A>C (ekzon 41/49, AC 824). RGC'nin
KIF1A için kullandığı transkript (ENST00000649096, MANE değil) bu alleli saymıyor: RGC'de KIF1A q =
1,16e-5, yani 50 kat düşük. Bu yüzden KIF1A kalibrasyonunda gnomAD q'su değil, **her kohortun kendi
sayımı** kullanıldı.

### 2.2 Beklenen sayı formülü

İki kopya kayıplı birey (homozigot ya da bileşik heterozigot) için:

```
beklenen = N * [ (1 - F) * q^2 + F * q ]
```

q gen düzeyindeki kümülatif LoF frekansı, F ortalama otozigotluk. Formül F'ye göre doğrusal olduğu
için ortalama F yeterli. Nadir allellerde ikinci terim baskındır: Genes & Health'te KIF13A için HWE
terimi 0,0027, otozigotluk terimi 0,17.

### 2.3 F nereden geldi: kohort içi kalibrasyon

Tahmin edilmiş bir F yerine **her kohortun kendi nakavt sayılarından** etkin F uyduruldu. Genes &
Health Ek Tablo 16, 18.766 genin her biri için heterozigot, homozigot ve bileşik heterozigot taşıyıcı
sayısını veriyor; buradan q = (het + 2 hom + 2 bileşik) / (2N) ve gözlenen nakavt sayısı hesaplandı.
Toleranslı genler (gnomAD v4.1 LOEUF > 0,8 ve DDG2P'de yok) üzerinde Poisson olabilirlik ile tek
parametre F uyduruldu. RGC-ME'de aynısı Ek Tablo 2 (gen başına LOFTEE HC SNV allel sayısı `n`,
16.710 gen) ve Ek Tablo 6 (nakavt sayısı `sum_KO`) birleştirilerek yapıldı.

| Kohort | uydurma kümesi | gen | gözlenen | beklenen | etkin F |
|---|---|---|---|---|---|
| Genes & Health | toleranslı, q 2e-5 ile 2e-3 | 7.726 | 2.336 | 2.326 | 0,0163 |
| Genes & Health | toleranslı, het 10 ile 50 (KIF13A'ya yakın) | 3.587 | 712 | 711 | **0,0154** |
| RGC-ME | toleranslı, q 2e-5 ile 1e-3 | 6.940 | 27.225 | 27.296 | 0,0212 |
| RGC-ME | toleranslı, n 230 ile 920 (KIF13A'ya yakın) | 2.406 | 12.269 | 12.264 | **0,0183** |

Sağlama: Genes & Health'in uydurulmuş F'si (0,0154) makalenin bağımsız ölçtüğü ortalama ROH payıyla
(%1,8) tutarlı. RGC'nin etkin F'si Avrupa ağırlıklı bir kohort için yüksek görünüyor; çünkü RGC'nin q'su
yalnız SNV sayıyor (frameshift'leri F emiyor) ve nakavtların çoğu akraba evliliği yüksek alt
kohortlardan geliyor (19.509 kişilik Pakistan kohortu, Meksiko şehri, Amish). Bu yüzden RGC'de F bir
fiziksel ölçü değil, bir kalibrasyon sabitidir.

### 2.4 Aşırı dağılım: Poisson yetmez, ampirik dağılım kullanıldı

Gen başına nakavt sayıları Poisson'dan çok daha dağınık. RGC'de KIF13A ile aynı frekans bandındaki
1.116 toleranslı genin ortalaması 4,6 nakavt, ama **%42'sinde hiç nakavt yok**. Saf Poisson bu
genlerin yalnız %1'inde sıfır öngörürdü. Bu yüzden asıl test şu: KIF13A'nın gözlemi, aynı frekans
bandındaki (KIF13A'nın değerinin %70 ile %130'u) (a) toleranslı genlerde ve (b) DDG2P'nin biallelik
LoF gelişimsel bozukluk genlerinde (definitive veya strong, "absent gene product" veya LoF mekanizması,
DD paneli; 1.127 gen) ne sıklıkla görülüyor. İkisinin oranı, H_C benzeri ile H0 benzeri arasındaki
olabilirlik oranıdır (LR > 1 H_C lehine).

---

## 3. Beklenen ve gözlenen

`data/interim/nakavt/KIF13A_expected_observed_20261003.tsv`,
`nakavt_summary_20261003.json` içinde `empirical_overdispersion`.

| Kohort | tanım | q | beklenen (H0) | gözlenen | Poisson P | ampirik: aynı banddaki toleranslı genler | ampirik: aynı banddaki çekinik DD genleri | LR |
|---|---|---|---|---|---|---|---|---|
| Genes & Health, 44.028 | kohortun pLoF tanımı | 2,50e-4 (22 het) | 0,172 | **0** | P(0) 0,84 | %85,7'sinde 0 (1.582 gen) | %97,2'sinde 0 (253 gen) | 1,14 |
| Genes & Health | gerçek LoF (SAS payı 0,58) | 1,45e-4 | 0,099 | 0 | P(0) 0,91 | %92,2 (1.648) | %97,3 (257) | 1,05 |
| RGC-ME, 983.578 | RGC pLOF, AAF < %1 | 2,80e-4 (n 460) | 5,12 | **2** | P(X ≤ 2) 0,115 | tam 2: %12,1; en çok 2: %76,9 (1.116) | tam 2: %4,8; en çok 2: %93,4 (166) | 0,40 |
| RGC-ME | gerçek LoF, iki homozigot NMD kaçışıysa | 1,66e-4 (n eşdeğeri 272) | 3,01 | 0 | P(0) 0,049 | %59,2'sinde 0 (1.436) | %87,9'unda 0 (232) | 1,49 |
| gnomAD v4.1, yalnız homozigot | gerçek LoF | 1,95e-4 | 0,46 | 0 | P(0) 0,63 | | | |
| Chundru 2024, 29.745 DD triosu | LoF/LoF, yazarların nötr beklentisi | | 0,0063 | 0 | P(0) 0,99 | | | |
| Chundru 2024 | LoF/LoF, H_C altında (ARDD çarpanı 8,2) | | 0,052 | 0 | P(0) 0,95 | | | |

**Birleşik LR** (G&H ile RGC örtüşmez; G&H, RGC-ME'nin 34 alt kohort listesinde yok):
RGC homozigotları sayılırsa 1,14 x 0,40 = **0,45**; NMD kaçışı sayılırsa 1,05 x 1,49 = **1,57**.

**Okuma:**

1. Saf Poisson RGC'nin gerçek LoF satırında P = 0,049 veriyor ve bu sahte bir sinyaldir. Aynı
   frekanstaki toleranslı genlerin %59'unda da sıfır nakavt var.
2. RGC'de KIF13A'nın 2 nakavtı, toleranslı genlerin dağılımına (%12 tam 2) çekinik DD genlerininkinden
   (%5) daha iyi uyuyor; ama bu ancak o iki kişi gerçek LoF taşıyorsa H_C aleyhine okunabilir, ve
   büyük olasılıkla taşımıyorlar (bölüm 5).
3. Genes & Health'te sıfır nakavt beklenen bir sonuç; 0,17'lik bir beklentiyle sıfırın olasılığı
   %84. Yetişkin bir kohortta sıfır gözlem ancak beklenen 3'ü geçtiğinde H0'ı %5 düzeyinde
   reddedebilir; KIF13A için hiçbir kohort buna yaklaşmıyor.

---

## 4. Kalibrasyon genleri

### 4.1 KIF1A (biallelik LoF hastalık geni, pozitif kontrol)

| Kohort | KIF1A gözlemi | beklenen (H0) | ampirik LR |
|---|---|---|---|
| Genes & Health | 6 het, 0 nakavt | 0,046 | 1,00 (toleranslıların %97,3'ünde, DD genlerinin %97,5'inde 0) |
| RGC-ME | n 19, Ek Tablo 6'da yok (0 nakavt) | 0,21 | 1,03 (%92,2'ye karşı %95,0) |
| gnomAD v4.1 | gerçek LoF homozigotu 0 | 1,52 | P(0) 0,22; beklentinin çoğu tek bir şüpheli allelden (c.4257-2A>C) |
| Chundru 2024 | Ek Tablo 5'te "Consensus" ARDD geni; biallelik 1 proband (Functional/Functional), beklenen 0,39, p 0,32; LoF/LoF beklenen 0,0 | | |

**Sonuç:** tasarım, biallelik LoF ile hastalık yaptığı kesin olan KIF1A'yı da toleranslı bir genden
ayırt edemiyor. KIF1A'nın LoF allelleri KIF13A'nınkinden bile nadir; nakavt kohortları bu frekansta tek
gen düzeyinde kördür. KIF13A'daki "bilgi vermiyor" hükmü bu yüzden yöntemin bir başarısızlığı değil,
bugünkü kohort boyutunun sınırıdır.

### 4.2 Gen kümesi düzeyinde yöntem çalışıyor

Aynı kalibrasyonla DDG2P biallelik LoF DD genlerinin toplamı:

| Kohort | gen | gözlenen nakavt | beklenen (toleranslı F ile) | oran | p (gözlenen kadar az) |
|---|---|---|---|---|---|
| Genes & Health | 1.004 | 97 | 276 | **0,35** | 1,2e-35 |
| RGC-ME | 964 | 943 | 3.645 | **0,26** | < 1e-300 |

Yani bilinen çekinik DD genleri yetişkin kohortlarda 3-4 kat tükenmiş; sıfıra inmiyorlar (hafif
fenotipler, NMD kaçışı, alternatif transkriptler). Bu tükenme Genes & Health makalesinin kendi
bulgusuyla (Ek Not 8: AR ve biallelik DD genlerinde homozigot pLoF tükenmiş) aynı yönde. Kalibrasyon
çalışıyor; KIF13A tek başına yetersiz güçte.

### 4.3 Toleranslı genler ve paralog

- Genes & Health, het 15-29 bandı: toleranslı genlerin %14,3'ünde en az bir nakavt var (ortalama
  0,18; model 0,17 öngörüyor), DD genlerinde %2,8.
- RGC-ME, n 322-598 bandı: toleranslıların %57,7'sinde, DD genlerinin %26,5'inde en az bir nakavt.
- **KIF13B** (paralog, LOEUF 0,663): Genes & Health'te 33 het, 0 nakavt (beklenen 0,26); RGC-ME'de
  **3 homozigot** (beklenen 3,22, kümülatif homozigot AAF 1,32e-5), Chundru'da 0 (Tümü sınıfı
  beklenen 1,20). KIF13B'nin üç homozigotu da varyant kimliği olmadan geliyor.

### 4.4 Chundru boru hattının sıfır kalibrasyonu

Bütün genlerde biallelik sinonim genotipler: gözlenen 9.870, beklenen 10.214 (oran **0,97**).
ARDD konsensüs genlerinde LoF/LoF: 299'a karşı 37,7 (7,9 kat); DDG2P biallelik LoF ile kesişen 850
konsensüs geninde 287'ye karşı 34,9 (**8,2 kat**). Beklentisi KIF13A'nınkine yakın (0,002-0,019) 296
konsensüs ARDD geninin yalnız **%8,4'ünde** en az bir LoF/LoF proband var.

---

## 5. Bulunan nakavtların tek tek dökümü

### 5.1 RGC-ME, KIF13A: 2 homozigot

| Alan | Değer (Ek Tablo 6) |
|---|---|
| transkript | ENST00000259711 (MANE) |
| `allhetRA` (bütün pLOF'larda heterozigot, AAF filtresiz) | 604 |
| `sum_homAA` (nadir homozigot pLOF taşıyıcısı) | **2** |
| `homAlt_aaf` (homozigot varyantların kümülatif AAF'si) | 7,73e-5 |
| `comp_het` | 0 |
| `inExternalDataset` | False |
| faz | gerekmez (homozigot); bileşik heterozigot yok |
| varyant kimliği, ekzon | **tabloda yok** |
| fenotip | **yayımlanmamış** |

**Hangi varyant olabilir** (`KIF13A_rgc_hom_candidates_20261003.tsv`): gnomAD v4.1 frekansları RGC'nin
köken dağılımına (AFR 55.281, EAS 5.849, EUR 722.477, IAM 85.125, MEA 2.229, SAS 30.752, UNK 81.865)
ağırlıklandırılıp 7,73e-5 ile karşılaştırıldı. Tek varyant olarak 3 kat aralığa yalnız ikisi giriyor:

| varyant | LOFTEE | NMD | ağırlıklı AF | oran | SAS / AMR / NFE AF | gnomAD hom |
|---|---|---|---|---|---|---|
| **c.4581+1G>C** (6-17771113-C-G), son intron vericisi | HC | kaçar | 8,54e-5 | **1,10** | 1,44e-4 / 1,18e-4 / 8,1e-5 | 0 |
| c.5399G>A p.Trp1800Ter (6-17764129-C-T) | LC, END_TRUNC | kaçar | 1,20e-4 | 1,55 | 0 / 0 / 1,49e-4 | 1 (NFE) |

Gerçek LoF adaylarının en sık olanı p.Arg1395Ter'in bile ağırlıklı AF'si 1,9e-5; iki gerçek LoF alleli
toplamı 7,7e-5'e ulaşamaz. Bu tabloyu bozabilecek tek senaryo, gnomAD'de olmayan bir RGC alt
kohortuna (Meksiko şehri, Pakistan kohortu) özgü bir kurucu alleldir; dışlanamaz.

c.4581+1G>C'nin olası sonucu: ekzon 38 CDS 4477-4581, **105 nt = 35 kodon**, yani ekzon atlanırsa
kalıntı 1493-1527 çerçeve içi silinir ve son ekzon korunur; intron tutulursa durdurucu kodon son
birleşmenin akışaşağısında kalır. İki durumda da NMD tetiklenmez. Splicing sonucu deneysel olarak
bilinmiyor (**DOĞRULANMADI**). RGC'nin pKO tanımına LOFTEE HC filtresinin RGC-ME için de uygulanıp
uygulanmadığı Ek Bilgi metninde muğlak; uygulandıysa p.Trp1800Ter (LC) zaten dışarıdadır ve aday
c.4581+1G>C'ye iner.

### 5.2 Genes & Health, KIF13A: 0

22 heterozigot pLoF taşıyıcısı, 0 homozigot, 0 bileşik heterozigot (Ek Tablo 16). Fazlanmış 39.148
ekzomluk alt kümede (Ek Tablo 14) pLoF heterozigot 14. Aynı tabloda yan bilgi: KIF13A'da **6 homozigot
"zararlı missense" (pDM) genotipi** var (CADD > 20, PolyPhen2 > 0,445, SIFT zararlı); bunlar LoF değil
ve varyant kimlikleri kamuya açık tablolarda yok. Ek Tablo 16 KIF13A için "In gnomAD v4 = True" diyor;
gnomAD v4.1'deki tek KIF13A LoF homozigotu p.Trp1800Ter'dir (NMD kaçışı), dolayısıyla bu bayrak
büyük olasılıkla o kayda dayanıyor (**DOĞRULANMADI**, G&H'nin harici küme ölçütü yazılmamış).

### 5.3 gnomAD v4.1

Gerçek LoF homozigotu 0. Tek homozigot p.Trp1800Ter (NFE), son 5 aminoasidi keser. Nötr kalibrasyon:
üç gendeki 2.994 nadir sinonim varyantta gözlenen homozigot 1.009, varsayılan F'lerle beklenen 878;
uydurulan F çarpanı 2,17. Bu çarpanla gerçek LoF için beklenen homozigot 0,46 (P(0) = 0,63);
p.Trp1800Ter için beklenen 0,13, gözlenen 1 (P(1 veya fazla) = 0,12).

### 5.4 Chundru 2024: KIF13A

Ana analiz (29.745 proband, kıtalar arası karışım yok): Tümü 0 (beklenen 0,359, p 1),
LoF/LoF 0 (0,0063), LoF/Functional 0 (0,031), Functional/Functional 0 (0,328). KIF13A "Known
recessive genes" sütununda NA. Güney Asya alt gruplarında LoF/LoF beklentisi SAS3 0, SAS4 0,
SAS5 0,0011.

---

## 6. popEVE bulgusu

`data/interim/nakavt/KIF1A_calibration_popEVE_20261003.tsv`, `nakavt_summary_20261003.json` içinde
`popeve`. Ağır eşik makaleden: popEVE ≤ **-5,056** (yüzde 99,99 olasılıkla zararlı bileşende).

**Gen düzeyi (makalenin Ek Tabloları):**

- KIF13A, 442 genlik listede de, Ek Tablo 7'deki yeni adaylarda da **yok**. Ek Tablo 6: ağır DD
  kohortunda 1 de novo missense, 0 ağır, toplama p 0,203, anlamlı değil. (Ek Tablo 5 KIF13A için iki
  SDD satırı listeliyor: D65E ve I403T; Ek Tablo 6'nın neden 1 saydığı yazılmamış, DOĞRULANMADI.)
- Ek Tablo 2: KIF13A missense'lerinin %1,30'u ağır, %3,87'si orta; en düşük skor -5,43. KIF1A'da
  %9,72 ağır, en düşük -6,66.

**Varyant düzeyi (popEVE v1.1 portal dosyaları, NP_071396.4 ve NP_004312.2):**

| Küme | n | ağır (≤ -5,056) | medyan popEVE |
|---|---|---|---|
| KIF13A katman A (gnomAD'de görülen, `KIF13A_candidates_20260915.tsv`) | 52 | **1** (S278W, -5,18) | -4,21 |
| KIF13A katman B / C / D | 73 / 100 skorlu / 33 | 0 / 0 / 0 | -4,02 / -3,82 / -3,81 |
| Doygunluk tablosu katman A, bütün değişimler | 929 | 119 | -4,51 |
| Doygunluk tablosu katman A, tek baz değişimiyle ulaşılabilir | 234 | 4 | -4,34 |
| **KIF1A'da birebir aynı değişimi patojenik olanlar (SNV)** | 68 | **0** | -4,23 |
| KIF13A motor domaininin bütün değişimleri | 6.593 skorlu | 223 | -3,92 |
| KIF1A çekirdek pozitifleri (motor, ClinVar P/LP, proje kürasyonu) | 111 | **86 (%77,5)** | -5,37 |
| KIF1A motor popülasyon negatifleri (gnomAD AC ≥ 3) | 61 | 2 (%3,3) | -3,79 |

De novo ve bilinen KIF13A varyantları: D65E -3,62, R173Q -4,13, Q237H -3,42, I403T -3,77; hiçbiri ağır
değil. Proje doygunluk tablosunun öne çıkardığı K108T -4,62, C97R -4,52, C97Y -4,61, G120V -4,67,
Y94D -4,69, A90D -4,13, G83R -4,32.

**Okuma:** popEVE aynı kalıntıda aynı değişimi KIF1A'da ağır, KIF13A'da orta-hafif sayıyor. Bu,
`docs/bolgesel_kisitlilik.md`'deki "KIF13A'nın motoru KIF1A'nınkinin dörtte biri kadar derin bir seçilim
izi taşıyor" bulgusunun proteom ölçekli, bağımsız bir modelden gelen tekrarıdır. popEVE organizma
düzeyinde ağırlığı ölçmeye çalıştığı için, katman A'nın "ağır gelişimsel bozukluk" yapacağı beklentisine
karşı bir işarettir; mekanizmayı (baskın-negatif gibi) dışlamaz.

**Kapsam uyarısı:** popEVE'nin filtreli skoru KIF13A'da yalnız **1-939** kalıntılarını kapsıyor (912
konum; motor domainin 347/348 konumu dahil). 940-1805 arası kuyruk "unfiltered" sütununda; K1122E
(-3,53), A1146V (-3,32), S1233L (-3,60) gibi kuyruk varyantları bu yüzden yukarıdaki sayımlara girmedi.
Sağlama: filtreli skorlardan hesaplanan KIF13A ağır payı (225 / 17.328 = %1,2985) Ek Tablo 2 ile
**birebir** tutuyor. KIF1A'da ise portal dosyası %12,2 veriyor, Ek Tablo 2 %9,72; fark açıklanamadı
(sürüm farkı olabilir, DOĞRULANMADI). Orta eşik makalede sayı olarak verilmiyor; KIF13A'nın Ek Tablo 2
oranlarından türetilen değer -4,62, KIF1A'nınkinden -4,86; tutarsız olduğu için yalnız ağır eşik
kullanıldı.

**Döngüsellik uyarısı:** popEVE KIF1A motorunda çekirdek pozitifleri popülasyon negatiflerinden AUROC
0,957 ile ayırıyor, AlphaMissense'in 0,925'inden yüksek. Ama popEVE'nin insan varyasyonu bileşeni
UKBB/gnomAD frekanslarını kullanır ve negatifler gnomAD'de görülmüş varyantlardır; bu sayı projenin
kalibrasyon tablosuna **eklenemez**.

---

## 7. Güç: ne kadar veri gerekirdi

- **Yetişkin nakavt kohortu:** Genes & Health ölçeğinde (F = 0,0154) gerçek LoF için kişi başına
  beklenen yaklaşık 2,2e-6. H0 altında beklenenin 3 olması için yaklaşık 1,3 milyon otozigot yetişkin
  gerekirdi. Ama RGC bandı gösteriyor ki ortalaması 4,6 olan toleranslı genlerin bile %42'sinde sıfır
  var; aşırı dağılım yüzünden bu boyut bile tek gen düzeyinde kesin bir ayrım vermez.
- **Gelişimsel bozukluk triosu:** H_C altında beklenen = R x 0,0063. R = 8,2 (ARDD ortalaması) ile
  0,052. Kohortun %95 güçle bir vakayı görmesi için R'nin 474 olması gerekirdi; Chundru'da HECTD4 bu
  düzeyde (2 gözlenen, 0,0042 beklenen), ortalama bir ARDD geni değil.
- **Pratik yol:** tek gen düzeyinde güç kohort büyütmekten değil, **vaka bulmaktan** gelir: RGC'deki
  iki homozigotun varyantı ve fenotipi (geri çağırma), ve akraba evliliği yüksek DD ekzomlarında
  biallelik KIF13A taraması (GeneMatcher türü). Bunların ikisi de insan hesabı ve başvuru gerektirir.

---

## 8. Aşırı iddia uyarıları

1. **"KIF13A nakavtları sağlıklı yetişkinlerde var."** Yanlış. RGC'nin iki homozigotunun varyantı
   yayımlanmamış, büyük olasılıkla NMD kaçışı; fenotipleri bilinmiyor. RGC-ME ayrıca hastalık
   kohortları da içeriyor (MS, glokom, lupus, sedef gibi).
2. **"Genes & Health'te KIF13A nakavtı olmaması çekinik hastalığı destekliyor."** Yanlış. Beklenen
   0,17; sıfırın olasılığı H0 altında %84.
3. **"Chundru KIF13A'yı çekinik DD geni olarak dışladı."** Yanlış. Beklenti 0,0063; aynı beklentideki
   kanıtlanmış ARDD genlerinin %91,6'sı da sıfır proband veriyor.
4. **"RGC'de gerçek LoF nakavtı yok, P = 0,049."** Yanlış. Bu saf Poisson; ampirik dağılımda aynı
   frekanstaki toleranslı genlerin %59'u da sıfır.
5. **"popEVE KIF13A varyantlarının zararsız olduğunu gösterdi."** Fazla. popEVE ağır DD düzeyinde
   sinyal görmüyor; orta etkili, eksik penetranslı ya da DD dışı (ak madde, yetişkin başlangıçlı)
   fenotipleri dışlamaz, ve kuyruğun yarısını filtreli skorda hiç kapsamıyor.
6. Varsayımlar: gnomAD homozigot hesabındaki popülasyon F'leri DD probandlarından (Chundru Tablo 1)
   alındı ve sinonimlerle 2,17 kat düzeltildi; RGC'nin etkin F'si bir kalibrasyon sabitidir;
   Genes & Health'in Ek Tablo 16 pLoF tanımı (HC mi HC+LC mi, hangi MAF eşiği) tablo açıklamasında
   yazmıyor (DOĞRULANMADI); Chundru'nun kullandığı VEP 94 kanonik transkriptinin KIF13A için MANE ile
   aynı olduğu doğrulanmadı.

---

## 9. Bu eksende öğrenilen tuzaklar

1. **Nakavt sayıları Poisson değildir.** Gen başına sayılar kurucu alleller ve akraba evliliği alt
   kohortları yüzünden aşırı dağınık; aynı frekans bandındaki genlerin gerçek dağılımıyla karşılaştır.
2. **RGC Ek Tablo 6 varyant vermez, kümülatif homozigot AAF verir.** Kimlik ancak lisanslı tarayıcıda.
3. **RGC Ek Tablo 2, KIF1A için MANE dışı ENST00000649096 kullanır**; gnomAD MANE'deki AC 824'lük
   c.4257-2A>C bu transkriptte sayılmaz (q 50 kat farklı). Bir genin LoF frekansını kohortlar arası
   taşırken transkripti kontrol et.
4. **popEVE'nin filtreli skoru boşluk oranı yüksek konumları boş bırakır**; KIF13A'da 940-1805 kuyruğu.
   Yüzde hesabında payda boş olmayan skorlardır (Ek Tablo 2 ile böyle tutuyor).
5. **Chundru Ek Veri Dosyası'nda sinonim bloğunun "Expected" başlığı "Exepected" yazılmış**; başlıkla
   sütun arayan kod sessizce boş döner.
6. **Chundru Ek Veri Dosyası'nın zip'i bir `__MACOSX/._*.xlsx` kaynak çatalı taşır**; Excel okuyucu onu
   açmaya çalışırsa çöker.
7. **Genes & Health kamu kovası anonim erişime kapalı** (VPC Service Controls, HTTP 403), makalenin
   "publicly available" ifadesine rağmen.

---

## 10. Yeniden üretim

Hepsi TRUBA'da, `/arf/scratch/avural/kif13a_nakavt`, debug bölümü (internet yalnız orada):

```bash
# yerelden script'leri gönder
scp scripts/truba/nakavt/{fetch_sources.py,inspect_sources.py,nakavt_analysis.py,01_fetch_inspect.slurm,02_analysis.slurm} \
    truba:/arf/scratch/avural/kif13a_nakavt/scripts/
# proje tablolarını gönder (popEVE karşılaştırması için)
scp data/processed/KIF13A_candidates_20260915.tsv data/processed/KIF13A_saturation_lookup_20260916.tsv \
    data/interim/clinvar_curation/KIF1A_calibration_sets_20260915.tsv \
    truba:/arf/scratch/avural/kif13a_nakavt/interim/project/
# 1) indirme + Excel/PDF/DOCX'i düz metne çevirme (openpyxl ve pypdf yalnız pylib/ altına kurulur)
ssh truba 'cd /arf/scratch/avural/kif13a_nakavt && sbatch scripts/01_fetch_inspect.slurm'
# 2) sinonim dahil gnomAD çekimi + bütün hesaplar
ssh truba 'cd /arf/scratch/avural/kif13a_nakavt && sbatch scripts/02_analysis.slurm'
# 3) küçük çıktıları geri al
scp 'truba:/arf/scratch/avural/kif13a_nakavt/interim/nakavt/*' data/interim/nakavt/
scp truba:/arf/scratch/avural/kif13a_nakavt/raw/MANIFEST.jsonl data/interim/nakavt/SOURCES_MANIFEST.jsonl
```

Login düğümünde (arf-ui) hiçbir Python ya da hesap çalıştırılmadı; orada yalnız görüntüleme ve iş
yönetimi komutları (`ls`, `mkdir`, `cat`, `head`, `grep`, `sed -n`, `zcat`, `scp`, `sbatch`, `squeue`) kullanıldı. Scratch dönemsel silinir; ham dosyalar (yaklaşık 190 MB) yeniden
indirilebilir, sha256'ları manifestte.

Üretilen dosyalar (`data/interim/nakavt/`):

| Dosya | İçerik |
|---|---|
| `KIF13A_expected_observed_20261003.tsv` | kohort başına beklenen ve gözlenen |
| `nakavt_summary_20261003.json` | bütün sayılar (F uydurmaları, ampirik bantlar, Chundru, popEVE) |
| `KIF13A_gnomad_lof_by_pop_20261003.tsv` (ve KIF1A, KIF13B) | gnomAD v4.1 LoF, popülasyon bazlı AF ve homozigot, sınıf bayrakları |
| `KIF13A_rgc_hom_candidates_20261003.tsv` | RGC homozigot frekansına uyan gnomAD varyantları |
| `gh_ST16_target_genes_20261003.tsv`, `rgc_ST2_ST6_target_genes_20261003.tsv` | hedef genlerin kohort satırları |
| `rgc_band_genes_KIF13A_20261003.tsv` | KIF13A frekans bandındaki toleranslı ve DD genleri (ampirik test) |
| `chundru_gene_level_20261003.tsv` | Chundru gen düzeyi gözlenen/beklenen, 16.655 gen |
| `KIF1A_calibration_popEVE_20261003.tsv` | KIF1A kalibrasyon setleri + popEVE |
| `ddg2p_target_genes_20261003.tsv` | KIF1A'nın iki DDG2P kaydı (KIF13A ve KIF13B DDG2P'de yok) |
| `SOURCES_MANIFEST.jsonl` | her indirmenin URL, UTC zaman, bayt, sha256 kaydı |
| `logs/` | indirme ve hesap işlerinin çıktıları |
