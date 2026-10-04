# De novo ekseni: KIF13A trio kohortlarinda vuruluyor mu

**Olcum tarihi: 16 Eylul 2026.** Uretici: `scripts/denovo_axis.py`. Tablolar: `data/interim/denovo/`.
Ham cekimler ve saglama toplamlari: `data/raw/denovo/` ve `data/raw/PROVENANCE.jsonl`.
Dogrulayamadigim her sey acikca **DOGRULANMADI** olarak isaretli. Canli getirilmemis hicbir PMID,
DOI, katilim numarasi veya URL bu belgede yer almiyor.

**Sorulan soru:** KIF13A'nin kabul edilmis bir hastaligi yok. Insandaki uc zayif gozlemin ucu de
fonksiyon kaybi sekilli, ama gnomAD KIF13A'nin fonksiyon kaybini ClinGen'in haployetmezligi
kanitlanmis 12 norogelisimsel geninin hepsinden daha iyi tolere ettigini soyluyor (gozlenen/beklenen
LoF 0,372'ye karsi 0,037-0,287). Yani agir, tam penetransli haployetmezlik zaten olasi degil.
Ayakta kalan hipotez **hafif veya eksik penetransli doz duyarliligi**. Bunu kesen test sudur:
norogelisimsel bozukluk triolarinda de novo varyantlar. KIF13A'nin bir kopyasini kaybetmek
gelisimsel bozukluk yapiyorsa, KIF13A buyuk trio kohortlarinda de novo fonksiyon kaybi
varyantlariyla **sanstan daha sik** vurulmalidir.

---

## 0. Tek paragrafta cevap

**Hukum: KARSI KANIT** (guclu bicimi icin kesin, hafif bicimi icin yaklasik 2,6 kata kadar).
Ulasilabilen her de novo kataloguna gore KIF13A, etkilenmis hicbir bireyde tek bir de novo fonksiyon
kaybi varyanti tasimiyor: Kaplanis 2020'nin 31.058 gelisimsel bozukluk triosunda 0 (beklenen 0,367),
Gene4Denovo2'nin 130.439 bireylik toplaminda 0 (beklenen 1,16), denovo-db'de 0; bulunan **tek** de
novo LoF, etkilenmemis bir kardeste (kontrol) ve son ekzonda, yani NMD kacis bolgesindedir. Ayni
boru hattinda pozitif kontroller calisiyor (ANKRD11 148 kat, CHD8 47 kat zenginlesme) ve KIF1A
tam beklenen deseni veriyor (53 de novo missense, 0 de novo LoF), dolayisiyla sifir bir olcum
basarisizligi degil, gercek bir sonuctur; bu kadar kohort ANKRD11 veya CHD8 tipi bir doz
duyarliligini kesin olarak disliyor, ama yaklasik 2,6 katin altindaki hafif bir zenginlesmeyi
hala disleyemiyor.

**En onemli tek sayi: 31.058 gelisimsel bozukluk triosunda gozlenen de novo LoF sayisi 0, beklenen
0,367.** Ikinci en onemlisi: butun kataloglarin toplaminda (yaklasik 97.700 etkin trio) yine 0,
beklenen 1,16, bu da %95 guvenle 2,6 katin uzerindeki her zenginlesmeyi disliyor.

---

## 1. Yontem

### 1.1 Kaynaklar (hepsi 16 Eylul 2026'da canli dogrulandi)

| Kaynak | Kapsam | Uc nokta | Durum |
|---|---|---|---|
| Kaplanis 2020 *Nature*, PMID 33057194 | 31.058 gelisimsel bozukluk triosu (DDD + GeneDx + Radboud), 45.221 de novo varyant | `ebi.ac.uk/europepmc/webservices/rest/PMC7116826/supplementaryFiles` | acik, indirildi |
| Gene4Denovo2, PMID 41024706 | 1.626.050 DNM, 130.439 birey, 96 fenotip | `java.genemed.tech/admin-api/backend/data/hg38/...` | acik, API ile cekildi |
| denovo-db v1.6.1, PMID 27907889 | 415.515 varyant, 49 calisma, 11.699 tekil orneklem (SSC disi kisim) | `denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz` | acik, indirildi |
| gnomAD v2.1.1 kisitlama | gen basina trinukleotid mutasyon hizlari (`mu_lof`, `mu_mis`, `mu_syn`) | `storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/...` | acik, indirildi |
| SFARI Gene, surum 2026 Q2 | otizm gen skorlari | `gene.sfari.org/database/human-gene/<GEN>` | acik, okundu |
| ClinVar (E-utilities) | VCV000442594'un kalitim alani | `eutils.ncbi.nlm.nih.gov` | acik, okundu |
| DECIPHER hasta CNV'leri | de novo CNV kalitimi | `deciphergenomics.org/about/downloads/data` | **HTTP 403, hesap gerekiyor, DURULDU** |
| denovo-db SSC kismi | Simons Simplex Collection | `denovo-db.gs.washington.edu/denovo-db.ssc-samples.variants.tsv.gz` | **kullanim kisiti beyan edilmis, INDIRILMEDI** |

Kural geregi hesap, odeme veya veri erisim sozlesmesi isteyen iki kaynakta durduk ve bunu burada
belirtiyoruz. DECIPHER'in gen sayfasi acilyor ama CNV tablosunu besleyen API anahtar istiyor, yani
de novo CNV kalitim bilgisi bu belgede **elde edilemedi**.

### 1.2 Sayim nasil yorumlanabilir hale getirildi

Ham sayi tek basina hicbir sey ifade etmez. Beklenen de novo sayisi Samocha cercevesinin gen basina
mutasyon hizlariyla hesaplandi:

```
beklenen(sinif) = 2 * trio_sayisi * mu(sinif)
```

`mu` degerleri gnomAD v2.1.1 kisitlama dosyasindan, o dosyanin **bizim MANE transkriptimizle ayni**
transkript uzerinde raporladigi satirdan alindi. Poisson p degeri ust kuyruktan, log uzayinda
toplanarak hesaplandi; sade `1 - CDF` formu pozitif kontrollerde tasarak sahte bir `p = 0`
uretiyordu (ANKRD11'de 111 gozlem, 0,31 beklenti).

### 1.3 Uzunluk ve transkript saglamasi (1805'e karsi 1770 tuzagi)

gnomAD kisitlama dosyasinin CDS uzunlugu, protein uzunlugunun **bagimsiz** bir dogrulamasidir:

| Gen | Hiz tablosundaki transkript | MANE (config) | Tutuyor mu | CDS | Protein | config | Tutuyor mu |
|---|---|---|---|---|---|---|---|
| KIF13A | ENST00000259711 | ENST00000259711 | evet | 5415 nt | 1805 aa | 1805 | evet |
| KIF1A | ENST00000498729 | ENST00000498729 | evet | 5373 nt | 1791 aa | 1791 | evet |
| KIF13B | ENST00000524189 | ENST00000524189 | evet | 5478 nt | 1826 aa | 1826 | evet |

Yani bu belgedeki butun KIF13A protein konumlari **1805 aa'lik MANE numaralandirmasindadir**.
1770 aa'lik izoform (Q9H1H9-2) hicbir yerde kullanilmadi.

### 1.4 Izoform ve yapi (build) tuzaklarinin nasil kapatildigi

Uc kaynak da farkli numaralandirma kullaniyor ve ucu de acikca donusturuldu:

- **Kaplanis 2020** tablosunda yapi (assembly) yazmiyor. Varsaymak yerine **test edildi**: dort
  KIF13A satirinin her biri hem GRCh38 hem GRCh37 uzerinde VEP'ten gecirildi ve tablonun belirttigi
  sonuc sinifiyla hangi yapinin tuttugu sayildi. Sonuc **4/4 GRCh37**, 0/4 GRCh38. Sonra hepsi
  Ensembl ile GRCh38'e tasindi ve MANE uzerinde yeniden anotlandi.
  - Bagimsiz dogrulama: bu tasima, Gene4Denovo2'nin ayni varyantlar icin bagimsiz olarak tuttugu
    hg38 koordinatlarinin **birebir aynisini** uretti (ornegin chr6:17873402A>T, p.Asp65Glu).
- **Gene4Denovo2** iki anotasyon sutunu tasiyor. `AAChange_refGene` KIF13A icin NM_001105567 /
  NM_001105568 kullaniyor ve numaralari **13 kalinti dusuk** veriyor (ornegin ayni varyant icin
  p.S1220L). `AAChange_mrefGene` ise MANE RefSeq'i NM_022113'u kullaniyor (p.S1233L). Bu belgede
  yalniz `mrefGene` okundu.
- **denovo-db** her varyanti **her transkript icin ayri bir satir** olarak yaziyor. KIF13A'nin 101
  ham satiri aslinda 25 tekil varyanttir. Protein konumu yalniz NM_022113 satirindan alindi;
  NM_001105567/8 satirlari ayni 13 kalintilik kaymayi gosteriyor ve bu, tuzagin ucuncu bagimsiz
  teyididir.

### 1.5 Ortusme (ayni bireyin iki kez sayilmasi)

Kaynaklar birbirini kapsiyor. KIF13A birlesim tablosu bu yuzden **(GRCh38 konum, ref, alt)**
uzerinden tekillestirildi ve her satirda hangi kaynak, calisma, PMID ve orneklem kimliginin
destekledigi yazildi (`data/interim/denovo/KIF13A_denovo_union_20260916.tsv`).

---

## 2. Kaynak kaynak ham satirlar

### 2.1 Kaplanis 2020 (PMID 33057194), 31.058 gelisimsel bozukluk triosu

URL: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7116826/supplementaryFiles`
(Europe PMC, acik erisim, 16 Eyl 2026'da indirildi, 7.841.523 bayt).
Supplementary Table 1 = 45.221 de novo varyant (GDX 26.679, DDD 15.211, RUMC 3.331).
Supplementary Table 2 = 19.654 gen, sayfa adi `kaplanis_samocha_denovoWEST_res`.
Sagalama: Table 1'den saydigimiz bes genin sinif dagilimi, Table 2'nin yayinlanmis sayilariyla
**birebir tutuyor**.

**KIF13A satirlari (hepsi, dort tane):**

| id | calisma | GRCh37 | GRCh38 | tablo sinifi | MANE HGVSc | MANE HGVSp | konum | domain |
|---|---|---|---|---|---|---|---|---|
| 98963 | GDX | 6:17873633 A>T | chr6:17873402 A>T | missense | c.195T>A | p.Asp65Glu | **65** | **motor** |
| DDD13k.07429 | DDD | 6:17834250 A>G | chr6:17834019 A>G | missense | c.1208T>C | p.Ile403Thr | 403 | CC1 |
| DDD13k.09819 | DDD | 6:17788092 G>A | chr6:17787861 G>A | sinonim | c.3276C>T | p.Asn1092= | 1092 | kuyruk |
| DDD13k.09076 | DDD | 6:17764668 G>A | chr6:17764437 G>A | sinonim | c.5091C>T | p.Gly1697= | 1697 | kuyruk |

Transkript: ENST00000259711.11 = NM_022113.6 (MANE Select), 1805 aa.
**De novo fonksiyon kaybi varyanti sayisi: 0.** Cerceve kaymasi 0, stop_gained 0, splice_acceptor 0,
splice_donor 0.

Kaplanis'in kendi istatistikleri KIF13A icin: `PTV_enrichment_p` = 0,448, `denovoWEST_p_full` = 0,709,
`significant` = False, `diagnostic_category` = none, `s_het` = 0,0337.

Tablolar: `data/interim/denovo/kaplanis_five_gene_dnms_20260916.tsv` (ham satirlar),
`kaplanis_KIF13A_mane_rows_20260916.tsv` (MANE anotasyonu),
`kaplanis_per_gene_poisson_20260916.tsv` (Poisson).

### 2.2 Gene4Denovo2 (PMID 41024706)

Arayuz: `https://genemed.tech/gene4denovo2/#/genemed/search?reference=hg38&queryby=geneSymbol&terms=KIF13A`
API (arayuzun kendi cagirdigi uc nokta, 16 Eyl 2026):
`POST https://java.genemed.tech/admin-api/backend/data/hg38/gene4denovo2_all/geneSymbol/page`
govde: `{"reference":"hg38","queryby":"geneSymbol","terms":"KIF13A","page":1,"pageNo":1,"currentPage":1,"pageSize":2000,"filterParam":"[]","mode":"variant"}`
Gen duzeyi icin ayni govde `mode":"gene"` ve tablo `gene4denovo2_Candidate_Gene2`.
Beyan edilen kapsam (ana sayfa): 1.626.050 DNM, 133.738'i kodlayan bolgede, 130.439 birey, 96 fenotip.

**KIF13A: toplam 174 de novo varyant.** Islevsel dagilim: 130 intronik, 26 intergenik, 9 kodlayan
missense, 4 sinonim, 3 UTR3, 1 cerceve kaymasi, 1 upstream. Yani **kodlayan yalniz 14 tanesi**.

KIF13A'nin butun kodlayan satirlari (GRCh38, MANE NM_022113):

| GRCh38 | MANE degisimi | konum | domain | sinif | fenotip | calisma / PMID | orneklem |
|---|---|---|---|---|---|---|---|
| chr6:17873402 A>T | c.T195A p.D65E | 65 | **motor** | missense | UDD | Kaplanis 2020 / 33057194 | 98963 |
| chr6:17852019 C>T | c.G518A p.R173Q | 173 | **motor** | missense | ASD | Yuen 2017 / 28263302 | AU1894304 |
| chr6:17852019 C>T | c.G518A p.R173Q | 173 | **motor** | missense | ASD | Cirnigliaro 2023 / 37506195 | mAGRE5564 |
| chr6:17850329 C>G | c.G711C p.Q237H | 237 | **motor** | missense | ALS | Chesi 2013 / 23708140 | Chesi50 |
| chr6:17837037 G>A | c.C996T p.A332A | 332 | motor | sinonim | ASD | Feliciano 2019 / 31452935 | SP0018753 |
| chr6:17837037 G>A | c.C996T p.A332A | 332 | motor | sinonim | ASD | Zhou 2022 / 35982159 | SP0018753 |
| chr6:17834019 A>G | c.T1208C p.I403T | 403 | CC1 | missense | UDD | Kaplanis 2020 / 33057194 | DDD13k.07429 |
| chr6:17787861 G>A | c.C3276T p.N1092N | 1092 | kuyruk | sinonim | UDD | Kaplanis 2020 / 33057194 | DDD13k.09819 |
| chr6:17785566 G>A | c.C3437T p.A1146V | 1146 | kuyruk | missense | ASD | Koko 2025 / 39954678 | SP0225804 |
| chr6:17780878 G>A | c.C3698T p.S1233L | 1233 | kuyruk | missense | **Kontrol** | Zhou 2022 / 35982159 | SP0130215 |
| chr6:17779096 T>C | c.A3943G p.T1315A | 1315 | kuyruk | missense | ASD | Koko 2025 / 39954678 | SP0456506 |
| chr6:17764612 G>A | c.C4916T p.S1639L | 1639 | kuyruk | missense | **Kontrol** | Koko 2025 / 39954678 | SP0419078 |
| **chr6:17764578 TTTG>-** | **c.4947_4950del p.N1649Kfs\*3** | **1649** | **kuyruk** | **LoF** | **Kontrol** | Koko 2025 / 39954678 | SP0246168 |
| chr6:17764437 G>A | c.C5091T p.G1697G | 1697 | kuyruk | sinonim | UDD | Kaplanis 2020 / 33057194 | DDD13k.09076 |

Iki nokta kritik:

1. **Tek de novo LoF, etkilenmemis bir kardeste.** `p.N1649Kfs*3`, SPARK'in kontrol (saglikli
   kardes) kolunda. Ustelik **ekzon 39'da**. Ensembl'den canli dogrulandi: ENST00000259711 tam
   **39 ekzon** ve 1805 aa tasiyor, yani ekzon 39 son ekzondur. Fonksiyon kaybi ekseninin
   belgeledigi uzere bu son ekzon kodon 1528'den itibaren proteinin son 278 kalintisini (%15,5)
   kodluyor ve oradaki kesilmeler NMD'den kacar; kalinti 1649 bu bolgenin icindedir. Yani bu
   varyant hem hasta degil hem de muhtemelen gercek bir fonksiyon kaybi degil. Iki ayri nedenle
   doz kaniti sayilmaz.
2. **A332A ayni bireyde (SP0018753) iki kez.** Feliciano 2019 ve Zhou 2022 ayni SPARK bireyini
   rapor ediyor. R173Q ise iki farkli kimlikte (AU1894304 ve mAGRE5564); AGRE ve MSSNG aile
   kohortlari ortustugu icin **bunun bir mi iki mi birey oldugu bu veriden ayirt edilemiyor**
   (DOGRULANMADI).

**Gene4Denovo2'nin kendi gen duzeyi ozeti (ALL Disorders satiri, yani kontroller haric):**

| Gen | LoF | Dmis | Tmis | Syn | Pval | Qval | Guven |
|---|---|---|---|---|---|---|---|
| **KIF13A** | **0** | 4 | 3 | 4 | **0,975374** | 0,951588 | Grey |
| KIF1A | 2 | 64 | 7 | 6 | 0 | 0 | Red |
| KIF13B | 2 | 2 | 10 | 4 | 0,835271 | 0,942419 | Grey |
| ANKRD11 | 144 | 15 | 17 | 14 | 0 | 0 | Red |
| CHD8 | 74 | 20 | 11 | 6 | 0 | 0 | Red |

KIF13A'nin fenotip kirilimi: ASD p = 0,900 (LoF 0), UDD p = 0,931 (LoF 0), ALL Disorders p = 0,975
(LoF 0). Kontrol satiri: LoF 1.

**KIF1A'nin 2 "LoF" satiri anotasyon hatasi:** biri chr2:240788192 `-`>`GTC`, yani 3 bazlik bir
**cerceve ici insersiyon**; ANNOVAR bunu `stopgain` diye etiketlemis, oysa Kaplanis'in kendi
tablosunda ayni varyant `inframe_insertion` olarak geciyor ve KIF1A'nin de novo LoF sayisi 0'dir.
Digeri PMID 38105687'den (Lin 2023, alti pediatrik KIF1A hastasi) bir splice bolgesi varyanti.
Bu, Gene4Denovo2'nin LoF etiketinin tek basina guvenilmez oldugunu gosterir ve neden birincil
istatistigi Kaplanis uzerinden kurdugumuzu aciklar.

Tablolar: `data/interim/denovo/gene4denovo2_five_gene_dnms_20260916.tsv`,
`gene4denovo2_gene_level_20260916.tsv`.

### 2.3 denovo-db v1.6.1 (PMID 27907889)

URL: `https://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz`
(16 Eyl 2026, 7.734.071 bayt, dosya ici surum etiketi `denovo-db.v.1.6.1`).
Kapsam (bizim saydigimiz): 415.515 varyant satiri, 49 calisma, 11.699 tekil orneklem kimligi.
Koordinatlar GRCh37.

**SSC kismi bilerek indirilmedi:** indirme sayfasi Simons Simplex Collection ve Simons VIP
verilerinin kullanimini otizm ve iliskili gelisimsel bozukluk arastirmasiyla sinirlayan bir
kullanim kisiti beyan ediyor. Kural geregi durduk ve bunu rapor ediyoruz.

**KIF13A: 101 ham satir, 25 tekil varyant, bunlarin 3'u kodlayan.** (Ham satir sayisi
yaniltici, cunku denovo-db her transkript icin ayri satir yazar.)

| GRCh37 | GRCh38 | transkript | cDNA | protein | sinif | fenotip | calisma / PMID | orneklem |
|---|---|---|---|---|---|---|---|---|
| chr6:17788092 G>A | chr6:17787861 | NM_022113.5 (MANE) | c.3276C>T | p.(N1092=) | sinonim | developmentalDisorder | DDD_2017 / 28135719 | DDD4K.03538 |
| chr6:17850560 C>G | chr6:17850329 | NM_022113.5 (MANE) | c.711G>C | p.(Q237H) | missense | amyotrophic_lateral_sclerosis | Chesi2015 (denovo-db etiketi) / 23708140 | Chesi50 |
| chr6:17852250 C>T | chr6:17852019 | NM_022113.5 (MANE) | c.518G>A | p.(R173Q) | missense | autism | Yuen2017 / 28263302 | AU1894304 |

**De novo fonksiyon kaybi varyanti sayisi: 0** (ne hastada ne kontrolde).
Kalan 22 tekil varyantin hepsi intronik veya 3'-UTR: 54 satir "mixed", 34 satir otizm fenotipinde.

Ayni tabloda NM_001105567.2 ve NM_001105568.2 satirlari ayni varyantlar icin p.(N1092=) yerine
p.(N1079=) veriyor, yani 13 kalintilik izoform kaymasinin bagimsiz teyidi.

Tablo: `data/interim/denovo/denovodb_five_gene_dnms_20260916.tsv`.

### 2.4 SFARI Gene (surum 2026 Q2)

| Gen | Durum |
|---|---|
| **KIF13A** | **"KIF13A is not in the SFARI:Gene database"** (acik, sozu gecen ifade sayfanin kendisinde) |
| KIF1A | Gene Score 2S, Strong Candidate, Syndromic |
| KIF13B | Gene Score 2, Strong Candidate, toplam 7 rapor |
| ANKRD11 | Gene Score 1S, High Confidence, Syndromic |
| CHD8 | Gene Score 1S, High Confidence, Syndromic |

Not: paralog KIF13B otizm listesinde, hedefimiz KIF13A degil.

### 2.5 De novo CNV'ler ve 253 kb'lik delesyon

- **ClinVar VCV000442594** (KIF13A'yi tek basina silen 253 kb'lik delesyon, etkilenmis cocuk):
  VCV XML'i canli okundu, **`Origin` alani "unknown"**. Yani bu delesyonun de novo oldugu
  **kayitli degil**. Projenin en temiz tek doz gozlemi, de novo argumanina hic giremiyor.
  Kayitli fenotip: "Delayed speech and language development"; `ReviewStatus`: "no assertion
  criteria provided"; yontem: mikroarray ve FISH.
- **DECIPHER hasta CNV'leri:** `deciphergenomics.org/about/downloads/data` HTTP 403 donuyor ve
  CNV tablosunu besleyen API anahtar istiyor. Hesap acmadik. Dolayisiyla **de novo CNV kalitim
  verisi bu belgede elde edilemedi (DOGRULANMADI)**.
- Gene4Denovo2 ve denovo-db yapisal varyant tasimiyor; ikisi de SNV ve kucuk indel veritabani.

---

## 3. Bes gen karsilastirmasi

Birincil tablo. Kaplanis 2020, 31.058 gelisimsel bozukluk triosu, hepsi ayni boru hattindan gecti.
Beklenen degerler `2 x 31.058 x mu`, mu gnomAD v2.1.1'den ve her gen icin MANE transkripti uzerinden.

| Gen | Rol | goz. LoF | bekl. LoF | p (LoF) | goz. missense | bekl. missense | p (missense) | goz. sinonim | bekl. sinonim | denovoWEST p | Kaplanis anlamli mi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **KIF13A** | **hedef** | **0** | **0,367** | **1** | **2** | **3,163** | **0,824** | 2 | 1,253 | **0,709** | **Hayir** |
| KIF1A | paralog, missense mekanizmasi | 0 | 0,298 | 1 | **53** | 4,806 | **2,9e-36** | 0 | 2,260 | 5,2e-21 | Evet |
| KIF13B | paralog, null beklenir | 2 | 0,337 | 0,046 | 6 | 3,763 | 0,179 | 2 | 1,618 | 0,286 | Hayir |
| ANKRD11 | pozitif kontrol (ClinGen HI = 3) | **111** | 0,308 | **8,2e-238** | 17 | 7,065 | 0,0011 | 7 | 4,098 | 1e-14 | Evet |
| CHD8 | pozitif kontrol (ClinGen HI = 3) | **19** | 0,497 | **8,8e-24** | 9 | 5,138 | 0,078 | 0 | 1,577 | 4,3e-16 | Evet |

Ayni bes gen, butun kataloglarin toplaminda (Gene4Denovo2 "ALL Disorders", etkin trio sayisi
sinonim de novo'lardan kalibre edildi, bkz. 4.2):

| Gen | goz. LoF | bekl. LoF | zenginlesme | Poisson p | %95 ust sinir (kat) | G4D2 kendi p | G4D2 guven | SFARI |
|---|---|---|---|---|---|---|---|---|
| **KIF13A** | **0** | **1,156** | **0,0x** | 1 | **2,6x** | 0,975 | Grey | listede degil |
| KIF1A | 2 (ikisi de anotasyon hatasi) | 0,937 | 2,1x | 0,241 | 6,7x | 0 | Red | 2S |
| KIF13B | 2 | 1,060 | 1,9x | 0,286 | 5,9x | 0,835 | Grey | 2 |
| ANKRD11 | 144 | 0,970 | **148x** | 9,0e-253 | 103x | 0 | Red | 1S |
| CHD8 | 74 | 1,564 | **47x** | 1,6e-94 | 57x | 0 | Red | 1S |

**Bu iki tablonun okunmasi:**

1. **Boru hatti calisiyor.** ANKRD11 ve CHD8 iki kohortta da ezici bicimde zenginlesmis. Sifir bir
   yontem hatasi degil.
2. **KIF1A tam beklenen deseni veriyor:** 53 de novo missense, **0 de novo LoF**. ClinGen'in KIF1A'ya
   haployetmezlik skoru 0 vermesi ve hastaligin missense kaynakli olmasi burada dogrudan goruluyor.
   Bu, "gen hastalik yapiyor ama LoF ile degil" durumunun nasil gorundugunun kanitidir ve
   KIF13A'nin sifirindan **acikca farklidir** (KIF13A'da ne LoF ne missense sinyali var).
3. **KIF13A hicbir sinifta sinyal vermiyor.** LoF 0, missense beklentinin **altinda** (2'ye karsi
   3,16), sinonim beklendigi gibi. Uc bagimsiz istatistik ayni yone bakiyor: bizim Poisson'umuz,
   Kaplanis'in denovoWEST'i (p = 0,709) ve Gene4Denovo2'nin kendi testi (p = 0,975).
4. **Durust bir rahatsizlik:** "null beklenen" paralog KIF13B, Kaplanis'te nominal olarak
   KIF13A'dan **daha zenginlesmis** gorunuyor (2 stop_gained, p = 0,046). Bu coklu test
   duzeltmesinden gecmez (19.654 gen icin Bonferroni esigi 2,5e-6) ve Kaplanis'in kendi denovoWEST'i
   0,286 veriyor. Ama kaydedilmeye deger: tasarimda null olmasi beklenen gen bile hedef genden
   daha iyi bir nominal p uretiyor, bu da KIF13A'nin sifirinin ne kadar bos oldugunu gosterir.

---

## 4. Gozlenen ve beklenen: hesap ve varsayimlari

### 4.1 Kaplanis kohortu (olculmus, N kesin biliniyor)

**Formul:** `beklenen = 2 * N_trio * mu(sinif)`. Iki carpani her cocugun iki ebeveyn kromozomu
tasimasindan gelir.

**Girdiler (hepsi dosyadan, elle girilmedi):**

| Nicelik | Deger | Kaynak |
|---|---|---|
| N_trio | 31.058 | Kaplanis 2020 ozeti, Europe PMC'den canli okundu |
| mu_lof (KIF13A) | 5,9159e-06 | gnomAD v2.1.1, ENST00000259711 satiri |
| mu_mis (KIF13A) | 5,0924e-05 | ayni satir |
| mu_syn (KIF13A) | 2,0177e-05 | ayni satir |

**Sonuc:**

| Sinif | Beklenen | Gozlenen | Poisson p (tek yonlu, P(X >= gozlenen)) |
|---|---|---|---|
| LoF | 2 x 31.058 x 5,9159e-06 = **0,3675** | **0** | 1,000 |
| missense | 2 x 31.058 x 5,0924e-05 = **3,1632** | 2 | 0,824 |
| sinonim | 2 x 31.058 x 2,0177e-05 = **1,2533** | 2 | 0,357 |

Sifir gozlem icin tek yonlu %95 Poisson ust siniri **2,996 olay**, yani gercek LoF hizinin
beklentinin **8,15 katindan** buyuk olmasi %95 guvenle disleniyor.

**Varsayimlar, acikca:**

1. Mutasyon hizlari gnomAD v2.1.1'in trinukleotid modelinden geliyor; bunlar Samocha cercevesinin
   hizlaridir ama gnomAD v2 (GRCh37, ekzom) uzerinde hesaplanmistir. Kaplanis'in kendi
   `expected_full` sutunu farkli bir buyukluktur (denovoWEST'in agirliklandirilmis skoru, ham sayim
   degil), o yuzden karsilastirilabilir bir "beklenen sayim" icin mu'lar kullanildi.
2. Trio sayisinin tamami her gende tam kapsanmis sayiliyor. Gercekte ekzom kapsami gene gore
   degisir; bu, beklenen degeri **yukari** dogru yanlilar, yani sifir gozlemi olduğundan daha
   carpici gosterebilir. Karsi yonde bir yanlilik yok.
3. Poisson bagimsizlik varsayiyor (mutasyonel sicak nokta yok). LoF sinifi icin makul.
4. Baba yasi dagilimi kohortlar arasi farkli olabilir; bu butun genleri ayni yonde etkiler,
   dolayisiyla gen karsilastirmasini bozmaz.

### 4.2 Toplam kohort (etkin N tahmin edildi, olculmedi)

Gene4Denovo2 kac **trio** oldugunu gen basina vermiyor, yalnizca 130.439 birey diyor. Bunun yerine
etkin trio sayisi **sinonim de novo varyantlardan kalibre edildi**: sinonim varyantlar notrdur,
dolayisiyla sayilari yalniz maruziyeti (kac trio dizilendi) olcer.

```
N_etkin = (bes genin ALL Disorders sinonim toplami) / (2 x ayni bes genin mu_syn toplami)
        = 34 / (2 x 1,73972e-04)
        = 97.717   (%95 Poisson araligi 71.880 - 130.095)
```

Bu tahmin Gene4Denovo2'nin kendi beyan ettigi 130.439 bireyle tutarlidir (aralik onu iceriyor) ve
tahmin edicinin dogrulugu pozitif kontrollerde gorulebilir: ANKRD11 148x, CHD8 47x zenginlesme
veriyor, yani dogru buyukluk sirasinda.

**KIF13A sonucu:** beklenen LoF = 2 x 97.717 x 5,9159e-06 = **1,156**; gozlenen **0**;
sifir gozlemin %95 ust siniri 2,996 olay, yani **2,6 katin uzerindeki her zenginlesme disleniyor**.

**Bu bolumun varsayimlari (4.1'dekilere ek olarak):** sinonim cagri duyarliliginin bes gende ayni
oldugu, Gene4Denovo2'nin sinonim de novo'lari eksiksiz topladigi ve "ALL Disorders" kumesindeki
fenotip karisiminin (UDD, ASD, NDD, ID, SCZ, CHD, EE ve digerleri) KIF13A icin anlamli bir test
kumesi oldugu. Ucuncusu en tartismali olani: kume saf norogelisimsel degil. Ama bu da sifiri
**kolaylastirir** degil zorlastirir, cunku daha genis fenotip daha cok firsat demektir.

Tablolar: `data/interim/denovo/gene4denovo2_aggregate_poisson_20260916.tsv`,
`power_sensitivity_20260916.tsv`.

---

## 5. Motor domain ozel kontrolu

Projenin aday listesi (katman A ve B) motor domain uzerine kurulu, dolayisiyla soru yerinde:
bulunan de novo varyantlarin kaci 5-352 arasina dusuyor?

Motor domain 5-352 = **348 kalinti** = KIF13A'nin **%19,3'u** (gorev metnindeki 352/1805 = %19,5
yaklasimiyla ayni buyukluk).

Uc kaynagin tekillestirilmis birlesiminde **9 kodlayan (missense + LoF) de novo varyant** var;
bunlarin **3'u motor domainde**:

| MANE degisimi | konum | domain | sinif | fenotip | destekleyen kaynak |
|---|---|---|---|---|---|
| p.Asp65Glu | 65 | **motor** | missense | UDD / gelisimsel bozukluk | Kaplanis 2020 + Gene4Denovo2 |
| p.Arg173Gln | 173 | **motor** | missense | ASD / otizm | Gene4Denovo2 + denovo-db |
| p.Gln237His | 237 | **motor** | missense | **ALS** | Gene4Denovo2 + denovo-db |
| p.Ile403Thr | 403 | CC1 | missense | UDD / gelisimsel bozukluk | Kaplanis 2020 + Gene4Denovo2 |
| p.Ala1146Val | 1146 | kuyruk | missense | ASD | Gene4Denovo2 |
| p.Ser1233Leu | 1233 | kuyruk | missense | Kontrol | Gene4Denovo2 |
| p.Thr1315Ala | 1315 | kuyruk | missense | ASD | Gene4Denovo2 |
| p.Ser1639Leu | 1639 | kuyruk | missense | Kontrol | Gene4Denovo2 |
| p.Asn1649Lysfs\*3 | 1649 | kuyruk | **LoF** | **Kontrol** | Gene4Denovo2 |

**Yorum, fazlasiyla degil:** 9 varyanttan 3'u (%33) motorda, uzunluk oranina gore beklenen %19,3.
Bu bir zenginlesme testi degildir cunku sayilar cok kucuk: uzunluk oranina gore **en az bir** motor
isabeti gorme olasiligi zaten 0,855'tir. Ustelik ucundan biri (Q237H) bir **ALS** triosundan
geliyor, norogelisimsel degil. Yani motor domainde de novo varyant **var**, ama bu varyantlar ne
sayica ne dagilim olarak sanstan ayrilmiyor ve ucu de fonksiyon kaybi degil, missense.

Tablo: `data/interim/denovo/KIF13A_denovo_union_20260916.tsv`.

---

## 6. Bilinen uc KIF13A gozlemi bu kohortlarda var mi

| Gozlem | Gosterim | Kohortlarda | Not |
|---|---|---|---|
| Lam 2021 mozaik cerceve kaymasi (PMID 33526817) | c.2357dupA, kalinti ~786 | **HAYIR** | Hicbir kaynakta 786 civarinda de novo yok. Zaten **postzigotik mozaik**, yani klasik de novo trio cagrisina girmesi beklenmez. |
| 253 kb 6p22.3 delesyonu | VCV000442594 | **HAYIR** | De novo katalogu degil, CNV. ClinVar'daki `Origin` alani **"unknown"**, yani de novo oldugu kayitli degil. |
| p.Lys1122Glu kuyruk missense VUS'u | VCV002570672, c.3364A>G, kalinti 1122 | **HAYIR** | Birlesim tablosunda 1122'de de novo yok; en yakin kodlayan de novo 1146 (p.A1146V, ASD) ve 1092 (sinonim). |

Ucunun de trio kohortlarinda cikmamasi **beklenen** bir sonuctur ve bir celiski degildir: ilki
mozaik, ikincisi CNV, ucuncusu arastirma gonderimi tekil bir hasta. Ama sonuc su: **projenin uc
insan gozleminden hicbiri de novo kanitiyla desteklenmiyor**, ve ikincisinin kalitimi bilinmedigi
icin "de novo delesyon" diye anilmasi da yanlis olur.

Tablo: `data/interim/denovo/known_observations_crosscheck_20260916.tsv`.

---

## 7. Bu eksen ne kadar ileri goturulebilir

### 7.1 Kanitin izin verdigi **en guclu** cumle

> Ulasilabilen butun kamuya acik de novo kataloglarinda (Kaplanis 2020'nin 31.058 gelisimsel bozukluk
> triosu, Gene4Denovo2'nin 130.439 bireyi ve denovo-db v1.6.1 dahil, sinonim de novo'larla kalibre
> edilmis yaklasik 97.700 etkin trio), KIF13A'da etkilenmis bir bireyde **tek bir de novo fonksiyon
> kaybi varyanti bile yoktur**; beklenen 1,16'ya karsilik gozlenen 0'dir, bulunan tek de novo LoF
> etkilenmemis bir kardeste ve NMD kacis bolgesindedir, ve ayni boru hatti ayni kohortlarda ANKRD11'i
> 148 kat, CHD8'i 47 kat zenginlesmis olarak dogru sekilde tespit etmektedir. Bu, KIF13A icin
> ANKRD11 veya CHD8 tipi bir haployetmezlik mekanizmasini **disliyor** ve doz duyarliligi kalacaksa
> bunun de novo LoF zenginlesmesi uretmeyecek kadar hafif veya nadir olmasi gerektigini soyluyor
> (%95 guvenle 2,6 kattan az).

### 7.2 Asiri iddia olacak cumle

> "KIF13A'nin haployetmezligi hastalik yapmaz" ya da "de novo verisi KIF13A'nin doz duyarli
> olmadigini kanitliyor."

Neden asiri iddia: sifir gozlem, sifir etki demek degildir. Beklenen sayi zaten 1'in civarinda
oldugu icin, gercek zenginlesme 2 kat olsaydi da buyuk olasilikla 0 veya 1 gozlem gorurduk. Bu
kohortlar yalnizca **2,6 katin uzerini** disliyor, altini degil. Ayrica: (a) eksik penetransli bir
allel tanisiz kalmis ailelerde hic dizilenmemis olabilir, (b) KIF13A'nin hucre tipi ekseni fenotipin
norogelisimsel degil **ak madde / miyelin** tarafinda olabilecegini soyluyor ve bu kohortlar
gelisimsel gecikme ile secilmistir, (c) mozaik ve CNV mekanizmalari bu kataloglarda zaten
gorunmez, (d) DECIPHER'in de novo CNV verisi alinamadi.

### 7.3 Bu eksenin projenin geri kalanina etkisi

- **Haployetmezlik hipotezi bir kademe daha zayifladi.** Fonksiyon kaybi ekseni "olasi doz
  duyarliligi, eksik penetransla" demisti. De novo ekseni bunun ustune, doz duyarliligi varsa
  **2,6 kattan kucuk** olmak zorunda oldugunu ekliyor. Bu, ClinGen'in haployetmezlik genleri
  icin gordugumuz 47-148 kat araligindan iki buyukluk mertebesi uzaktir.
- **Missense hipotezi bu eksenden ne destek ne darbe aliyor.** Gozlenen de novo missense (2)
  beklenenin (3,16) altinda, ama bu sayilarla hicbir sey soylenemez. Motor domainde uc de novo
  missense var (D65E, R173Q, Q237H) ve bunlar doygunluk tablosuyla karsilastirilabilecek gercek
  insan gozlemleridir; ama ucu de fenotipi KIF13A'ya atfedilmemis vakalardan geliyor ve biri ALS.
- **Baskin-negatif hipotezi bu eksenle celismiyor.** De novo LoF'un sifir olmasi, tam boy hatali
  proteinin zehirleyici olabilecegi fikrini curutmez; tersine, "doz kaybi mekanizma degil"
  yonundeki delili guclendirdigi olcude, katman A'nin tek tutarli gerekcesi olarak baskin-negatifi
  geride birakir. Bu mekanizma hala **hic olculmedi**.

### 7.4 Bu eksende sirada ne var

1. **D65E, R173Q ve Q237H'i doygunluk tablosuyla kesistirin.** Ucu de motor domainde, ucu de gercek
   insanlarda gorulmus. Doygunluk tablosu (`KIF13A_saturation_lookup_20260916.tsv`) bunlar icin ne
   diyor, KIF1A'da karsilik kalintida hastalik varyanti var mi? Bu, tablonun ilk gercek sinavidir.
2. **R173Q'nun bir mi iki birey oldugunu cozun.** AGRE ve MSSNG ortusuyor. Iki bagimsiz otizm
   triosunda ayni kalinti degisimi anlamli olurdu; tek birey iki kez sayilmissa degildir.
3. **Kaplanis'in kendi expected sutunlariyla capraz kontrol.** Biz mu tabanli beklenen kullandik;
   denovoWEST'in `expected_full` degeri farkli bir olcek. Ikisi ayni sonuca gotururse hesap saglam.
4. **DECIPHER erisimi.** De novo CNV kalitim verisi bu eksenin tek kapali kapisi. Hesap
   gerektiriyor, yani insan kullanici acmali.
5. **Daha yeni ve daha buyuk kohortlar.** Duyarlilik tablosuna gore 200.000 trioda beklenen LoF
   2,37 olur ve sifir gozlem 1,3 katin uzerini disler, yani soru pratik olarak kapanir.

---

## 8. Uretilen dosyalar

| Dosya | Icerik |
|---|---|
| `data/interim/denovo/five_gene_comparison_20260916.tsv` | bes gen yan yana, butun kaynaklar |
| `data/interim/denovo/kaplanis_per_gene_poisson_20260916.tsv` | Kaplanis gozlenen/beklenen ve Poisson |
| `data/interim/denovo/kaplanis_five_gene_dnms_20260916.tsv` | Kaplanis ham de novo satirlari |
| `data/interim/denovo/kaplanis_KIF13A_mane_rows_20260916.tsv` | KIF13A satirlarinin MANE anotasyonu |
| `data/interim/denovo/gene4denovo2_five_gene_dnms_20260916.tsv` | Gene4Denovo2 varyant duzeyi |
| `data/interim/denovo/gene4denovo2_gene_level_20260916.tsv` | Gene4Denovo2 gen duzeyi ozet |
| `data/interim/denovo/gene4denovo2_aggregate_poisson_20260916.tsv` | toplam kohort, sinonim kalibrasyonlu |
| `data/interim/denovo/denovodb_five_gene_dnms_20260916.tsv` | denovo-db tekillestirilmis satirlar |
| `data/interim/denovo/KIF13A_denovo_union_20260916.tsv` | uc kaynagin tekillestirilmis birlesimi |
| `data/interim/denovo/power_sensitivity_20260916.tsv` | kohort buyudukce ne disleniyor |
| `data/interim/denovo/known_observations_crosscheck_20260916.tsv` | uc bilinen gozlemin taranmasi |
| `data/interim/denovo/denovo_axis_summary_20260916.json` | makine okunur ozet |

Ham cekimler (`data/raw/denovo/`, her biri `.prov.json` yan dosyasi ve `PROVENANCE.jsonl` satiri
ile): Kaplanis ek dosyalari (zip), Gene4Denovo2 API cikti (json), denovo-db non-SSC (tsv.gz),
gnomAD v2.1.1 kisitlama (txt.bgz), SFARI durumu (json), liftover onbellekleri (json).

---

## 8b. Dogrulanmis atiflar

Bu belgedeki her PMID 16 Eylul 2026'da canli getirildi (Europe PMC arama uc noktasi ya da NCBI
E-utilities `esummary`). Getirilmemis hicbir tanimlayici kullanilmadi.

| PMID | Yil | Ilk yazar | Dergi | Ne icin kullanildi |
|---|---|---|---|---|
| 33057194 | 2020 | Kaplanis J | *Nature* | birincil kohort, 31.058 DD triosu, ek tablolar |
| 41024706 | 2026 | Zhu Z | *Nucleic Acids Res* | Gene4Denovo2 veritabani |
| 31642496 | 2020 | Zhao G | *Nucleic Acids Res* | Gene4Denovo (ilk surum) |
| 27907889 | 2017 | Turner TN | *Nucleic Acids Res* | denovo-db v1.6.1 |
| 28135719 | 2017 | DDD Study | *Nature* | denovo-db'deki DDD_2017 satirlarinin kaynagi |
| 28263302 | 2017 | Yuen RKC | *Nat Neurosci* | R173Q'yu bildiren MSSNG calismasi |
| 37506195 | 2023 | Cirnigliaro M | *PNAS* | R173Q'yu ikinci kez bildiren AGRE calismasi |
| 23708140 | 2013 | Chesi A | *Nat Neurosci* | Q237H'nin geldigi **ALS** trio calismasi |
| 39954678 | 2025 | Koko M | *Am J Hum Genet* | kontroldeki cerceve kaymasi ve uc kuyruk missense'i |
| 31452935 | 2019 | Feliciano P | *NPJ Genom Med* | A332A'yi bildiren SPARK calismasi |
| 35982159 | 2022 | Zhou X | *Nat Genet* | SPARK/ASC, S1233L ve A332A tekrari |
| 35982160 | 2022 | Fu JM | *Nat Genet* | ASC ekzom calismasi (Gene4Denovo2 uzerinden kapsandi) |
| 38105687 | 2023 | Lin J | *Zhejiang Univ Med Sci* | KIF1A splice de novo'sunun kaynagi |
| 33526817 | 2021 | Lam CW | *J Hum Genet* | KIF13A mozaik cerceve kaymasi (c.2357dupA) |

**Dikkat:** denovo-db calismayi "Chesi2015" diye etiketliyor, ama PMID 23708140 **2013** tarihli
*Nat Neurosci* makalesidir ("Exome sequencing to identify de novo mutations in sporadic ALS trios").
Yil olarak denovo-db'nin etiketi degil, canli dogrulanan 2013 gecerlidir.

---

## 9. Yeniden uretim

Python ortami notu: projenin `mcp\atlas_server\.venv\Scripts\python.exe` yorumlayicisi bu makinede
**calismiyor**; Windows Smart App Control (Uygulama Denetimi) uv'nin yonettigi
`%APPDATA%\uv\python\cpython-3.12.14-windows-x86_64-none\python.exe` dosyasini blokluyor
("Uygulama Denetimi ilkesi bu dosyayi engelledi"). Script bu yuzden sistem Python 3.13'u ile
calistirildi; `atlas_mcp.core` bu yorumlayicida sorunsuz ice aktariliyor (httpx ve pyyaml mevcut)
ve butun provenans yazimi normal calisiyor.

```powershell
# calisan yol (16 Eyl 2026'da bu sekilde uretildi)
python -u scripts\denovo_axis.py

# her seyi yeniden indirmek icin
python -u scripts\denovo_axis.py --refresh

# venv onarildiginda tercih edilen yol
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\denovo_axis.py
```

Tek tek dogrulama komutlari:

```powershell
# Kaplanis PMID ve trio sayisi
curl.exe -s "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:33057194&format=json&resultType=core"

# Kaplanis ek dosyalari
curl.exe -sL "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7116826/supplementaryFiles" -o kaplanis.zip

# Gene4Denovo2, KIF13A butun de novo varyantlari
curl.exe -s -X POST "https://java.genemed.tech/admin-api/backend/data/hg38/gene4denovo2_all/geneSymbol/page" `
  -H "Content-Type: application/json" -H "Origin: https://genemed.tech" `
  -d '{\"reference\":\"hg38\",\"queryby\":\"geneSymbol\",\"terms\":\"KIF13A\",\"page\":1,\"pageNo\":1,\"currentPage\":1,\"pageSize\":2000,\"filterParam\":\"[]\",\"mode\":\"variant\"}'

# denovo-db (SSC olmayan kisim)
curl.exe -sL "https://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz" -o denovodb.tsv.gz

# mutasyon hizlari
curl.exe -sL "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" -o rates.txt.bgz

# 253 kb delesyonun kalitim alani
curl.exe -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=442594&retmode=xml"
```

---

## 10. Bu eksende ogrenilen tuzaklar

16. **Kaplanis 2020'nin ek tablosu GRCh37'dir ve bunu hicbir yerde yazmaz.** Koordinatlari GRCh38
    sanip VEP'e verirseniz dordunun de sonuc sinifi degisir (missense'ler intron olur). Varsaymak
    yerine test edin: her satiri iki yapida da anotlayip tablonun kendi `consequence` sutunuyla
    hangisinin tuttugunu sayin (4/4 GRCh37).
17. **Gene4Denovo2'nin iki anotasyon sutunu farkli izoformlardadir.** `AAChange_refGene` KIF13A icin
    NM_001105567/8 kullanir ve MANE'den **13 kalinti dusuk** numara verir. Dogrusu
    `AAChange_mrefGene` (NM_022113). Bu, 1770 tuzagindan **farkli, ucuncu** bir izoform kaymasidir.
18. **denovo-db her varyant icin her transkriptten bir satir yazar.** KIF13A'nin "101 satiri"
    aslinda 25 varyanttir. Satir saymak, sayiyi dort katina kadar sisirir.
19. **Gene4Denovo2'nin LoF etiketi tek basina guvenilmez.** KIF1A'da 3 bazlik bir cerceve ici
    insersiyon `stopgain` olarak etiketlenmis; Kaplanis'in kendi tablosunda ayni varyant
    `inframe_insertion`. LoF sayilari daima orijinal calismanin tablosuyla capraz kontrol edilmeli.
20. **`1 - CDF` formundaki Poisson p degeri pozitif kontrollerde tasar** ve sahte bir `p = 0`
    basar (ANKRD11: 111 gozlem, 0,31 beklenti). Ust kuyrugu log uzayinda toplayin.
21. **VCV000442594'un `Origin` alani "unknown"dur.** 253 kb'lik delesyon projede sik sik doz kaniti
    olarak aniliyor, ama **de novo oldugu kayitli degil**; de novo argumanina giremez.
22. **Sifir gozlem, sifir etki degildir.** KIF13A icin beklenen de novo LoF zaten 0,37 (Kaplanis) ve
    1,16 (toplam). Bu kohortlar yalniz 8,2x ve 2,6x uzerini disler. Sifiri "kanitlanmis yokluk"
    diye yazmak, bu belgedeki en kolay hatadir.
