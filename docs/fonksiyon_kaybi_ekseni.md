# Fonksiyon kaybi (haployetmezlik) ekseni

**Olcum tarihi: 16 Eylul 2026.** Uretici: `scripts/lof_axis.py`. Tablolar: `data/interim/lof/`.
Sekil: `figures/KIF13A_lof_axis.svg`. Butun kisitlama degerleri 16 Eylul 2026'da yeniden cekildi ve
15 Eylul cekimiyle ayni cikti. Dogrulayamadigim her sey acikca **DOGRULANMADI** olarak isaretli.

**Sorulan soru:** Proje motor alan missense varyantlari uzerine kuruldu ve modelini KIF1A'dan aldi.
Buna karsilik insanda gozlenen her sey fonksiyon kaybi (LoF) sekilli. Bu belge, haployetmezlik
hipotezini acikca **curutebilecek** sekilde test eder.

---

## 0. Tek paragrafta cevap

KIF13A **gercekten LoF intoleransi gosteriyor** (pLI 0,9998, LOEUF 0,449, gozlenen/beklenen LoF 0,372,
DECIPHER pHaplo 0,990), ve bu intolerans **paralogu KIF13B'de yok** (LOEUF 0,663, pLI 1,4e-20).
Insandaki uc gozlemin **ikisi** (Lam 2021 mozaik frameshift'i ve yalniz KIF13A'yi silen 253 kb'lik
delesyon) LoF sekilli; ucuncusu (kuyruk missense'i p.Lys1122Glu) degil.
**Ama haployetmezlik hipotezinin guclu bicimi
desteklenmiyor:** ClinGen'in haployetmezligi kanitlanmis 12 norogelisimsel genine gore KIF13A
**hepsinden daha az kisitli** (LOEUF 0,449 vs 0,094-0,369); gnomAD v4.1'de **158 tane NMD'ye duyarli,
yuksek guvenli LoF alleli** hala populasyonda dolasiyor (yaklasik her 2.570 kisiden birinde bir tane)
ve **bunlarin 29'u dogrudan motor alanda**; ClinVar'da siniflandirilmis tek bir truncating KIF13A
varyanti **yok**; farede heterozigotlar **hic fenotiplenmedi**, sadece canlilik bakildi.
Durust hukum: **LoF ekseni missense ekseninden daha iyi destekleniyor, ama kanit "haployetmezlik
hastalik yapar" demeye degil, "eksik penetransli veya hafif etkili bir doz duyarliligi olabilir"
demeye yeter.**

---

## 1. Yontem

### 1.1 Veri kaynaklari (hepsi canli, kayitsiz, ucretsiz)

| Kaynak | Ne icin | Uc nokta |
|---|---|---|
| gnomAD v4.1 GraphQL | LoF varyantlari (MANE transkript sorgusu), kisitlama | `gnomad.broadinstitute.org/api` |
| Ensembl REST | MANE ekzon yapisi, CNV bolgelerindeki gen icerigi | `rest.ensembl.org`, `grch37.rest.ensembl.org` |
| ClinGen doz duyarliligi kurasyon listesi | haployetmezlik skoru (HI) | `ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv` |
| NCBI E-utilities (ClinVar) | CNV kayitlari, VCV XML fenotipleri | `eutils.ncbi.nlm.nih.gov` |
| DECIPHER gen sayfasi (sunucu tarafli yuk) | pHaplo, pTriplo, s_het, mekanizma tahmini | `deciphergenomics.org/gene/KIF13A/...` |
| Europe PMC | yayinlanmis insan vakalari | `ebi.ac.uk/europepmc/webservices/rest` |
| IMPC Solr | fare nakavt fenotipleri | `ebi.ac.uk/mi/impc/solr/...` |

### 1.2 "Guvenilir tasiyici" tanimi

Bir LoF satiri **guvenilir** sayilir ise: `joint_ac >= 1` **ve** joint filtreleri arasinda `AC0` yok.
gnomAD, sifir allel sayimli ama dosyada duran varyantlari da dondurur; AC0 filtresi tam olarak
bunlari isaretler. KIF13A'nin 382 LoF satirinin 178'i bu nedenle disari duser.

> **Tuzak (kaydedilmeye deger):** projenin mevcut cekimi (`variants.gnomad_variants`)
> `lof_filter` ve `lof_flags` alanlarini **atiyor**. Oysa gercek bir haployetmezlik allelini
> 3' ucta kesilmis bir yalanci pozitiften ayiran tam da o iki alandir. `scripts/lof_axis.py`
> kendi sorgusunu yazar ve bu alanlari saklar.

### 1.3 NMD kacisi nasil hesaplandi

Nonsense-mediated decay (NMD), erken durdurucu kodonlari tasiyan mRNA'lari yikar. Kodon **son
ekzondaysa** veya **son ekzon-ekzon birlesmesinin 50 nukleotid oncesinden sonraysa** mRNA yikilmaz;
kisalmis protein uretilir. Boyle bir varyant **haployetmezlik kaniti degildir**.

Kural, tahmin edilmedi, **MANE ekzon yapisindan hesaplandi**: Ensembl REST `lookup/id?expand=1;utr=1`
ile ekzonlar ve translasyon sinirlari alindi, her ekzonun kodlayan uzunlugu CDS koordinatina
cevrildi, son ekzon-ekzon birlesmesi CDS koordinatinda bulundu. Varyantin CDS konumu HGVS `c.`
dizisinden okundu (`c.1234+1G>A` icin 1234, `c.1235-1del` icin 1235), boylece splice varyantlari da
konumlanabiliyor.

| Gen | MANE | ekzon | CDS (nt) | son birlesme (CDS) | NMD kacisi baslangici |
|---|---|---|---|---|---|
| KIF13A | ENST00000259711.11 | 39 | 5.418 | 4.581 | CDS 4.532 = **kodon 1511** |
| KIF1A | ENST00000498729.9 | 49 | 5.376 | 5.333 | CDS 5.284 = kodon 1762 |
| KIF13B | ENST00000524189.6 | 40 | 5.481 | 5.218 | CDS 5.169 = kodon 1723 |

KIF13A'nin son ekzonu **837 nt kodlar**: kodon 1528'den stop kodonuna kadar, yani proteinin son
278 kalintisi (**%15,5'i**). Karsilastirma: KIF1A'nin son ekzonu 43 nt, KIF13B'ninki 263 nt kodluyor.
Bu, asagida gorulecegi gibi tabloyu tek basina carpitiyor.

### 1.4 "Gercek LoF adayi" tanimi

Haployetmezlik argumaninin dayandigi sinif su uc kosulu birden saglayan varyantlardir:

1. guvenilir tasiyici (AC >= 1, AC0 yok),
2. LOFTEE `HC` (yuksek guvenli) **ve** hicbir LOFTEE filtresi yok (END_TRUNC, GC_TO_GT_DONOR, ...),
3. 50-nt kuralina gore **NMD'ye duyarli**.

Bu belgede buna **NMD-duyarli HC** deniyor.

### 1.5 Iki bagimsiz NMD olcutunun karsilastirilmasi (yontem denetimi)

Kendi 50-nt kuralimi LOFTEE'nin `END_TRUNC` filtresiyle capraz kontrol ettim:

| Gen | iki olcut de kacis diyor | sadece benim kuralim | sadece END_TRUNC |
|---|---|---|---|
| KIF13A | 41 | 4 | 0 |
| KIF1A | 0 | 1 | 0 |
| KIF13B | 0 | 18 | 0 |

KIF13A'daki 4 uyusmazlik **son intronun splice bolgesindeki varyantlar** (`c.4581+1G>A/T/C`,
`c.4582-1G>A`); LOFTEE bunlara END_TRUNC uygulamiyor ama son birlesmenin bozulmasi da NMD kacisi
uretir, dolayisiyla benim kuralim daha muhafazakar. KIF13B'de LOFTEE **hic** END_TRUNC uygulamiyor
(son ekzon kisa ve korunumlu). *LOFTEE'nin END_TRUNC'u NMD kacisini korunum olcutuyle birlestirdigi
icin bu farkin ortaya ciktigi yorumu bir **CIKARIMDIR**, LOFTEE belgesinden alinti degildir:
DOGRULANMADI.* Her iki olcut de tablolarda ayri sutun olarak duruyor
(`nmd_escape_rule`, `lof_filter`), okuyucu istedigini secebilir.

---

## 2. gnomAD v4.1 LoF dokumu

Tablolar: `data/interim/lof/gnomad_lof_variants_20260916.tsv` (1.223 satir, uc gen),
`data/interim/lof/gnomad_lof_summary_20260916.tsv`.

### 2.1 Ana karsilastirma

| Olcut | **KIF13A** | KIF1A | KIF13B |
|---|---|---|---|
| protein (aa) / ekzon | 1805 / 39 | 1791 / 49 | 1826 / 40 |
| LoF satiri (toplam) | 382 | 335 | 506 |
| **guvenilir tasiyicili** | **204** | 107 | 301 |
| AC0 veya sifir allel | 178 | 228 | 205 |
| frameshift | 89 | 41 | 146 |
| stop_gained | 69 | 30 | 82 |
| splice_donor | 26 | 17 | 41 |
| splice_acceptor | 20 | 19 | 32 |
| LOFTEE HC / LC | 162 / 42 | 103 / 4 | 301 / 0 |
| LOFTEE filtresi END_TRUNC | **41** | 0 | 0 |
| son ekzonda | **38** | 0 | 14 |
| son 50 kodonda | 8 | 1 | 9 |
| **NMD kacisi (50-nt kurali)** | **45** | 1 | 18 |
| **NMD-duyarli HC ("gercek LoF adayi")** | **158** | 102 | 283 |
| guvenilir LoF allel toplami | 765 | 1.188 | 902 |
| bunun NMD kacisinda olan kismi | **456 (%59,6)** | 1 (%0,1) | 49 (%5,4) |
| NMD-duyarli HC allel toplami | 308 | 1.173 | 853 |
| NMD-duyarli HC en yuksek AC | **29** | 824 | 225 |
| singleton (AC = 1) | 142 | 78 | 215 |
| homozigot | **1** | 0 | 0 |
| kumulatif NMD-duyarli LoF AF | 1,95e-4 | 7,71e-4 | 5,95e-4 |
| yaklasik tasiyici sikligi | **1 / 2.570** | 1 / 648 | 1 / 840 |
| gnomAD gozlenen/beklenen LoF | **79 / 212,2 = 0,372** | 46 / 227,3 = 0,202 | 127 / 222,0 = 0,572 |
| LOEUF (persentil) | **0,449 (11)** | 0,259 (3) | 0,663 (21) |
| pLI | 0,99975 | 1,0000 | 1,43e-20 |
| missense Z / o-e missense | 2,79 / 0,861 | 6,87 / 0,677 | 1,80 / 0,913 |

> **Kendi sayimim ile gnomAD'in `obs_lof` sayisi neden farkli?** Ben transkript uzerindeki her LoF
> satirini sayiyorum (204 guvenilir). gnomAD'in kisitlama modeli daha dar bir alt kume kullanir
> (yuksek kapsamli bolgeler, dusuk frekans esigi, yalniz LOFTEE HC) ve 79 verir. *Kisitlama
> modelinin tam filtre listesi bu oturumda gnomAD metodoloji belgesinden okunmadi: **DOGRULANMADI**.*
> Iki sayi farkli seyler olctugu icin **birbiriyle karsilastirilmamalidir**; tabloda ikisi de ayri
> satirda duruyor ve genler arasi karsilastirma her satirda kendi icinde yapiliyor.

### 2.2 Okumasi

1. **KIF13A'nin "yuksek LoF sayisi" buyuk olcude bir NMD kacisi yanilsamasi.** Guvenilir LoF
   allellerinin **%59,6'si** (456/765) NMD'den kacan bolgede. Bu KIF1A'da %0,1 (1/1.188),
   KIF13B'de %5,4 (49/902). Sebep yapisal: KIF13A'nin son ekzonu proteinin %15,5'ini kodluyor.
2. **En sik dort "LoF" alleli gercek LoF degil.** p.Trp1800Ter AC 180 (proteinin son 6 kalintisini
   kesiyor, **tek homozigot bu varyantta**), `c.4581+1G>C` AC 128 (son intronun donoru),
   p.Tyr1599Ter AC 37, p.Glu1530Ter AC 30. Ilk gercek aday ancak besinci sirada geliyor:
   **p.Arg1395Ter, AC 29, ekzon 35**.
3. **Gercek LoF alleli suzuldukten sonra bile KIF13A yalniz degil.** 158 NMD-duyarli HC alleli,
   toplam 308 allel. Populasyonda yaklasik **her 2.570 kisiden biri** bunlardan birini tasiyor.
   KIF13B'de bu oran 1/840, KIF1A'da 1/648.
4. **Kritik gozlem: motor alanda 29 tane NMD-duyarli HC LoF var.** Hepsi LOFTEE HC, hepsi
   NMD'ye duyarli, toplam 39 allel, en yuksegi p.Arg268Ter (AC 5). Yani **proteinin motor alanini
   tumuyle kaybettiren heterozigot varyantlar populasyonda mevcut ve bildirilmis bir fenotipleri
   yok.** Bu, projenin motor alan missense adaylarina dair en dogrudan itiraz.

### 2.3 KIF13A LoF varyantlarinin alan dagilimi

| Alan (MANE) | guvenilir LoF | NMD-duyarli HC | allel toplami |
|---|---|---|---|
| motor (5-352) | 29 | 29 | 39 |
| boyun baglayici (353-380) | 5 | 5 | 19 |
| CC1 (381-468) | 10 | 10 | 11 |
| FHA (469-519) | 2 | 2 | 5 |
| kuyruk (520-1805) | 111 | 70 | 491 |
| protein konumu yok (splice) | 46 | 41 | 199 |

### 2.4 Ekzon duzeyinde

`data/interim/lof/KIF13A_exon_lof_20260916.tsv` 39 ekzonun tamamini veriyor. Uc satir yeterince
anlatiyor:

| Ekzon | kodonlar | guvenilir LoF | NMD-duyarli HC | allel toplami |
|---|---|---|---|---|
| 38 | 1493-1527 | 7 | **0** | 138 |
| 39 (son) | 1528-1806 (stop dahil) | 38 | **0** | 318 |
| 1-37 | 1-1492 | 159 | 158 | 309 |

Son iki ekzon, guvenilir LoF allellerinin %60'ini tasiyor ve **hicbiri** haployetmezlik kaniti
degil.

---

## 3. Kisitlama baglami: ClinGen haployetmezlik paneli

Tablo: `data/interim/lof/constraint_panel_20260916.tsv`.
Panel, ClinGen doz duyarliligi kurasyon listesinden (**dosya tarihi: 16 Eylul 2026**, canli indirildi)
**haployetmezlik skoru 3 = "sufficient evidence for dosage pathogenicity"** olan norogelisimsel,
otozomal genler secilerek kuruldu. Skorlar sabit yazilmadi, her calistirmada dosyadan yeniden okunuyor.
ClinGen olcegi: 0 = kanit yok, 1 = az, 2 = bir miktar, 3 = yeterli, 30 = otozomal resesif fenotiple
iliskili gen. **KIF13A ve KIF13B bu dosyada hic gecmiyor**, yani hic degerlendirilmemisler.

| Gen | ClinGen HI | LOEUF | persentil | pLI | gozl/bekl LoF | missense Z | o-e missense |
|---|---|---|---|---|---|---|---|
| TCF4 | 3 | 0,094 | 0 | 1,0000 | 3 / 82,1 | 4,89 | 0,607 |
| STXBP1 | 3 | 0,099 | 0 | 1,0000 | 3 / 78,5 | 6,67 | 0,442 |
| SCN1A | 3 | 0,107 | 0 | 1,0000 | 13 / 193,7 | 8,78 | 0,584 |
| ADNP | 3 | 0,108 | 0 | 1,0000 | 4 / 84,9 | 2,27 | 0,854 |
| SYNGAP1 | 3 | 0,124 | 0 | 1,0000 | 9 / 126,3 | 8,76 | 0,515 |
| ANKRD11 | 3 | 0,130 | 0 | 1,0000 | 17 / 195,7 | **-1,21** | 1,047 |
| NSD1 | 3 | 0,149 | 0 | 1,0000 | 26 / 242,7 | 5,80 | 0,762 |
| CHD8 | 3 | 0,150 | 0 | 1,0000 | 31 / 279,2 | 8,14 | 0,664 |
| POGZ | 3 | 0,189 | 1 | 1,0000 | 17 / 134,6 | 4,70 | 0,738 |
| ARID1B | 3 | 0,192 | 1 | 1,0000 | 30 / 212,3 | 2,36 | 0,902 |
| DYRK1A | 3 | 0,197 | 1 | 1,0000 | 9 / 79,7 | 4,61 | 0,651 |
| *KIF1A* | *0 (kanit yok, 2013-08-08)* | *0,259* | *3* | *1,0000* | *46 / 227,3* | *6,87* | *0,677* |
| SETD5 | 3 | 0,369 | 7 | 1,0000 | 44 / 153,3 | 2,27 | 0,873 |
| **KIF13A** | **kurate edilmemis** | **0,449** | **11** | **0,9998** | **79 / 212,2** | **2,79** | **0,861** |
| *KIF13B* | *kurate edilmemis* | *0,663* | *21* | *1,43e-20* | *127 / 222,0* | *1,80* | *0,913* |

### 3.1 Bu tablodan cikan uc sey

1. **KIF13A, panelin 12 geninin hepsinden daha az LoF-kisitli.** En yakin komsusu SETD5 (0,369) ve
   aradaki fark kucuk degil: gozlenen/beklenen LoF KIF13A'da 0,372, panelde 0,037-0,287.
   **LOEUF 0,449 "kisitli" demektir, "haployetmezlik geni" demek degildir.**
2. **Missense kisitlamasi haployetmezlik geni olmanin kosulu degil.** ANKRD11 missense Z = -1,21
   (yani missense'e hic kisitli degil) ve yine de ClinGen'in haployetmezlik skoru 3. ARID1B 2,36,
   ADNP 2,27, SETD5 2,27, hepsi KIF13A'nin 2,79'undan dusuk. Yani **KIF13A'nin mutevazi missense
   kisitlamasi LoF mekanizmasina karsi bir kanit degildir**; tersine, haployetmezlik genlerinin
   tipik profiline uyar.
3. **KIF1A bu is icin uygun bir karsilastirma degil.** ClinGen KIF1A icin haployetmezlik skoru
   **0, "No evidence available"** veriyor (son degerlendirme 2013-08-08). KAND missense baskinli bir
   hastalik; KIF1A'nin truncating varyantlari ise `docs/KIF1A_ClinVar_kurasyonu.md`'de resesif
   HSN2C referansi olarak ayri tutulmustu. Dolayisiyla **KIF1A, "haployetmezligin hastalik yaptigi
   gen" olcegi degildir**; o olcegi yukaridaki 12 genlik panel verir. KIF1A karsilastirmasi yine de
   degerli, ama LoF eksenini degil, **truncating varyantlarin heterozigot halde sorunsuz
   tasinabildigini** gosterdigi icin degerli.

> **KIF1A tuzagi:** KIF1A'nin en sik "gercek LoF" alleli `c.4257-2A>C`, **AC 824** (AF 5,1e-4).
> gnomAD varyant sorgusuyla kontrol edildi (16 Eyl 2026): bu varyant KIF1A'nin **17 transkriptinde**
> splice_acceptor, yani transkripte ozgu bir yapay sonuc degil. KIF1A'nin yuksek kumulatif LoF
> frekansi (1/648) buyuk olcude bu tek alleleden geliyor; onu cikarinca 2,60e-4 (yaklasik 1/1.920)
> kaliyor. Bu deger bile KIF13A'nin 1,95e-4'unun ustunde.

---

## 4. Insandaki fonksiyon kaybi gozlemleri

### 4.1 ClinVar dizi varyantlari: siniflandirilmis truncating varyant yok

`KIF13A[gene] AND single_gene[prop]` -> 293 kayit, patojenik/muhtemel patojenik **0**
(bkz. `docs/vaka_arama.md` bolum 4e). Bu 293 kaydin LoF sonuclu olanlari:

| VCV | Degisiklik | Sonuc | Siniflandirma | Yildiz |
|---|---|---|---|---|
| VCV004323838 | c.5C>A (p.Ser2Ter) | nonsense | **yok** | 0 |
| VCV004323825 | c.160-1G>T | splice acceptor | **yok** | 0 |
| VCV004323812 | c.1266+1G>A | splice donor | **yok** | 0 |
| VCV004323809 | c.1741G>T (p.Glu581Ter) | nonsense | **yok** | 0 |
| VCV004323806 | c.1787-1G>T | splice acceptor | **yok** | 0 |
| VCV004289694 | c.2942+4del | splice donor | **yok** | 0 |
| VCV003114419 | c.5375T>C (p.Val1792Ala) | missense + splice donor | Uncertain significance | 1 |

Siniflandirmasi olmayan **alti kaydin tamami** VCV XML'inden tek tek dogrulandi (efetch,
16 Eyl 2026): gonderici **Dr. Peter K. Rogan Lab, Western University**, yontem **in vitro**,
siniflandirma `evidence_only`, gozlem aciklamasi TCGA tumor orneklerinde (sirasiyla TCGA-DJ-A2PW,
TCGA-DJ-A1QM, TCGA-CM-5861, TCGA-ZF-A9R9, TCGA-35-4123, TCGA-DD-A39Z) Veridical ile yapilmis
hesaplamali splicing anotasyonu. **Hicbiri hasta gozlemi degildir**
(`docs/vaka_arama.md` 4b ayni gondericinin baska iki varyanti icin ayni tuzagi kaydetmisti).

**Sonuc: ClinVar'da bugune kadar siniflandirilmis tek bir KIF13A truncating varyanti yatirilmamis.**
Lam 2021'in `c.2357dupA` varyanti da ClinVar'da **yok** (293 kaydin hicbiri frameshift degil).

### 4.2 ClinVar CNV'leri ve icerdikleri protein kodlayan genler

Canli E-utilities sorgusu (16 Eyl 2026): `KIF13A[gene] NOT "single nucleotide variant"[Type of
variation]` -> **24 kayit**; 1 kb'den buyuk olan **20 CNV**. Her birinin protein kodlayan gen icerigi
Ensembl REST `overlap/region` ile (5 Mb'lik pencerelerle) sayildi.
Tablo: `data/interim/lof/KIF13A_cnv_gene_content_20260916.tsv`.

| VCV | Buyukluk | Tip | Siniflandirma | Protein kodlayan gen |
|---|---|---|---|---|
| **VCV000442594** | **253 kb** | **kayip** | **Uncertain significance** | **1: KIF13A** |
| VCV000057174 | 569 kb | kazanc | Uncertain significance | 4: CAP2, FAM8A1, KIF13A, NUP153 |
| VCV000253699 | 1,29 Mb | kazanc | Uncertain significance | 7 |
| VCV000152197 | 3,62 Mb | kazanc | Uncertain significance | 17 |
| VCV000689105 | 4,84 Mb | kazanc | Uncertain significance | 22 |
| VCV000688193 | 8,09 Mb | kayip | **Pathogenic** | 23 |
| VCV000442457 | 10,53 Mb | kayip | **Pathogenic** | 29 |
| VCV004755369 | 10,83 Mb | kayip | **Pathogenic** | 35 |
| VCV001809386 | 10,84 Mb | kayip | **Pathogenic** | 36 |
| VCV000253430 | 12,60 Mb | kayip | **Pathogenic** | 59 |
| VCV000058124 | 18,25 Mb | kazanc | **Pathogenic** | 91 |
| VCV001810398 | 20,88 Mb | kazanc | **Pathogenic** | 98 |
| VCV000563140 | 21,80 Mb | kazanc | **Pathogenic** | 102 |
| VCV000687409 | 23,07 Mb | kazanc | **Pathogenic** | 104 |
| VCV000560067 | 23,85 Mb | duplikasyon | **Pathogenic** | 104 |
| VCV000253526 | 23,85 Mb | kazanc | **Pathogenic** | 104 |
| VCV000150044 | 27,43 Mb | kazanc | **Pathogenic** | 237 |
| VCV000155430 | 46,63 Mb | kazanc | **Pathogenic** | sayilmadi (kromozom kolu olcegi) |
| VCV000443496 / 443497 | 170,76 Mb | kazanc | **Pathogenic** | sayilmadi (kromozom olcegi) |

**Bu tablonun tek cumlelik ozeti:** ClinVar'da KIF13A'yi kapsayan **her patojenik CNV en az 23 protein
kodlayan gen iceriyor**; en kucugu 8,09 Mb. Patojeniklik atfi KIF13A'ya degil, delesyona/duplikasyona
aittir. **Patojenik etiketli hicbir CNV KIF13A'yi tek basina silmiyor.**

### 4.3 Tek gercek doz kaniti: VCV000442594

VCV XML'inden (efetch, 16 Eyl 2026) **kendim dogruladim**, esummary yalniz "See cases" diyor:

| Alan | Deger (VCV XML'den birebir) |
|---|---|
| Varyant | `GRCh37/hg19 6p22.3(chr6:17758616-18011991)x1`, yaklasik 253 kb heterozigot kayip |
| Fenotip | **Delayed speech and language development** |
| Etkilenme durumu | **present** |
| Yontem | clinical testing |
| Kalitim | Origin: unknown (**de novo mu, ebeveynden mi bilinmiyor**) |
| Gonderici | ISCA site 1 / International Standards For Cytogenomic Arrays Consortium |
| Siniflandirma | Uncertain significance, **no assertion criteria provided** |
| Atifli PMID | 20466091 |
| Gen icerigi | Ensembl GRCh37 overlap: **tek protein kodlayan gen KIF13A** |

> `docs/vaka_arama.md` bu kayit icin PMID 21844811 (Kaminsky 2011) kaydediyor. Benim efetch'imde
> kayitta gecen tek PubMed kimligi **20466091**. Iki kaynak farkli alanlari okumus olabilir;
> ikisinin ayni kayda ait olup olmadigi **DOGRULANMADI**.

### 4.4 Ikinci fenotipli kayit (LoF degil, ama LoF tartismasinin bir parcasi)

VCV002570672, `NM_022113.6(KIF13A):c.3364A>G (p.Lys1122Glu)`. VCV XML'den dogrulandi:
fenotipler **Apraxia; Chiari malformation; Moderate intellectual disability; Seizure; Speech
apraxia**; etkilenme durumu **present**; yontem **research**; gonderici HudsonAlpha Institute for
Biotechnology; atifli PMID 25741868. **Bu bir missense VUS'tur, kuyruk bolgesindedir (1122) ve
haployetmezlik kaniti degildir**; LoF ekseninin degil, "KIF13A norogelisimle ilgili olabilir"
iddiasinin kanitidir.

### 4.5 Yayinlanmis vakalar (Europe PMC, 16 Eyl 2026)

| Sorgu | Hit | Sonuc |
|---|---|---|
| `EXT_ID:33526817` | 1 | Lam 2021 dogrulandi |
| `KIF13A AND (haploinsufficiency OR "loss of function" OR frameshift OR nonsense)` | 144 | yeni insan vakasi yok |
| `KIF13A AND ("6p22.3" OR microdeletion) AND (patient OR proband OR case)` | 25 | yalniz PMID 33526817 ve 34303382 |
| `"Kif13a" AND (knockout OR "null mice" OR heterozygous)` | 195 | insan vakasi yok |

**Lam 2021 (PMID 33526817), ozetinden birebir dogrulandi** (*J Hum Genet* 66(8), 2021,
doi 10.1038/s10038-020-00883-w; Europe PMC'de **tam metin yok, abonelik gerekiyor**, dolayisiyla
ozet disindaki ayrintilar **DOGRULANMADI**):

> "a heterozygous de novo frameshift variant in the KIF13A gene, i.e., NM_022113.6: c.2357dupA,
> resulting in nonsense-mediated decay. The low mutant allelic ratio suggested that the mutation
> has occurred postzygotically leading to embryonic mosaicism."

ve ayni ozet acikca "supporting the association of KIF13A haploinsufficiency with the novel
neuroectodermal syndrome" diyor. Hasta: 3 yasinda kiz, deri hipopigmentasyonu (hypomelanosis of Ito)
ve gelisim geriligi.

**Bagimsiz capraz kontrol (kendi ekzon haritamla):** `c.2357dupA` CDS 2357'de, **ekzon 20/39**,
CDS aralik 2305-2454, kodon 786. NMD kacis sinirim CDS 4532. Yani **benim kuralim da bu varyanti
NMD'ye duyarli sayiyor** ve makalenin "resulting in nonsense-mediated decay" ifadesiyle uyusuyor.
Bu, bugune kadarki **tek NMD-duyarli, fenotipli, de novo KIF13A truncating varyanti**.

Ikinci yayin PMID 34303382 (*Mol Cytogenet* 2021): 2,5 Mb de novo 6p22.3 delesyonu, **13 gen**,
gelisim geriligi. Bitisik-gen sendromu; KIF13A'ya dair tek basina kanit degil.

### 4.6 DECIPHER: gen duzeyi doz metrikleri (hesapsiz okunabildi)

DECIPHER'in JSON API'si anahtar istiyor (HTTP 403), ama **gen sayfasinin sunucu tarafli yuku
(SSR payload) hesapsiz okunabiliyor** ve gen duzeyi metrikleri tasiyor.
`data/raw/decipher/KIF13A_gene_ssr_20260916.json`:

| Metrik | Deger | Okumasi |
|---|---|---|
| `p_hi` (pHaplo) | **0,990** | haployetmezlik olasiligi; yuksek |
| `p_ts` (pTriplo) | 0,604 | triplo duyarlilik icin belirgin degil |
| `s_het` | 0,0489 | heterozigot secilim katsayisi, orta duzey |
| `hi_index` (HI index) | %42,3 | **haployetmezlige egilimli degil** (dusuk persentil beklenirdi) |
| DECIPHER LOEUF / pLI | 0,471 / 0,992 | bizim 0,449 / 0,9998 ile yakin, ayni degil (farkli surum) |
| Mekanizma tahmini | pLoF 0,616; pDN 0,560; pGOF 0,598; **rank_max = "lof"** | uc mekanizma arasinda **LoF basta** |
| ClinGen doz kaydi | **yok** | bolum 3 ile ayni sonuc |

**Iki metrik celisiyor:** `p_hi` 0,990 haployetmezligi destekliyor, `hi_index` %42,3 desteklemiyor.
Bunu saklamak yerine yan yana raporluyorum.

> Alan adlari DECIPHER'in kendi SSR yukundeki adlardir. `p_hi`/`p_ts`'nin Collins ve ark. 2022
> pHaplo/pTriplo skorlari, `hi_index`'in Huang ve ark. haployetmezlik indeksi, mekanizma
> skorlarinin Badonyi ve Marsh tahmin edicisi oldugu **DECIPHER arayuzundeki adlandirmadan
> cikarilan bir yorumdur; kaynak yayin eslemesi bu oturumda DOGRULANMADI.** Sayilar olculdu,
> etiketleri degil.

> **DECIPHER hasta CNV'leri:** `docs/vaka_arama.md` (16 Eyl 2026) gen sayfasindan okunan
> **55 eslesen hasta varyanti, bunlarin 47'si kopya sayisi degisimi** rakamini kaydediyor.
> Ben bunu **tabloya cekemedim**: `/api/` uc noktalarinin tamami HTTP 403
> (`"Client key not provided"`) donuyor ve SSR yuku `public_patients: 0` diyor (bu alan acik erisimli
> **dizi** varyantlarini sayiyor). **47 CNV rakami bu belgede bagimsiz olarak DOGRULANMADI**;
> tarayicidan okunan bir kayit olarak kalir ve her birinin gen icerigi bilinmedigi icin
> **doz kaniti olarak kullanilamaz.**

---

## 5. Fare kaniti

Kaynak: IMPC Solr, canli sorgu 16 Eyl 2026. `data/raw/impc/KIF13A_impc_20260916.json`.

### 5.1 Ozet

| | Kif13a | Kif13b | Kif1a |
|---|---|---|---|
| MGI | MGI:1098264 | - | - |
| Istatistiksel olarak cagrilan fenotip | **21** | 4 | **0** |
| Bunlarin zigotlugu | **hepsi homozigot** | hepsi homozigot | - |
| Toplam istatistik sonucu | 500 (498 hom, 1 het, 1 hemi) | - | - |
| **Heterozigot testi** | **1** | 45 | 0 |
| Heterozigotta anlamli bulgu | **0** | 0 | - |

**Kif1a, IMPC tarafindan hic fenotiplenmemis** (0 sonuc). Bu eksende KIF1A karsilastirmasi yapilamaz.

### 5.2 Homozigot Kif13a nakavtinin fenotipi

Allel: `Kif13a<em1(IMPC)Mbp>` (Mouse Biology Program, UC Davis), tum cagrilar ayni merkezden.
21 anlamli cagrinin norolojik ve duyusal olanlari (tamami JSON dosyasinda):

| Fenotip (MP) | Parametre | p | Etki | Cinsiyet |
|---|---|---|---|---|
| decreased prepulse inhibition | % PPI2 | **5,19e-08** | -30,0 | disi |
| decreased startle reflex | Response amplitude - S | **1,41e-07** | -487,0 | erkek |
| hyperactivity | Center distance travelled | **3,84e-07** | +1028,9 | erkek |
| decreased prepulse inhibition | % PPI Global | 1,07e-06 | -23,9 | disi |
| decreased prepulse inhibition | % PPI1 | 1,68e-06 | -32,6 | disi |
| abnormal auditory brainstem response | Click-evoked ABR | 2,88e-06 | +0,97 | ikisi |
| hyperactivity | Number of center entries | 1,04e-05 | +49,3 | erkek |
| cataract | Lens Opacity | 1,03e-05 | +0,59 | ikisi |
| abnormal vitreous body morphology | Vitreous | 7,66e-06 | +0,59 | ikisi |

Ayrica homozigotlarda bobrek (kucuk bobrek, anormal morfoloji), testis, uterus (hydrometra) ve deri
morfolojisi anormallikleri (p = 0,0, Gross Pathology). IMPC'nin ust duzey terim listesi:
**nervous system, behavior/neurological, vision/eye, hearing/vestibular/ear, integument,
reproductive, renal/urinary, endocrine/exocrine** anlamli; **mortality/aging anlamli degil**.

### 5.3 Canlilik

Viability Primary Screen sayimlari (IMPC experiment core): **170 yavru; 42 WT, 92 heterozigot,
36 homozigot** (22 erkek + 14 disi). Homozigot orani %21,2, Mendel beklentisi %25.
Test edilen dort canlilik parametresinin dordunde de p = 1,0, anlamli degil.
**Homozigot Kif13a nakavti yasiyor.**

### 5.4 En onemli nokta: heterozigot hic bakilmamis

Insanda haployetmezlik, farede **heterozigot** fenotip bekletir. Kif13a icin IMPC'de heterozigot
uzerinde yapilmis **tek bir istatistiksel test var** ve o da canlilik
("Heterozygous animals viability", p = 1,0, anlamli degil). PPI, startle, open field, ABR gibi
homozigotta anlamli cikan hicbir bataryada heterozigot **olculmemis**.

Bunun "kanit yoklugu" oldugunu, "yokluk kaniti" olmadigini vurgulamak gerekiyor. Kontrol:
**Kif13b'de 45 heterozigot testi var** (hicbiri anlamli degil), yani IMPC bazi hatlarda
heterozigotlari gercekten fenotipliyor; Kif13a hattinda yapmamis.

### 5.5 Allen Brain Cell Atlas ile capraz kontrol

`C:\Users\pc\Claude\kif13a-allenbrain` (README ve `results/derived_*.csv`, TRUBA uzerinde uretildi):

- Kif13a, 2.341.350 hucre ve 265 yeterince orneklenmis alt sinif icinde **olgun oligodendrositlerde
  birinci** (ortalama log2(CPM+1) = 7,77, ifade eden hucre orani 0,93), **OPC'lerde ikinci** (6,76 / 0,85).
  Ucuncusu mikroglia (3,88).
- Yaslanmayla **yalniz olgun oligodendrositlerde** artiyor: 53 erisin ve 44 yasli fare, fare basina
  tek deger, delta = +0,250 log2 (1,19 kat), Mann-Whitney p = 2,42e-04, Cliff's delta = +0,435.
  OPC'de anlamli degil (p = 0,219).
- Ayni repo'nun kendi uyarisi: siralama tespit oranini ifade duzeyiyle harmanliyor; "oligodendrosit
  spesifiktir" iddiasi bu veriyle **kanitlanmis degil**.

**LoF ekseniyle iliskisi:** eger mekanizma haployetmezlik ise, tek kopya kaybinin en cok etkileyecegi
hucre tipi bu gene en cok dayanan tiptir. Allen verisi bunun **oligodendrosit soyu** oldugunu
soyluyor. IMPC'nin homozigot fenotipi (PPI dusuklugu, artmis aktivite, anormal ABR) bu ile **tutarli
ama ispatlayici degil**: hicbiri oligodendrosit-spesifik bir okuma degildir ve
**miyelin veya beyaz cevher fenotipi IMPC bataryasinda olculmemis** (test edilen prosedurler listesi
`data/raw/impc/KIF13A_impc_20260916.json` icinde). Bu bir **CIKARIMDIR**, olculmus bir baglanti degil.

---

## 6. Hukum

### 6.1 LoF hipotezi missense hipotezinden daha iyi mi duruyor?

**Evet, ama mutlak anlamda degil, goreli olarak.**

**LoF lehine (olculmus):**

1. Insanda bugune kadarki dogrulanmis gozlemlerin **cogu** LoF veya doz sekilli: Lam 2021'in
   NMD-duyarli de novo mozaik frameshift'i (ekzon 20/39, kodon 786), yalniz KIF13A'yi silen
   253 kb'lik delesyon (etkilenmis cocuk), ve bitisik-gen 6p22.3 delesyonlari. Tek missense
   gozlemi (p.Lys1122Glu) kuyruk bolgesinde ve tek gonderenli bir arastirma VUS'u.
2. Gen duzeyi kisitlama LoF'a isaret ediyor: pLI 0,9998, LOEUF 0,449 (11. persentil),
   gozlenen/beklenen LoF 0,372, DECIPHER pHaplo 0,990, mekanizma tahmininde `rank_max = "lof"`.
3. **Paralog kontrolu temiz:** KIF13B ayni boyda, ayni ailede, ayni bataryadan geciyor ve
   **LoF'a hic kisitli degil** (LOEUF 0,663, pLI 1,4e-20, gozlenen/beklenen 0,572). Yani KIF13A'nin
   LoF kisitlamasi bir kinezin artefakti degil, gene ozgu.
4. Missense tarafinda **hicbir sey yok**: 52 tier A adayimizin hicbiri hicbir yayinda gecmiyor
   (`docs/vaka_arama.md` 5c), ClinVar'da fenotipli tek missense kuyruk bolgesinde (p.Lys1122Glu) ve
   bizim listemizde degil.
5. Missense kisitlamasinin mutevazi olmasi (Z = 2,79) **LoF mekanizmasina karsi bir kanit degil**:
   panelin haployetmezlik genlerinin cogu bundan daha dusuk missense Z tasiyor (ANKRD11 -1,21).

**LoF aleyhine (olculmus):**

1. **158 NMD-duyarli, yuksek guvenli LoF alleli** populasyonda sag kaliyor; yaklasik 1/2.570 tasiyici.
   Bunun ayni bicimde olculmus karsiligi panel genleri icin cikarilmadi, ama gnomAD'in kendi
   modellenmis olcusu ayni yone isaret ediyor: gozlenen/beklenen LoF KIF13A'da 0,372, panelde
   0,037-0,287.
2. **Motor alanda 29 tane NMD-duyarli HC LoF var** ve bildirilmis fenotipleri yok.
3. **LOEUF 0,449, panelin 12 geninin hepsinden yuksek.** En yakin komsu SETD5 (0,369).
4. ClinGen KIF13A'yi **hic kurate etmemis**; ClinVar'da siniflandirilmis tek bir truncating varyant
   yok; patojenik etiketli hicbir CNV KIF13A'yi tek basina silmiyor.
5. Tek temiz insan doz kaniti (VCV000442594) **tek bir cocuk, kriter uygulanmamis VUS, kalitim
   bilinmiyor**. Aile calismasi olmadan bu kayit tasarimi geregi zayif.
6. Farede heterozigotlar **fenotiplenmemis**; homozigotlar yasiyor. Yani en dogrudan model testi
   yapilmamis durumda.
7. DECIPHER'in `hi_index` degeri (%42,3) kendi `p_hi` degeriyle (0,990) **celisiyor**; iki metrik
   ayni gen icin zit sonuc veriyor.

### 6.2 Karar

**Proje icin onerilen konum:** KIF13A'yi "motor alan missense hastaligi adayi" olarak degil,
**"doz duyarli olabilen, ama haployetmezligi henuz kanitlanmamis bir norogelisimsel aday"** olarak
konumlandirmak. Somut olarak:

- Motor alan missense tier A listesi **birakilmamali**, cunku dominant-negatif bir mekanizma
  LoF'tan farkli ve daha siddetli bir etki uretebilir (DECIPHER mekanizma tahmini pDN 0,560,
  dn_pctl %73,3; KIF13A'nin dimerize olup olmadigi bu oturumda **DOGRULANMADI**).
  **Ama tier A listesi artik "haployetmezlik yapar" gerekcesiyle savunulamaz:**
  populasyonda zaten 29 heterozigot motor-alan truncation var. Basit bir "motilite kaybi" missense'i
  en fazla bir truncation kadar zararli olabilir, o da fenotipsiz gorunuyor.
  **Yani tier A listesinin mantikli tek gerekcesi artik dominant-negatif etkidir ve bu henuz
  olculmedi.**
- Vaka arama onceligi `docs/vaka_arama.md`'nin R0 rotasi olarak kalmali (iki fenotipli ClinVar
  kaydinin gondericisi), ama arama sorgusu **truncating ve delesyon varyantlarini da kapsamali**.

### 6.3 Hipotezi ne curutur

Her madde bugun olculebilir veya yakin vadede olculebilir hale gelir:

1. **Kisitlamanin ornek buyudukce dusmemesi.** Haployetmezlik geni olsaydi, ornek buyudukce
   gozlenen/beklenen LoF orani panel genlerinin araligina (0,04-0,29) yaklasmasi beklenirdi.
   gnomAD v5 ve sonrasinda KIF13A hala **0,37 civarinda** kalirsa ve panel genleri 0,10-0,20
   araliginda kalirsa, **hipotezin guclu bicimi curur.**
2. **Saglikli aktaran ebeveynler.** Iki veya daha fazla ailede NMD-duyarli bir KIF13A truncating
   alleli **etkilenmemis bir ebeveynden** aktarilmis olarak bulunursa, haployetmezlik penetransi
   pratikte cokmustur. (Bugun veri yok: 158 alelin hicbirinde segregasyon bilgisi yok.)
3. **Populasyon CNV setlerinde KIF13A-tek delesyonu.** gnomAD-SV veya benzeri bir saglikli-kohort
   yapisal varyant setinde KIF13A'yi tek basina silen bir heterozigot delesyon bulunursa,
   VCV000442594'un tasidigi tum agirlik gider. *(Bu belgede aranmadi: dbVar'da KIF13A'yi ortusen
   945 kayit var, bunlarin cogunlugu populasyon calismalari; ayristirilmasi ayri bir istir.
   DOGRULANMADI.)*
4. **Farede heterozigot bataryasi.** IMPC'nin homozigotta anlamli cikardigi bataryada (PPI, startle,
   open field, ABR) heterozigotlar test edilir ve **hicbiri anlamli cikmazsa**, model tarafinda
   haployetmezlik desteklenmez. Anlamli cikarsa en guclu tek kanit bu olur.
   Bugun bu test **yapilmamis** (tek het testi canlilik).
5. **Uc veya daha fazla bagimsiz de novo truncating proband.** GeneMatcher / DDD / ClinVar
   uzerinden tutarli norogelisimsel fenotipli uc bagimsiz de novo truncating KIF13A varyanti
   bulunursa, ClinGen doz duyarliligi olceginde en azindan "some evidence" (skor 2) esigine ulasilir
   ve hipotez **dogrulanma yoluna girer**. Bugun sayi **bir** (Lam 2021) ve o da mozaik.
6. **Mekanizma testi (missense vs LoF ayrimi).** Bir motor-alan tier A missense varyanti hucrede
   veya in vitro motilitede **truncating bir allelden daha kotu** cikarsa (baskin-negatif),
   LoF-yalniz model yetersiz demektir ve missense ekseni geri doner. Ayni missense en fazla
   truncating kadar kotuyse, missense ekseninin ayri bir hipotez olarak var olma sebebi kalmaz.
7. **Oligodendrosit testi.** Allen verisi Kif13a'yi oligodendrosit soyuna yerlestiriyor. Heterozigot
   farede (veya insan iPSC oligodendrositlerinde) miyelinasyon veya beyaz cevher olcumu yapilir ve
   fark cikmazsa doz duyarliligi hipotezinin en makul hucresel mekanizmasi kaybolur.

---

## 7. Yeniden uretim

```powershell
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\lof_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\lof_axis.py --refresh --no-figure
```

Cikan dosyalar (hepsi `.prov.json` yan dosyasiyla, `data/raw/PROVENANCE.jsonl` icinde kayitli):

| Dosya | Icerik |
|---|---|
| `data/interim/lof/gnomad_lof_variants_20260916.tsv` | uc genin butun LoF varyantlari, ekzon numarasi, NMD verdikti, LOFTEE alanlari |
| `data/interim/lof/gnomad_lof_summary_20260916.tsv` | gen basina ozet (bolum 2.1'in kaynagi) |
| `data/interim/lof/KIF13A_exon_lof_20260916.tsv` | 39 ekzonun kodon araliklari ve LoF yuku |
| `data/interim/lof/constraint_panel_20260916.tsv` | ClinGen HI paneli + gnomAD kisitlama (bolum 3) |
| `data/interim/lof/KIF13A_cnv_gene_content_20260916.tsv` | ClinVar CNV'leri ve protein kodlayan gen icerigi |
| `data/interim/lof/KIF13A_clinvar_lof_sequence_20260916.tsv` | ClinVar LoF sonuclu dizi varyantlari |
| `data/interim/lof/KIF13A_phenotyped_records_20260916.tsv` | VCV XML'den dogrulanmis iki fenotipli kayit |
| `data/interim/lof/lof_axis_summary_20260916.json` | butun eksenin makine okunur ozeti |
| `figures/KIF13A_lof_axis.svg` | uc panel: konum haritasi, LoF bilesimi, LOEUF paneli |

## 8. Kaynaklar

- gnomAD v4.1, GraphQL transkript ve gen sorgulari, 15-16 Eyl 2026.
- Ensembl REST (GRCh38 ve GRCh37), MANE ekzon yapisi ve bolge ortusmesi, 16 Eyl 2026.
- ClinGen Gene Curation Results, `ClinGen_gene_curation_list_GRCh38.tsv`, dosya tarihi 16 Eyl 2026.
- ClinVar, NCBI E-utilities esearch/esummary/efetch, 16 Eyl 2026.
- DECIPHER gen sayfasi sunucu tarafli yuku, 16 Eyl 2026 (JSON API kapali, HTTP 403).
- Europe PMC REST, 16 Eyl 2026.
- IMPC Solr (genotype-phenotype, statistical-result, experiment, gene), 16 Eyl 2026.
- Lam CW, Chan CY, Wong KC, Chang ST. *J Hum Genet* 2021;66(8). **PMID 33526817**,
  doi 10.1038/s10038-020-00883-w. Ozet Europe PMC'de dogrulandi; tam metin abonelik arkasinda.
- **PMID 34303382**, *Mol Cytogenet* 2021, 2,5 Mb de novo 6p22.3 delesyonu, 13 gen.
- Allen Brain Cell Atlas analizi: `C:\Users\pc\Claude\kif13a-allenbrain`, release `20260711`,
  `results/derived_donor_level_age_test.csv` ve `results/derived_wholebrain_subclass_ranking.csv`.
- Ilgili proje belgeleri: `docs/vaka_arama.md` (vaka arama rotalari ve ClinVar profili),
  `docs/KIF1A_ClinVar_kurasyonu.md` (KIF1A etiketleri), `docs/yapi_ekseni.md` (yapi ekseni).
