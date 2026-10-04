# HudsonAlpha KIF13A kaydının adversaryal doğrulaması

**Kayıt:** VCV002570672 / SCV004009653 · `NM_022113.6(KIF13A):c.3364A>G (p.Lys1122Glu)`
**Doğrulama tarihi:** 16 Eylül 2026 · tüm sayılar bu tarihte canlı API'lerden çekildi
**Yöntem:** NCBI E-utilities (ClinVar), gnomAD GraphQL v4.1, Ensembl REST (VEP + Compara), UniProt REST,
MyVariant.info (dbNSFP), Europe PMC REST, EBI G2P, GenCC, ClinGen, Genomics England PanelApp, LOVD REST

---

## Hüküm

**Bu kayıt KIF13A'yı bir nörogelişimsel gen olarak destekleme gücü: "suggestive" (en zayıf basamak).**
Kayıt, seçilmiş bir aday gen çağrısı değil, HudsonAlpha'nın CSER programından gelen 384 varyantlık toplu
bir araştırma yatırımının içindeki 124 VUS'tan biridir; aynı partideki 152 varyant **de novo** olarak
işaretlenmişken bizimkinin alel kökeni **unknown**, etkilenme durumu da **unknown**'dır ve gönderici
hiçbir serbest metin yorumu bırakmamıştır (aynı partide 43 kayıtta yorum var). Varyantın gnomAD v4.1'de
hiç görülmemiş olması gerçektir ama ayırt edici değildir: KIF13A'nın teorik olarak mümkün 11.918 missense
SNV'sinin yaklaşık %82'si gnomAD'de yoktur, yani "yok" olmak bu gen için varsayılan durumdur. Kalıntı
1122 gerçekten korunmuştur (203 ortologun 186'sında Lys) ve bu, kaydın tek güçlü tarafıdır; ancak kalıntı
motor alanda değil, UniProt'un 1100-1138 coiled coil'inde oturur ve o coiled coil içinde gnomAD'de AC>0
olan 44 farklı missense vardır, bunların arasında **aynı Lys'ten Glu'ya değişimi taşıyan p.Lys1104Glu** ve
AC=27 ile taşınan bir yük tersine çevirmesi olan p.Glu1128Lys de bulunur. AlphaMissense 0,8044 (doğrulandı)
bizim gen-kalibre eşiğimiz olan 0,968'in altındadır, MPC 0,573 ve REVEL 0,642 düşüktür, PrimateAI / FATHMM /
DEOGEN2 varyantı "tolere edilir" olarak çağırır. Fenotip listesindeki Chiari malformasyonu gelişimsel
gerilik nedeniyle görüntülenen çocukların yaklaşık %2,3'ünde bulunur ve o seride 6/6 vakada "etiyolojiyle
ilgisiz" olarak sınıflanmıştır. Varyant ya da kalıntı hakkında sıfır yayın vardır; KIF13A ClinGen gen
geçerliliği, GenCC, G2P DD paneli ve Genomics England PanelApp'in **hiçbirinde** yoktur ve ClinVar'da tek
bir patojenik KIF13A nokta varyantı bile yoktur. Dolayısıyla bu kayıt bir hipotez üretir, bir hipotezi
desteklemez.

---

## 1. KIF13A gerçekten aday mıydı, yoksa hastanın asıl tanısı başka yerde miydi?

### 1.1 Kaydın kendisi hangi çalışmadan geldiğini söylüyor

VCV XML'de, kendi SCV'sinin içinde, projenin daha önce okumadığı bir alan var:

```
<StudyName>CSER-HudsonAlpha</StudyName>
<SubmissionNameList><SubmissionName>SUB13646016</SubmissionName></SubmissionNameList>
<ClinVarSubmissionID localKey="6:17785639" localKeyIsSubmitted="1" submittedAssembly="GRCh38"/>
<Attribute Type="SampleLocalID">6:17785639:A00166</Attribute>
```

`localKey` bir gen adı ya da vaka kimliği değil, **genomik pozisyonun kendisidir**. Bu, elle kürate edilmiş
tek bir aday çağrısının değil, pozisyona göre anahtarlanmış bir boru hattı yatırımının imzasıdır.

Çekilen URL: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=2570672&from_esearch=true`

### 1.2 HudsonAlpha'nın ClinVar ayak izi: "tek KIF13A kaydı" en sık görülen desendir

ClinVar gönderici sayfasından (`https://www.ncbi.nlm.nih.gov/clinvar/submitters/505530/`) doğrulandı:
gönderici **1.293** kayıt, PI **Greg Cooper**, irtibat **Susan Hiatt**, Huntsville, Alabama.
E-utilities ile `HudsonAlpha[SBM]` üzerinden 1.276 ayrı varyasyon kaydı ve bunların içindeki
**1.293 HudsonAlpha SCV'si** XML'den ayrıştırıldı. Sayı iddiadaki 1.293 ile birebir tutuyor.

| Ölçüm | Değer |
|---|---|
| HudsonAlpha SCV sayısı | **1.293** |
| Ayrı varyasyon kaydı | 1.276 |
| Ayrı örnek (proband) kimliği | **1.118** |
| Toplama yöntemi "research" olan kayıt | **1.276 / 1.276 (%100)** |
| En az 1 HPO terimi taşıyan SCV | **1.028 / 1.293 (%79,5)** |
| Serbest metin yorumu taşıyan SCV | 496 / 1.293 |
| Tek gene düşen varyasyon kaydı | 1.068 |
| Bu tek gen kayıtlarının kapsadığı ayrı gen | **607** |
| Bunlardan **tam olarak 1** kaydı olan gen | **368 (%60,6)** |
| CNV'ler dahil dokunulan tüm genler | 4.112 |
| Bunlardan tam olarak 1 kaydı olan gen | **2.257** |
| KIF13A kaydı | 1 |

Çalışma adı dağılımı (SCV düzeyinde, 15 parti):
CSER-HudsonAlpha 384 · CSER-SouthSeq 370 · HudsonAlpha-AGHI-WGS 328 · AGHI-WGS-HudsonAlpha 50 ·
AGHI_GT 45 · isimsiz 42 · HudsonAlpha-PGEN 19 · AGHI-GT-HudsonAlpha 13 · HudsonAlpha-COAGS 12 ·
PGEN-HudsonAlpha 8 · UCC-HudsonAlpha 7 · AGHI GT 6 · CSER 5 · CSER-NYCKidSeq 2 · St. Jude-INSIGHT-HD 2

**Sonuç:** "HudsonAlpha'nın KIF13A'da tam olarak 1 kaydı var" cümlesi 368 geni (tek gen kayıtları içinde)
veya 2.257 geni (hepsi içinde) aynı şekilde tarif eder. Bu, verilerindeki **en yaygın** desendir. Projenin
"kasıtlı aday gen çağrısı" okuması bu sayının önünde ayakta kalamaz.

### 1.3 CSER-HudsonAlpha partisinin anatomisi (384 SCV, XML'den ayrıştırıldı)

| Ölçüm | Değer | Bizim kayıt |
|---|---|---|
| SCV sayısı | 384 | 1'i |
| Ayrı örnek kimliği | **358** (proband başına 1,20 varyant) | A00166, 1 varyant |
| Ayrı gen | 297 (tek gen kayıtlarında 228, bunların 136'sı tekil) | KIF13A tekil |
| Sınıflama (SCV düzeyi) | Patojenik 165 · **VUS 124** · LP 76 · LB 18 · B 1 | VUS |
| Alel kökeni | **de novo 152** · unknown 114 · maternal 68 · paternal 47 · biparental 3 | **unknown** |
| Etkilenme durumu | **yes 278** · unknown 94 · no 12 | **unknown** |
| >= 1 HPO terimi | 253 / 384 (ortalama 3,3, en çok 17) | 5 terim (ortalama civarı) |
| Serbest metin yorumu | 43 / 384 | **yok** |
| NumberTested = 1 | 338 / 384 | 1 |

Bu tablo iki yönlü okunmalı. HudsonAlpha **kalıtım bilgisini biliyorsa yazıyor**: aynı partide 152 varyant
açıkça "de novo" işaretli. Bizim varyantımız de novo değil, "unknown". Aynı şekilde 278 kayıtta etkilenme
durumu "yes" iken bizimkinde "unknown". Yani kaydın **kendi partisi içinde**, bilgi taşıyan her eksende
zayıf kuyrukta duruyor.

Sorgular: `"CSER-HudsonAlpha"[STNM]` (384), `... AND "de novo"[ORIG]` (162 varyasyon düzeyi),
`HudsonAlpha[SBM] AND research[COLLMETHOD]` (1.276), `"SouthSeq"[STNM]` (0, doğru etiket CSER-SouthSeq).

### 1.4 Proband A00166'nın başka bir kaydı var mı?

1.293 SCV'nin tamamı SampleLocalID'ye göre tarandı: **A00166 yalnızca bu tek kayıtta geçiyor.** Yani
HudsonAlpha bu proband için ClinVar'a rakip bir aday yatırmamış.

Bu, iddiayı **desteklemez**. ClinVar'a yatırılmamış olmak, alternatif tanı olmadığı anlamına gelmez;
HudsonAlpha yalnızca seçtiklerini yatırır ve bu probanda ait kaydı alel kökeni ve etkilenme durumu
"unknown" bırakarak yatırmıştır, yani **tanıyı kurmamıştır**. CSER kohortlarında probandların
yaklaşık %70'i tanısız kalır ve tanısız probandların aday VUS'ları tam olarak bu biçimde yatırılır.

### 1.5 Literatür: HudsonAlpha VUS'ları seçerek değil, sistematik olarak yatırdığını yazıyor

Europe PMC canlı sorguları (hitCount'lar gerçek yanıtlardan):

| Sorgu | Hit |
|---|---|
| `AUTH:"Cooper GM" AND HudsonAlpha` | 88 |
| `SouthSeq` | 36 |
| `CSER AND "HudsonAlpha"` | 64 |
| `"CSER-HudsonAlpha"` | **2** (ikisi de HudsonAlpha çalışma makalesi değil) |

**"CSER-HudsonAlpha" etiketinin arkasında bir yayın yok.** Bu bir ClinVar gönderim partisi etiketidir.

Kohort makaleleri (canlı doğrulandı):

| PMID | Yıl | Makale | Kohort |
|---|---|---|---|
| 28554332 | 2017 | Bowling ve ark., Genome Med | **371 birey**, 309 proband-ebeveyn üçlüsü; %27 P/LP, **%11,3 VUS** |
| 34930662 | 2022 | Bowling ve ark., Genet Med (= **SouthSeq**) | **367 bebek**; %30 tanı, %14 belirsiz |
| 29652076 | 2018 | Hiatt ve ark., Clin Genet | 494 birey |
| 37161864 | 2023 | Felker ve ark., Genet Med | "**2999 proband**" boyunca |

Bowling 2017 (PMC5448144) tam metninden, canlı çekildi:

> "systematically shared via ClinVar"

ve aynı makalede: "21 variants of uncertain significance ... in **19 genes not currently associated with
disease**". Yani VUS yatırımı ilan edilmiş bir **veri paylaşım politikasıdır**, nedensellik iddiası değil.

**KIF13A, HudsonAlpha'nın hiçbir yayınında geçmiyor.** Tam metin taramaları: PMC5448144 → 0 kez,
PMC8087252 → 0 kez, PMC11610584 → 0 kez. ~3.000 proband dizilenmiş, aday genler makalede tek tek
sayılmış, KIF13A hiçbirinde yok.

**1. sorunun cevabı: KIF13A seçilmiş bir aday değildi. Kayıt, ilan edilmiş bir toplu paylaşım
politikasının ürünüdür ve "kasıtlı çağrı" okuması çöker.**

---

## 2. Kaydın tamamı: her gözlem alanı

`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=2570672&from_esearch=true`

| Alan | Değer |
|---|---|
| VariationID / Accession | 2570672 / VCV002570672.**1** |
| AlleleID | 2735628 |
| RCV | RCV003311600.1 |
| SCV | **SCV004009653.1**, SubmissionDate **2023-07-07**, DateCreated 2023-07-08 |
| Gönderici | "HudsonAlpha Institute for Biotechnology, HudsonAlpha Institute for Biotechnology", **OrgID 505530**, OrganizationCategory **laboratory**, kısaltma HA |
| Sınıflama | Uncertain significance · criteria provided, single submitter |
| Assertion | "variation to disease" |
| **Assertion criteria** | **"ACMG Guidelines, 2015"**, atıf PubMed **25741868** (Richards ve ark.) |
| **Alel kökeni (Origin)** | **unknown** · de novo DEĞİL, kalıtsal olduğu da belirtilmemiş |
| **Etkilenme durumu (AffectedStatus)** | **unknown** |
| NumberTested | 1 |
| Method | research |
| ObservedData | VariantAlleles = 1 |
| SampleLocalID | `6:17785639:A00166` |
| **Zigotluk** | **YOK** (alan hiç doldurulmamış) |
| **Birlikte bulunan varyantlar** | **YOK** |
| **Aile öyküsü** | **YOK** |
| **Yaş** | **YOK** |
| **Cinsiyet** | **YOK** |
| **Etnisite / popülasyon** | **YOK** |
| **Gönderici yorumu / yorum metni** | **YOK** (aynı partideki 43 kayıtta var) |
| TraitSet | `multipleConditionExplanation="Novel disease"` (aynı partide 25 kayıtta aynı etiket) |
| Fenotipler | HP:0002186 Apraxia · HP:0011098 Speech apraxia · HP:0001250 Seizure · HP:0002342 Moderate intellectual disability · HP:0002308 Chiari malformation, hepsi `ClinicalFeaturesAffectedStatus="present"` |
| StudyName | **CSER-HudsonAlpha** |
| SubmissionName | SUB13646016 |
| Koordinat | GRCh38 chr6:17785639 T>C · SPDI `NC_000006.12:17785638:T:C` · GRCh37 chr6:17785870 |
| dbSNP / ClinGen | **rs2532452593** / CA362830916 |
| DateLastUpdated | 2026-03-21 (ama MostRecentSubmission hâlâ 2023-07-08, yani içerik güncellenmedi) |

**Kalıtım açıkça bilinmiyor.** Bu, tek vakalık bir nörogelişimsel iddia için en kritik eksiktir:
ACMG'nin PS2/PM6 kodları uygulanamaz, PP4 zayıftır, PM2 tek başına kalır.

### 2.1 Mevcut proje dokümanında bir hata

`docs/fonksiyon_kaybi_ekseni.md` satır 332-333 şunu yazıyor:

> "etkilenme durumu **present**"

Bu yanlış. XML'de `<AffectedStatus>unknown</AffectedStatus>`. "present" değeri, ayrı bir alan olan
`ClinicalFeaturesAffectedStatus`'tan geliyor ve o alan yalnızca "listelenen HPO terimleri mevcut" demektir,
"hasta etkilenmiştir" demek değildir. Bu satır düzeltilmelidir.

---

## 3. Varyant nedensel olacak kadar nadir mi?

gnomAD GraphQL, `https://gnomad.broadinstitute.org/api`, dataset `gnomad_r4` (v4.1),
transkript **ENST00000259711** (MANE Select), GRCh38.

### 3.1 Varyantın kendisi

| Ölçüm | Değer |
|---|---|
| Doğrudan varyant sorgusu `6-17785639-T-C` | **`{"errors":[{"message":"Variant not found"}]}`** |
| Allele count (AC) | **0** |
| Allele number (AN) | kayıt yok |
| Allele frequency (AF) | **0** |
| Homozigot sayısı | **0** |
| Popmax frekansı | uygulanamaz (kayıt yok) |
| Filtre durumu | uygulanamaz (kayıt yok) |

**Kapsam kontrolü (önemli):** yokluk, "okunmamış" demek değil. 17785639'un 40 bp çevresinde gnomAD'de
**92 varyant** var; hemen komşu ekzonik pozisyonlar 17785640 (c.3363G>T) ve 17785635 (c.3368C>A) listede
ve bu ekzonda AN yaklaşık **1.453.000**. Bölge okunmuş ve çağrılabilir durumda; varyant gerçekten hiç
görülmemiş.

Ayrıca gnomAD bu bölgede AC=0 olan kayıtları da listeliyor (585 missense kaydı AC=0), yani
"listede yok" gerçekten "hiç gözlenmedi" anlamına geliyor.

### 3.2 Bağlam: KIF13A'da "gnomAD'de yok olmak" varsayılan durumdur

| Ölçüm | Değer |
|---|---|
| Transkriptteki toplam varyant kaydı | 10.936 |
| Missense kaydı | 2.671 |
| ... bunlardan **AC > 0** (gerçekten gözlenmiş) | **2.086** |
| ... AC = 1 (singleton) | 1.047 |
| ... AC >= 10 | 308 |
| En az 1 homozigot taşıyan missense | 32 |
| En yüksek AN | 1.461.712 |
| **CDS'te teorik olarak mümkün missense SNV** (5.418 bp, canlı Ensembl CDS'inden hesaplandı) | **11.918** |
| **Gözlenenlerin oranı** | **2.086 / 11.918 = %17,5** |

**Yani mümkün KIF13A missense SNV'lerinin yaklaşık %82,5'i gnomAD v4.1'de yoktur.** Bir varyantın
gnomAD'de olmaması, bu gen için bilgi taşımayan bir gözlemdir. PM2 gerekli bir filtredir ama ayırt edici
bir kanıt değildir. Bu, iddiaya bir şey kazandırmaz; yalnızca iddiayı öldürmez.

### 3.3 Gen kısıtı: LoF sıkı, missense gevşek

| Ölçüm | Değer |
|---|---|
| pLI | **0,99975** |
| LOEUF (oe_lof_upper) | **0,449** |
| obs/exp LoF | 79 / 212,2 = **0,372** |
| lof_z | 7,76 |
| obs/exp **missense** | 1.901 / 2.209 = **0,861** |
| **mis_z** | **2,79** |

KIF13A güçlü bir **LoF** kısıtı taşıyor, ama missense kısıtı zayıf (oe_mis 0,86). Bir missense VUS'un
bu gende önsel olasılığı düşüktür. pLI 0,9998'i missense bir varyantın lehine kullanmak kategori hatası olur.

### 3.4 Aynı coiled coil'de daha yüksek frekanslı missense'ler

UniProt coiled coil 1100-1138 içinde, gnomAD'de **AC > 0** olan **44 farklı missense**, toplam AC 253:

| Varyant | AC | AN | AF |
|---|---|---|---|
| p.Asp1126Asn | **81** | 1.453.008 | 5,58e-05 |
| **p.Glu1128Lys** (yük tersine çevirmesi) | **27** | 1.453.488 | 1,86e-05 |
| p.Arg1129Trp | 25 | 1.453.352 | 1,72e-05 |
| p.Thr1120Ile | 15 | 1.443.290 | 1,04e-05 |
| p.Ser1117Ile | 13 | 1.453.680 | 8,94e-06 |
| p.Val1116Ile | 12 | 1.456.684 | 8,24e-06 |

Bu coiled coil'deki yük tersine çevirmeleri popülasyonda taşınıyor:
**p.Glu1128Lys AC=27**, p.Glu1124Lys AC=3, p.Glu1107Lys AC=1 ve en çarpıcısı
**p.Lys1104Glu AC=1, yani bizim varyantımızla tam olarak aynı Lys → Glu değişimi, aynı coiled coil'de,
18 kalıntı ötede, sağlıklı bir gnomAD taşıyıcısında.**

Sorusu sorulan "kaç başka missense benzer ya da daha yüksek frekansta" sorusunun cevabı: varyant AF=0
olduğu için gnomAD'de gözlenen **2.086 missense'in tamamı** daha yüksek frekanstadır.

---

## 4. Kalıntı önemli mi?

### 4.1 UniProt Q9H1H9 (1805 aa, MANE ile uyumlu), 1050-1200 aralığındaki tüm özellikler

`https://rest.uniprot.org/uniprotkb/Q9H1H9.json`

| Tip | Aralık | Açıklama |
|---|---|---|
| Chain | 1-1805 | Kinesin-like protein KIF13A |
| **Coiled coil** | **1100-1138** | (kalıntı 1122 **burada**) |
| Alternative sequence | 71-1805 | izoform 5 |
| Alternative sequence | 1071-1083 | izoform 3 ve 4 |
| Sequence conflict | 1185 | Ref. 2 |

Kalıntı 1122'de ya da yakınında **hiçbir modifikasyon bölgesi, bağlanma bölgesi, aktif bölge ya da
fonksiyonel motif yok**. Yalnızca bir coiled coil. Motor alan 5-352'dir; kalıntı 1122 motor alandan
**770 kalıntı uzakta**, kuyrukta. Yapı ekseni bu kalıntıyı değerlendiremez.

Coiled coil bir "alan" değil, düşük karmaşıklıklı bir sarmal sarım bölgesidir; 3.4'te gösterildiği gibi
bu bölge gnomAD'de missense'e açıktır. "Coiled coil içinde" olmak tek başına bir patojenite argümanı değildir.

### 4.2 AlphaMissense: 0,8044 doğrulandı, kalibre eşiğin ALTINDA

Canlı çekim: `https://rest.ensembl.org/vep/human/region/6:17785639-17785639/C?AlphaMissense=1&CADD=1&mane=1&hgvs=1`
ENST00000259711 için dönen blok:

```
"alphamissense": { "am_class": "likely_pathogenic", "am_pathogenicity": 0.8044 }
```

| | Değer |
|---|---|
| AlphaMissense (canlı) | **0,8044** |
| Projenin kaydettiği değer | 0,804 → **DOĞRULANDI** |
| AlphaMissense genel eşiği | 0,564 → geçiyor |
| **Projenin gen-kalibre eşiği** (`data/interim/clinvar_curation/KIF1A_calibration_metrics_20260915.json`, threshold 0,9681) | **0,968 → 0,164 farkla KALIYOR** |
| dbNSFP'deki transkript başına AM skorları | 0,4927 - 0,8044 arası; bir transkriptte sınıf "A" (ambiguous) |

Yani projenin kendi kalibre edilmiş ölçütüne göre bu varyant **aday listesine girmez**.

### 4.3 Diğer prediktörler bölünmüş (MyVariant.info / dbNSFP, canlı)

`https://myvariant.info/v1/variant/chr6:g.17785639T>C?assembly=hg38`

**Zararlı diyenler:** CADD phred **25,2** · PolyPhen-2 HVAR 0,996 D · SIFT 0,001 D · SIFT4G D ·
PROVEAN -3,56 D · MutationTaster D · ESM-1b -14,84 D (rankscore 0,968) · MetaRNN 0,719 D ·
ClinPred 0,970 D · MVP 0,843 · gMVP 0,861

**Tolere edilir diyenler / düşük olanlar:** **FATHMM T** · **PrimateAI 0,673 T** · **DEOGEN2 0,20-0,38 T** ·
**MPC 0,573** (rankscore 0,53; NDD de novo missense'lerde beklenen sinyal MPC > 2'dir) ·
**REVEL 0,642** (yaygın "supporting" eşiği 0,75'in altında) · VEST4 0,62-0,70

Not: REVEL, MetaRNN ve ClinPred etiketle eğitilmiş modellerdir ve projenin kuralına göre konsensüse
girmezler; burada yalnızca tabloyu tamamlamak için listeleniyorlar.

### 4.4 Korunma: kalıntı düzeyinde gerçekten güçlü

**Nükleotid düzeyi (dbNSFP):** phyloP100way **8,017** · phastCons100way **1,0** · GERP++ **5,8** ·
SiPhy 29way 15,81. Ancak phyloP 17way primat yalnızca **0,665**.

**Uyarı:** VEP ve SnpEff'in ikisi de varyantı `missense_variant & splice_region_variant` olarak çağırıyor.
c.3364, ekzon 28/39'un **3. bazıdır** (ekzon c.3362'de başlıyor, c.3362-1 son intronik baz). Bir ekzonun
ilk birkaç bazındaki yüksek phyloP, kısmen kırpılma (splicing) ve ESE kısıtını yansıtır, salt amino asit
kısıtını değil.
Bu, nükleotid korunma sayılarını kalıntı korunması olarak okumayı riskli kılar.

**Kalıntı düzeyi (Ensembl Compara, 203 ortolog, canlı):**
`https://rest.ensembl.org/homology/symbol/human/KIF13A?type=orthologues&sequence=protein&format=full`

İnsan K1122'ye hizalanan kalıntı:

| Kalıntı | Ortolog sayısı |
|---|---|
| **K (Lys)** | **186** |
| boşluk (-) | 11 |
| R (Arg) | 2 |
| E / A / Q / X | her biri 1 |

İnsan bağlamı 1116-1128: `VSNKTEKTEDDVE`. Lizin, şempanze ve fareden tavuğa, Xenopus'a, zebra balığına
(`ILNKPEKSEEDIE`), coelacanth'a ve lampreye (`LVAKPDKSEVEKE`) kadar korunmuş.

Not: hizalama, gnomAD'den bağımsız olarak numaralandırmayı da doğruluyor (T1120, E1124, D1126, V1127,
E1128 bağlamdaki harflerle birebir uyuşuyor).

**Bu, kaydın en güçlü tek verisidir:** kalıntı 1122 gerçekten korunmuş bir lizindir. Ama korunmuş olmak
ile hastalığa neden olmak arasındaki boşluğu 3.4'teki gnomAD tablosu dolduruyor: aynı coiled coil'de
korunmuş kalıntılarda da yük tersine çevirmeleri sağlıklı bireylerde taşınıyor.

---

## 5. Chiari malformasyonu

### 5.1 Genel popülasyon taban oranı (canlı doğrulanmış PMID'ler)

| PMID | Yıl | Çalışma | Bulgu |
|---|---|---|---|
| **21806364** | 2011 | Strahle ve ark., J Neurosurg Pediatr; 14.116 ardışık çocuk, **herhangi bir endikasyonla** MR | **509 (%3,6)** CM-I; yalnızca **%32'si** hekimce semptomatik sayıldı |
| **19433279** | 2009 | Aitken ve ark., Pediatr Neurol; popülasyon tabanlı, 5.248 MR | **51 (%1,0)**; "sıklıkla asemptomatik çocuklarda tesadüfen bulunur" |
| **10839250** | 2000 | Meadows ve ark., J Neurosurg; 22.591 kişi | 175'te >5 mm tonsiller herniasyon (**%0,77**); 25/175 (**%14**) asemptomatik |
| **32113729** | 2020 | Sadler ve ark., Pediatr Neurol | semptomatik ~**1/1.000**, radyolojik ölçütleri karşılayan ~**1/100** |

### 5.2 Gelişimsel gerilik nedeniyle görüntülenen çocuklarda: kritik sayı

**PMID 36548185** (2022), Sätilä ve ark., Neurol Int, PMC9785934: sınırda entelektüel işlevsellik olan
651 çocuk, 263'ü beyin MR'ı çekildi. Makalenin kendi Tablo 4'ünde bulgular
"etiyolojiye götüren / katkıda bulunan / **etiyolojiyle ilgisiz**" olarak sınıflanıyor:

> **Chiari 1: toplam n = 6, altısı da "Not Related to Etiology" sütununda.**

Yani gelişimsel gerilik kohortunda Chiari 1 görüntülenenlerin **%2,3'ünde** bulundu ve
**6/6 vakada (%100) tesadüfi** sayıldı. O kohortta MR'ın toplam etiyolojik verimi yalnızca %11,8'di.

Arka plan: PMID 38349653 (2024, JAMA Netw Open) çocukların **%9,6'sında** tesadüfi bulgu;
PMID 23811972 (2013, AJNR) gençlerin **%10,6'sında**.

### 5.3 Herhangi bir kinezin hastalığında Chiari tanınmış bir özellik mi? Hayır.

- Chiari-1 gen keşif çalışmaları **kromatin yeniden şekillendirme genlerine** işaret ediyor:
  **PMID 39907171** (2025, Birth Defects Res), "Variants in Chromatin Remodeling Genes Are Involved in
  Patients With Chiari Malformation Type 1" · tam metin (PMC11795737) canlı tarandı:
  **KIF13A 0 kez, "KIF" 0 kez, "kinesin" 0 kez**. İlgili genler KMT2A, SETD2, HP1BP3.
  Ayrıca **PMID 42640505** (2026, Brain), "De novo chromatin remodelling variants in sporadic Chiari 1".
- İki güncel kinezin-NDD derlemesi KIF13A'yı hastalık geni listesinden dışlıyor:
  **PMID 42651835** (2026, Curr Issues Mol Biol) hastalık genleri olarak KIF1A, KIF5A, KIF7, KIF11,
  KIF2A, KIF5C, KIF14, KIF15, KIF16B sayıyor, **KIF13A yok**; metinde KIF13A yalnızca bir *fare* fenotipi
  için 2 kez geçiyor, Chiari 0 kez.
  **PMID 40964093** (2025, Front Cell Neurosci) "Intellectual disability" satırında KIF16A, KIF4A, KIF5C
  sayıyor, **KIF13A yok** (yalnızca glioma ve şizofreni yatkınlığı altında geçiyor).
- Chiari, ilgisiz birçok NDD geninde özgül olmayan bir eşlik eden bulgu olarak çıkıyor: RERE, NFIX/Malan,
  SON/ZTTK, PTPN11, FAM111A, HRAS/Costello.

### 5.4 HudsonAlpha'nın kendi verisinde Chiari

1.293 SCV tarandı: yalnızca **3 kayıt** HP:0002308 taşıyor.

| Sınıflama | Çalışma | Örnek | Varyant |
|---|---|---|---|
| Uncertain significance | CSER-HudsonAlpha | A00166 | **KIF13A c.3364A>G** |
| **Likely pathogenic** | PGEN-HudsonAlpha | B000028C | **KMT2C c.7613dup (p.Gln2539fs)** |
| Uncertain significance | AGHI-WGS-HudsonAlpha | 130-G-C | NFKB2 c.2579C>G |

Yani aynı laboratuvarın verisinde Chiari, açıkça patojenik bir **kromatin yeniden şekillendirici**
(KMT2C) frameshift'ine de eşlik ediyor. Bu, 5.3'teki literatürle birebir tutarlı ve Chiari'nin
KIF13A'ya özgü bir imza olmadığını gösteriyor.

**Sonuç:** Chiari malformasyonu fenotip listesinde **vakanın görünür özgüllüğünü şişiriyor**. Bu
kohortta beklenen taban oranı zaten ~%2,3'tür ve o oranın %100'ü tesadüfi sayılmıştır. Listeden Chiari
düşürülürse geriye kalan "apraksi + konuşma apraksisi + nöbet + orta derecede ID" tablosu, herhangi bir
tanısız NDD probandı için son derece sıradan bir HPO setidir.

---

## 6. Bu varyantın ya da bu kalıntının başka bildirimi var mı?

### 6.1 Europe PMC (hepsi sıfır)

| Sorgu | hitCount |
|---|---|
| `KIF13A AND K1122E` | **0** |
| `"K1122E"` | **0** |
| `KIF13A AND Lys1122` | **0** |
| `"p.Lys1122Glu"` | **0** |
| `"c.3364A>G" AND KIF13A` | **0** |
| `KIF13A AND rs2532452593` | **0** |
| `rs2532452593` | **0** |
| `"c.3364A>G"` (çıplak) | 1, ama **PMID 17914064, LRRK2** ile ilgili; farklı gen, rastlantısal koordinat çakışması |

**Bu kalıntı ya da bu varyant hakkında hiçbir yayın yok.**

### 6.2 KIF13A ve nörogelişim literatürünün tamamı

| Sorgu | hitCount |
|---|---|
| `KIF13A` (çıplak) | 649 |
| `KIF13A AND (intellectual disability OR neurodevelopmental OR epilepsy OR seizure)` | 57 |
| `KIF13A AND Chiari` | **2** (ikisi de vaka bildirimi değil: biri FTLD/GRN GWAS, biri ESHG kongre özet kitabı) |
| `KIF13A AND apraxia` | **4** (üçü ESHG özet kitabı derlemeleri) |

**Apraksi + nöbet + orta ID + Chiari kombinasyonunu KIF13A'ya bağlayan sıfır makale var.**

Tek gerçek insan KIF13A nöro bildirimi **PMID 33526817** (2021, Lam ve ark., J Hum Genet):
3 yaşında bir kız, `NM_022113.6:c.2357dupA` **postzigotik mozaik frameshift** (missense değil),
tanımlayıcı fenotip **hypomelanosis of Ito**; **apraksi yok, nöbet yok, Chiari yok**. Bizim kaydımızla
ne varyant sınıfını ne de dört ayırt edici HPO terimini paylaşıyor.

### 6.3 Diğer varyant veri tabanları

| Kaynak | Sonuç |
|---|---|
| **LOVD** (`https://databases.lovd.nl/shared/api/rest.php/variants/KIF13A`) | **11 KIF13A girdisi, hiçbiri c.3364 / p.1122'de değil**; hepsi VKGL-NL "notClassified" |
| **DECIPHER** | Açık API istemci anahtarı istiyor (`{"title":"Forbidden","status":"403","detail":"Client key not provided"}`); varyant düzeyinde doğrulanamadı |
| **G2P / DDG2P** (`https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/`, 2.870 girdi, canlı indirildi) | **KIF13A YOK**. KIF1A iki kez var, ikisi de "definitive" |
| **GenCC** (`https://search.thegencc.org/genes/HGNC:14566`) | **"This gene doesn't have any submissions at this time"** |
| **ClinGen gen geçerliliği** (`https://search.clinicalgenome.org/kb/gene-validity/download`, 3.681 kürasyon) | **KIF13A YOK**. KIF1A = "syndromic intellectual disability, AD, **Definitive**" |
| **Genomics England PanelApp** (`https://panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name=KIF13A`) | `{"count":0}` → **hiçbir tanısal panelde değil** |
| **OMIM** | 605433 yıldızlı (`*`) **gen girdisi**, `#` fenotip girdisi yok |

### 6.4 KIF13A'nın ClinVar profili (kendi doğrulamam)

`KIF13A[GENE]` → **313 kayıt**, bunların **293'ü tek gen** (projenin 293 sayısı **doğrulandı**).

| Sınıflama | Tüm 313 kayıt | Yalnızca 289 nokta varyantı |
|---|---|---|
| Uncertain significance | 221 | **216** |
| (boş) | 62 | 58 |
| **Pathogenic** | **15** | **0** |
| Likely benign | 8 | 8 |
| Benign | 7 | 7 |

**15 "Pathogenic" kaydın hepsi çok genli 6p kopya sayısı varyantıdır** (23 genden 2.586 gene kadar;
örneğin `GRCh37/hg19 6p25.3-q27(chr6:156975-170919482)x3` = 1.030 gen). **ClinVar'da tek bir patojenik
ya da olası patojenik KIF13A nokta varyantı yoktur.**

Fenotip alanları (313 kayıt): "not specified" 217, "not provided" 19, "See cases" 14, artı bu tek kaydın
5 HPO terimi. Yani **"gerçek fenotipli tek kayıt" gözlemi doğrudur**; ama açıklaması "birisi KIF13A'yı
seçti" değil, "diğer 250+ kaydı yatıran tanısal laboratuvarlar fenotip yazmaz, HudsonAlpha ise
kayıtlarının %79,5'ine HPO ekler" olmasıdır.

---

## Bu kayıt ne kadar ileri götürülebilir

### Kanıtın gerçekten izin verdiği en güçlü cümle

> "ClinVar'da, bir araştırma kohortundan gelen tek bir proband, apraksi, konuşma apraksisi, nöbet, orta
> derecede entelektüel yetersizlik ve Chiari malformasyonu ile birlikte KIF13A kuyruğunda gnomAD v4.1'de
> hiç gözlenmemiş, 203 ortologun 186'sında korunmuş bir lizini etkileyen bir missense VUS taşımaktadır
> (VCV002570672, p.Lys1122Glu). Kalıtımı bilinmemektedir, etkilenme durumu kayıtta 'unknown'dır,
> gönderici yorum bırakmamıştır, varyant hakkında hiçbir yayın yoktur ve AlphaMissense skoru (0,804)
> bizim gen-kalibre ettiğimiz eşiğin (0,968) altındadır. Bu gözlem KIF13A'yı nörogelişimsel aday olarak
> **dışlamaz** ve hipotezle tutarlıdır; tek başına hiçbir kanıt ağırlığı taşımaz."

### Aşırı iddia olacak cümleler (kullanılmamalı)

> ~~"KIF13A'yı nörogelişimsel bir gen olarak destekleyen fenotipli bir hasta vardır."~~
> Kayıt bir *hastayı* değil, bir *veri paylaşım politikasının çıktısını* belgeliyor. Etkilenme durumu
> alanı bile "unknown".

> ~~"HudsonAlpha'nın 1.293 gönderiminden tam olarak 1'inin KIF13A'da olması kasıtlı bir aday gen
> çağrısıdır."~~
> Bu ifade HudsonAlpha'nın verisindeki 368 (tek gen bazında) ya da 2.257 (tüm genler bazında) geni aynı
> şekilde tarif eder. En sık desendir, ayrıcalık değildir. Bowling 2017 VUS'ların "systematically shared
> via ClinVar" olduğunu açıkça yazıyor.

> ~~"Chiari malformasyonu vakayı özgül kılıyor / KIF13A'yı posterior fossa gelişimine bağlıyor."~~
> Gelişimsel gerilik nedeniyle görüntülenen çocuklarda taban oran ~%2,3 ve o seride 6/6 tesadüfi.
> Hiçbir kinezin geni Chiari nedeni olarak tanınmıyor; Chiari-1 gen çalışmaları kromatin yeniden
> şekillendiricilere işaret ediyor ve tam metinlerinde "kinesin" hiç geçmiyor.

> ~~"Varyantın gnomAD'de olmaması nedenselliği destekliyor."~~
> KIF13A'da mümkün 11.918 missense SNV'nin ~%82'si gnomAD'de yok. Bu bir filtredir, kanıt değildir.

> ~~"pLI 0,9998 bu missense'in zararlı olduğunu destekliyor."~~
> pLI bir **LoF** ölçütüdür. Aynı genin missense kısıtı zayıftır (oe_mis 0,861, mis_z 2,79) ve varyantın
> bulunduğu coiled coil'de gnomAD'de AC>0 olan 44 missense, bunlar arasında aynı Lys→Glu değişimi
> (p.Lys1104Glu) ve AC=27 ile bir yük tersine çevirmesi (p.Glu1128Lys) vardır.

### Pratik öneri

1. Bu kaydı "fenotiplenmiş destekleyici vaka" olarak değil, **hipotez üreten tek gözlem** olarak
   etiketleyin; ağırlığı ACMG'de en fazla PM2 + PP3(zayıf) düzeyindedir ve PP3 bile gen-kalibre eşiğe
   göre karşılanmıyor.
2. `docs/fonksiyon_kaybi_ekseni.md` bölüm 4.4'teki "etkilenme durumu **present**" ifadesini
   "**AffectedStatus: unknown**; HPO terimleri ClinicalFeaturesAffectedStatus=present" olarak düzeltin.
3. Kalıtımı öğrenmenin tek yolu göndericiye yazmaktır. Susan Hiatt / Greg Cooper (OrgID 505530) örnek
   A00166 için trio verisi olup olmadığını söyleyebilir. de novo çıkarsa kaydın ağırlığı "suggestive"den
   "moderate"a sıçrar; bu, tek başına en yüksek getirili sorudur.
4. `p.Lys1104Glu` (gnomAD AC=1) bu varyant için doğal bir negatif kontroldür; coiled coil ekseninde
   kullanılabilir.

---

## Çekilen canlı uç noktalar

- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=2570672&from_esearch=true`
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=...` (`HudsonAlpha[SBM]`, `"CSER-HudsonAlpha"[STNM]`, `KIF13A[GENE]`, `[ORIG]`, `[COLLMETHOD]`)
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=...` (1.276 + 384 + 313 kayıt)
- `https://www.ncbi.nlm.nih.gov/clinvar/submitters/505530/`
- `https://gnomad.broadinstitute.org/api` (transcript ENST00000259711 + variant `6-17785639-T-C` + gnomad_constraint, dataset `gnomad_r4`)
- `https://rest.ensembl.org/sequence/id/ENST00000259711?type=cds`
- `https://rest.ensembl.org/vep/human/region/6:17785639-17785639/C?AlphaMissense=1&CADD=1&mane=1&hgvs=1`
- `https://rest.ensembl.org/homology/symbol/human/KIF13A?type=orthologues&sequence=protein&format=full`
- `https://rest.uniprot.org/uniprotkb/Q9H1H9.json`
- `https://myvariant.info/v1/variant/chr6:g.17785639T>C?assembly=hg38`
- `https://databases.lovd.nl/shared/api/rest.php/variants/KIF13A`
- `https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/`
- `https://search.thegencc.org/genes/HGNC:14566`
- `https://search.clinicalgenome.org/kb/gene-validity/download`
- `https://panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name=KIF13A`
- `https://www.deciphergenomics.org/api/variants?gene=KIF13A` (403, anahtar gerekli)
- `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=...` (Bölüm 1.5, 5, 6'daki tüm sorgular)

**Doğrulanamayanlar:** DECIPHER varyant düzeyi (istemci anahtarı gerekli); ESHG özet kitaplarındaki
(PMC11627200, PMC6778822) KIF13A ile Chiari eş geçişinin belge içi bağlamı (Europe PMC tam metin için
HTTP 500 döndü; ikisi de yüzlerce bağımsız özetten oluşan derlemeler olduğu için eş geçişin ilgisiz
özetler arasında olması çok olası, ama doğrudan doğrulanamadı).
