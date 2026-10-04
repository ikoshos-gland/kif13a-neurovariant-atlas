# KIF13A vaka arama yol haritasi

**Soru:** Tier A adaylarimizdan birini tasiyan bir insanin norolojik fenotipi var mi, bunu nasil ogreniriz?

**Hazirlanma tarihi:** 2026-09-16. Bu belgedeki her durum satiri o gun canli sorgulandi; sorgu URL'leri
her satirin yaninda. Dogrulayamadigim her sey acikca **DOGRULANMADI** olarak isaretli.

**Bu belge bir plandir, bir islem kaydi degildir.** Bu oturumda hicbir hesap acilmadi, hicbir form
gonderilmedi, hicbir e-posta atilmadi, hicbir yere kayit yapilmadi. Asagidaki taslak metinler bir
insanin (danisman klinisyen veya laboratuvar sorumlusu) gondermesi icin hazirlanmis taslaklardir.

---

## 0. Tek paragrafta cevap

KIF13A **kurumsal gen-hastalik kaynaklarinin neredeyse tamaminda yok** (ClinGen, G2P, PanelApp GE,
PanelApp Australia, OMIM fenotipi: hepsi bos); tek istisna **SysNDD, "Limited"** ve o da tek bir
yayina dayaniyor. Ama en onemli bulgu su: **aradigimiz hasta aslinda zaten ClinVar'da duruyordu ve
bizim tier A listemizde degildi.** 293 tek-gen KIF13A kaydinin 217'si "not specified", 13'u
"not provided", ve **tam ikisi gercek bir norolojik fenotip tasiyor**: (1) HudsonAlpha'nin arastirma
gonderimi **p.Lys1122Glu**, fenotip apraksi + konusma apraksisi + nobet + orta duzey zihinsel
yetersizlik + Chiari malformasyonu; (2) **chr6:17,758,616-18,011,991 delesyonu**, Ensembl ile
dogrulandi ki icinde **tek protein kodlayan gen KIF13A**, ve tasiyici cocuk **etkilenmis**, fenotip
"delayed speech and language development". Ikisinin de ClinVar'da yayinlanmis iletisim adresi var.
Buna karsilik bizim **52 tier A protein degisikligimizin hicbiri hicbir yayinda gecmiyor**, ve
ClinVar'daki 4 adayimizdan ikisi (A254V, R305L) **hasta bile degil, TCGA tumor orneklerinden
hesaplamali splicing anotasyonu**. Bugun izin istemeden **ucret odemeden ve kayit olmadan** gnomAD,
Genebass, AZ PheWAS, All of Us, Open Targets ve GWAS Catalog'un tamami sorgulanabiliyor ve hepsi
sorgulandi: **KIF13A guclu LoF kisitlamasi gosteriyor (pLI 0.9998, LOEUF 0.449) ama 394k eksom ve
484k genomda eksom capinda anlamliliga yaklasan hicbir nadir varyant yuku sinyali yok.**
Kayit gerektiren ama bize acik uc kapi var: **GeneMatcher** (ucretsiz, elle onay yok, yalin gen
sembolu yeterli), **ModelMatcher** (ucretsiz, "in silico" acikca kapsamda) ve **RD-Connect GPAP**
(ucretsiz ama kimlik taramasi ve 10 is gunluk komite karari). Hesap acma karari insanindir.

---

## 1. Rota karar tablosu

| # | Rota | Ne verir | Ne gerektirir | Ucret | Gercekci sure | Bize bugun acik mi? |
|---|---|---|---|---|---|---|
| **R0** | **ClinVar'daki iki fenotipli vakanin gondericisine yazmak** (HudsonAlpha ve ISCA site 1) | Halihazirda var olan, fenotiplenmis iki KIF13A vakasi. En somut hedef | Klinisyen/PI imzasi. HudsonAlpha bir **arastirma** grubu, ticari tani laboratuvari degil, isbirligine cok daha acik | Ucretsiz | Haftalar | **EVET, en yuksek getirili rota** (bkz. 4f) |
| R1 | **GeneMatcher** (genematcher.org) | Gen-bazli eslesme; gonderim sirasinda hangi MME veritabanlarinin sorgulanacagini sen secersin, yani tek hesapla tum federasyona ulasirsin | Self-servis hesap + e-posta aktivasyonu. Elle onay yok, kurumsal e-posta sarti yok. **Yalin gen sembolu yeterli** | Ucretsiz | Kayit dakikalar; eslesme saatler ile aylar | **EVET**, hesap acmak insan karari |
| R2 | **ModelMatcher** (modelmatcher.net) | Scientist Registry seni bulunabilir yapar; Match Dashboard giden MME sorgusunu baslatir | Iki ayri hesap. Yalniz e-posta dogrulamasi. "in silico" acikca uygunluk kapsaminda | Ucretsiz | Kayit dakikalar | **EVET** |
| R3 | **RD-Connect GPAP** | Gercek genom-fenotip platformunda arama; SSS'si hastasiz bilim insanini adiyla davet ediyor | Kimlik taramasi (pasaport), Code of Conduct imzasi, Data Access Committee karari (10 is gunu), 12 aylik yenileme | Ucretsiz | 2 ile 4 hafta | **EVET, ama denetimli** |
| R4 | **DECIPHER, okuyucu olarak** | KIF13A'yi kapsayan 47 CNV hastasi; her kayitta "Initiate contact" formu, DECIPHER ekibi dogrulayip klinisyene iletiyor | **Hesap gerekmiyor.** Yatirma icin klinisyen sifati gerekli, okuma ve temas icin degil | Ucretsiz | Hemen (gezinme), temas icin haftalar | **EVET** (bkz. 6e) |
| R5 | **Ambry'ye yazmak** (G283V, H175R) | Iki tier A varyantimizin klinik testte gorulmus olmasi | Klinisyen/PI imzasi. Toplu VUS yatirimi, fenotip yok, ticari laboratuvar | Ucretsiz | Haftalar, cevap garantisi yok | **Kismen**, R0'dan sonra gelir |
| R6 | **Kohort portallari** (gnomAD, Genebass, AZ PheWAS, All of Us, Open Targets, GWAS Catalog) | Kohort duzeyinde gen-fenotip sinyali | **Hicbiri kayit istemiyor** | Ucretsiz (Genebass toplu indirme haric: requester-pays GCP) | Saatler | **EVET, yapildi** (bkz. bolum 7) |
| R7 | **PhenomeCentral / MyGene2 / seqr MME / PubCaseFinder / PatientMatcher / IRUD / Monarch** | Hasta duzeyinde MME dugumleri | Hasta ve/veya klinik sifat; cogunda elle onay | Ucretsiz | - | **HAYIR** (bkz. 6f). MyGene2'nin yazili bir arka kapisi var |
| R8 | **MME API'sini dogrudan sorgulamak** | - | Kabul edilmis dugum olmak | - | - | **HAYIR**, 401 ile test edildi (bkz. 6g) |
| R9 | **Yayin cikarmak (in silico hipotez makalesi)** | KIF13A'yi haritaya koyar; sonradan gelen vakalarin bize yonelmesini saglar | Yazma emegi | Yok | Aylar | **EVET, en cok kontrolumuzde olan rota** |

---

## 2. KIF13A'nin kurumsal kaynaklardaki canli durumu

Gen kimligi: **KIF13A, HGNC:14566, ENSG00000137177, Entrez 63971, OMIM gen 605433, 6p22.3,
GRCh38 chr6:17,759,126-17,987,635**, MANE Select NM_022113.6 / ENST00000259711.11.
Her satirda pozitif kontrol olarak **KIF1A** kullanildi: ayni sorgu KIF1A icin dolu donuyorsa
KIF13A'nin bos donmesi gercek bir yoklugu gosterir, API hatasini degil.

| Kaynak | KIF13A durumu | KIF1A kontrolu | Sorgu | Erisim |
|---|---|---|---|---|
| **2.1 ClinGen gen-hastalik gecerliligi** | **YOK.** `has_validity=false`, `has_dosage=false`, `has_actionability=false`, `ncurated=0`. Tam gen-validity indirmesinde (dosya olusturma: 2026-09-16, 3682 satir) KIF13A **hic gecmiyor** | "Definitive", syndromic intellectual disability, MONDO:0000508, AD, SOP7, 2020-09-24 | `search.clinicalgenome.org/api/genes?query=KIF13A` ve `/kb/gene-validity/download` | 2026-09-16 |
| **2.2 ClinGen doz duyarliligi** | **YOK.** Haploinsufficiency / triplosensitivity skoru kurate edilmemis | - | ayni | 2026-09-16 |
| **2.3 Gene2Phenotype (EBI)** | **YOK.** `/api/gene/KIF13A/summary/` -> `records_summary: []`; arama -> "No matching Gene found for: KIF13A"; gen stub'inda `last_updated: null` | 2 "definitive" kayit: G2P00959 NESCAV sendromu (monoallelik), G2P00963 herediter duyusal noropati (biallelik), DD paneli | `ebi.ac.uk/gene2phenotype/api/gene/KIF13A/summary/` | 2026-09-16 |
| **2.4 PanelApp Genomics England** | **YOK. Sifir panel.** `count=0`. Hicbir GE panelinde yesil, kehribar veya kirmizi degil | **15 panel** (Herediter spastik paraparezi yesil, Ataksi yesil, Agri sendromlari yesil, ...) | `panelapp.genomicsengland.co.uk/api/v1/genes/?entity_name=KIF13A` | 2026-09-16 |
| **2.5 PanelApp Australia** | **YOK. Sifir panel.** `count=0` | Avustralya orneginde ayrica kontrol edilmedi | `panelapp-aus.org/api/v1/genes/?entity_name=KIF13A` | 2026-09-16 |
| **2.6 SysNDD** | **VAR, "Limited".** `entity_id 3406`, MONDO:0001071 zihinsel yetersizlik, HP:0000006 otozomal dominant, `ndd_phenotype=Yes`, giris tarihi 2021-04-08. Sinopsis aynen: *"de novo postzygotic fs variant in 1 girl with skin hypopigmentation and developmenta delay"*. PMID:33526817 SysNDD yayin tablosunda kayitli. **KIF13A'nin gorundugu tek kurate kaynak budur** | 3 kayit; NESCAV sendromu OMIM:614255 AD "Definitive" (son guncelleme 2026-05-28) | `sysndd.dbmr.unibe.ch/api/gene/?page_size=all` ve `/api/entity/?filter=equals(entity_id,3406)` | 2026-09-16 |
| **2.7 OMIM** | **Sadece gen kaydi, fenotip yok.** MIM 605433, `mim2gene.txt` (uretim 2026-09-15) icinde tip alani `gene`, `gene/phenotype` degil. KIF13A'ya bagli Mendelyan fenotip MIM numarasi bulunamadi | MIM 601255 gen; KIF1A fenotipleri ayri MIM'ler | `omim.org/static/omim/data/mim2gene.txt` | 2026-09-16 |
| **2.8 ClinVar** | **293 tek-gen kaydi, 0 patojenik/muhtemel patojenik.** Tum KIF13A kayitlari 313, missense 232, VUS 217, tek-gen P/LP **0**. Uzman panel yok, uygulama kilavuzu yok. Ayrintili CNV notu icin bkz. bolum 4e | - | E-utilities `esearch db=clinvar term=KIF13A[gene] AND single_gene[prop]` | 2026-09-16 |
| **2.9 DECIPHER** | **Gen sayfasi canli okundu (2026-09-16):** *"DECIPHER holds no open-access sequence variants in this gene"*, ama **55 eslesen hasta varyanti** var, bunlarin **47'si kopya sayisi degisimi**, 8'i "other". Acik erisimli CNV hasta kayitlari gezilebiliyor ve her birinde bir Contact sekmesi var. **JSON API kapali:** butun `/api/` uc noktalari HTTP 403 `{"title":"Forbidden","detail":"Client key not provided"}` donuyor | - | `deciphergenomics.org/gene/KIF13A` | 2026-09-16 |

**Cikarim:** KIF13A "tartismali" degil, **hic degerlendirilmemis**. Bu bizim icin kotu haber degil,
tam tersine projenin var olma sebebi. Ama pratikte sunu demek: bugun bir klinisyen KIF13A varyanti
gorse, hicbir panelde, hicbir G2P kaydinda, hicbir ClinGen kurasyonunda karsiligini bulamaz ve
buyuk ihtimalle raporlamaz. **Vakalar ortada olmayabilir cunku kimse bakmiyor.**

---

## 3. Sorgulanmaya deger tier A varyantlarimiz

`data/processed/KIF13A_candidates_20260915.tsv` icinde **52 tier A** aday var; 51'i motor alanda,
1'i boyun baglayicida (L370P). Hepsi KIF1A ile ayni rezidude ve o rezidude KIF1A'nin patojenik
varyanti mevcut. **25 tanesi gnomAD'da hic gorulmemis** (de novo adayi olarak en degerli grup).

### 3a. Eslestirme sorgusuna gonderilecek oncelikli liste

Bir eslestirme gonderiminde varyant listesi vermek isterseniz, once bu 12'yi verin: AlphaMissense
en yuksek, gnomAD'da yok veya cok nadir, ve KIF1A'da ayni rezidude patojenik varyant var.

| Protein | MANE HGVSc | gnomAD AC | AM | KIF1A rezidusu | KIF1A'da o rezidude patojenik |
|---|---|---|---|---|---|
| G283W | c.847G>T | 0 | 0.9997 | G279 | G279R |
| R173P | c.518G>C | 0 | 0.9995 | R169 | R169W, R169T (fonksiyonel), R169K |
| R11P | c.32G>C | 0 | 0.9994 | R11 | R11L, R11Q (fonksiyonel), R11W |
| G107V | c.320G>T | 0 | 0.9993 | G102 | G102S, G102D |
| **E19K** | c.55G>A | 0 | 0.9992 | E19 | **E19K, birebir ayni degisiklik** |
| R305L | c.914G>T | 0 | 0.9990 | R307 | R307P, R307G, R307Q |
| R11G | c.31C>G | 0 | 0.9985 | R11 | R11L, R11Q, R11W |
| S63P | c.187T>C | 0 | 0.9984 | S58 | S58W, S58L |
| R13G | c.37C>G | 0 | 0.9981 | R13 | R13L, R13C, R13S, R13H |
| S278W | c.833C>G | 0 | 0.9978 | S274 | S274L (fonksiyonel) |
| L370P | c.1109T>C | 1 | 0.9999 | L372 | L372V |
| G283V | c.848G>T | 3 | 0.9992 | G279 | G279R |

**En guclu tek argüman E19K'dir:** KIF1A'da p.Glu19Lys patojenik olarak biliniyor ve KIF13A'daki
birebir ayni aminoasit degisikligi gnomAD'da hic gorulmemis. Bir eslestirme metninde tek cumlede
anlatilabilecek tek nokta budur.

*(Tam 52'lik liste icin `data/processed/KIF13A_candidates_20260915.tsv`, `tier == "A"`.)*

---

## 4. ClinVar'daki 4 aday ve gonderici rotasi

Tablo: `data/interim/casefinding/KIF13A_clinvar_tierA_submitters_20260916.tsv`

### 4a. Ne bulduk

| Varyant | VCV / SCV | Gonderici | Yontem | Gercek hasta gozlemi mi? |
|---|---|---|---|---|
| **G283V** c.848G>T | VCV004092057 / SCV006347025 | **Ambry Genetics** (OrgID 61756) | clinical testing, germline, affected status **unknown**, kosul "not specified" | **EVET** ama fenotip yatirilmamis |
| **H175R** c.524A>G | VCV003533722 / SCV005608632 | **Ambry Genetics** (OrgID 61756) | clinical testing, germline, affected status **unknown**, kosul "not specified" | **EVET** ama fenotip yatirilmamis |
| **A254V** c.761C>T | VCV004323816 / SCV006649357 | **Dr. Peter K. Rogan Lab, Western University** (OrgID 505643) | **in vitro**, origin "not applicable", affected status "not applicable", siniflandirma yok (`evidence_only`) | **HAYIR** |
| **R305L** c.914G>T | VCV004323813 / SCV006649354 | **Dr. Peter K. Rogan Lab** (OrgID 505643) | ayni | **HAYIR** |

### 4b. Rogan Lab kayitlarindaki tuzak (onemli duzeltme)

A254V ve R305L "siniflandirilmamis arastirma gonderimi" degil, **hasta bile degil**. ClinVar XML'inde
gozlem aciklamasi aynen soyle:

> "In the sample (TCGA-BJ-A45K), the counts of splice junction spanning reads and abundance of the
> transcript adjacent to this mutation were compared with RNASeq data from matched tissues and tumors
> lacking the same mutation using Veridical (PMID:24741438)"

R305L icin ayni cumle, orneklem **TCGA-DE-A2OL**. Ikisi de PMID **31275557**'ye atif yapiyor:
Shirley, Mucaki & Rogan, *Pan-cancer repository of validated natural and cryptic mRNA splicing
mutations*, F1000Research 2018, DOI 10.12688/f1000research.17204.3. Bu, TCGA/ICGC tumor eksomlarinda
splicing etkisi tahmin edilen 341,486 mutasyonun toplu yatirimidir. Rogan Lab'in ClinVar'a toplam
**523,700** gonderimi, yalniz KIF13A'da **164** gonderimi var.

**Sonuclar:**
1. **Bu iki varyant icin aranacak bir hasta yoktur.** Norolojik fenotip bilgisi sifirdir.
2. Bunlari "ClinVar'da bulunuyor" diye kanit sutununa yazmak yanlis olur; yalnizca o pozisyonda bir
   tumorde bir nukleotid degisimi gozlendigini gosterir.
3. Bizim projede bu iki varyantin tier A olmasi degismez (AlphaMissense + KIF1A paralog argumani
   ayakta), ama **insan kaniti sutunu bos kalmalidir**.

### 4c. Ambry kayitlarinin degeri ve siniri

Ambry'nin yorumu iki varyantta da ayni sablon: *"Based on insufficient or conflicting evidence, the
clinical significance of this alteration remains unclear."* Fenotip yok, kosul "not specified",
etkilenme durumu "unknown". Ambry'nin ClinVar'a toplam gonderimi **1,801,188**, yalniz KIF13A'da
**217** gonderim (son guncelleme 2026-03-31). Yani bunlar hedefli bir KIF13A calismasi degil,
klinik eksom/panel testinden cikan toplu VUS yatirimidir.

Yine de **G283V ve H175R gercekten bir insanda goruldu.** Bu iki varyanti tasiyan kisiler bir sebeple
klinik genetik testine girmisti. Fenotipleri ClinVar'da degil, Ambry'nin kendi veritabaninda.

### 4d. ClinVar'in gondericiler hakkinda acikladigi seyler

ClinVar her gonderici icin bir kurum sayfasi yayinliyor
(`ncbi.nlm.nih.gov/clinvar/submitters/<OrgID>/`) ve orada **acik adres, web sitesi, GTR kaydi,
kullandigi siniflandirma kriterleri, toplam gonderim sayisi, gen bazinda gonderim sayisi ve isimli
bir irtibat kisisi (telefon + e-posta)** bulunuyor:

- **Ambry Genetics**, OrgID 61756, 7 Argonaut, Aliso Viejo CA 92656, ABD. Irtibat: "Variant HelpDesk",
  `(adres ClinVar başvuru sahibi sayfasında)`, +1 949-900-5500.
- **Dr. Peter K. Rogan Lab**, OrgID 505643, Western University, 1151 Richmond St, London ON N6A3K7,
  Kanada. Irtibat: Peter Rogan (Lab Director), `(adres ClinVar başvuru sahibi sayfasında)`; Ben Shirley (informatics),
  `(adres ClinVar başvuru sahibi sayfasında)`.

**Ama ClinVar'in kendi yonlendirme/aracilik servisi yoktur.** ClinVar FAQ'te gondericiye ulasma
konusunda hicbir yonerge yok; tek onerdikleri sey siniflandirmaya katilmiyorsaniz kendi kanitinizi
gondermenizdir: *"If you think a variant in ClinVar has been classified incorrectly, we encourage you
to submit your own classification of the variant along with your evidence, such as recent
publications."*

Var olan resmi mekanizma **ClinGen Sequence Variant Inter-Laboratory Discrepancy Resolution**'dir
(`clinicalgenome.org/tools/variant-interpretation-discrepancy-resolution-module/`), fakat o
**laboratuvarlar arasi** isler ve ayni varyanti farkli siniflandiran iki gonderici gerektirir.
Bizim durumumuzda tek gonderici var ve biz laboratuvar degiliz, dolayisiyla **bu kapi bize kapali**.

**Gercekci beklenti:** klinik bir tani laboratuvarina "su VUS'u hangi hastada gordunuz, fenotipi
neydi" diye ucuncu taraf bir hesaplamali arastirmaci olarak yazarsaniz, hasta duzeyinde bilgi
alamazsiniz; onay ve gizlilik cercevesi buna izin vermez. Anlamli cevap gelme ihtimali ancak
**bir klinisyen/kurum ortagi uzerinden, isbirligi teklifi seklinde** gidildiginde vardir.
*(Bu son cumle makul bir beklentidir, yayinlanmis bir Ambry politikasindan alinti degildir:
**DOGRULANMADI**.)*

### 4e. ClinVar'daki "patojenik" KIF13A kayitlari: hepsi buyuk CNV

Brifingteki "sifir patojenik missense" ifadesi dogru, ama tam resim sudur: **`KIF13A[gene]` sorgusu
15 adet patojenik/muhtemel patojenik kayit dondurur ve bunlarin tamami cok genli buyuk kromozomal
kopya sayisi degisimidir.** Tek-gen filtresi (`single_gene[prop]`) uygulandiginda patojenik sayisi
**sifira** duser.

| Sorgu | Sayi |
|---|---|
| `KIF13A[gene]` toplam | 313 |
| `KIF13A[gene] AND single_gene[prop]` | 293 |
| `... AND single_gene[prop]` missense | 232 |
| `... AND single_gene[prop]` VUS | 217 |
| `... AND single_gene[prop]` patojenik/muhtemel patojenik | **0** |
| `KIF13A[gene]` patojenik/muhtemel patojenik (tek-gen filtresiz) | **15** |
| `KIF13A[gene]` copy number loss | 6 |
| `KIF13A[gene]` copy number gain | 13 |

Ornek patojenik kayitlar: `VCV004755369` 6p24.1-22.3 chr6:12,020,008-22,846,454 x1 (10.8 Mb kayip),
`VCV001810398` 6p25.3-22.3 chr6:820,000-21,700,000 x3, `VCV000442457` 6p23-22.3
chr6:13,693,852-24,225,515 x1. Bunlarin hepsinde KIF13A onlarca genden biridir; patojeniklik
atfi KIF13A'ya degil, delesyona/duplikasyona aittir.

Iki kayit KIF13A'yi **kismen** kapsiyor ve daha ilginc olabilir, ama ikisi de VUS:
`VCV000057174` 6p22.3 chr6:17,397,932-17,966,557 x3 (GRCh38) ve `VCV000253699` 6p22.3
chr6:17,637,157-18,931,615 x3 (GRCh37).

**Cikarim:** ClinVar'da KIF13A'yi hastalikla iliskilendiren tek bir tek-gen patojenik kayit yoktur.
Butun "patojenik" etiketler bitisik-gen olaylarindan gelir. **Tek istisna asagidaki 4f'dir ve o
patojenik degil, VUS'tur, ama fenotipi vardir.**

### 4f. ASIL BULGU: ClinVar'da fenotipi olan iki KIF13A vakasi var

Tablo: `data/interim/casefinding/KIF13A_phenotyped_leads_20260916.tsv`

293 tek-gen KIF13A kaydinin kosul dagilimini cikardik. Sonuc:

| Kosul | Kayit sayisi |
|---|---|
| "not specified" | 217 |
| "not provided" | 13 |
| "See cases" | 1 |
| **Apraxia / Speech apraxia / Seizure / Moderate intellectual disability / Chiari malformation** | **1** |

Ayrica tek-gen filtresi disinda, **KIF13A'dan baska protein kodlayan gen icermeyen bir delesyon**
daha var ve onun da fenotipi kayitli. Bu iki kayit, projenin sordugu sorunun bugun var olan en iyi
cevabidir ve **ikisi de bizim tier A listemizde degildi.**

#### Lead 1: p.Lys1122Glu, HudsonAlpha arastirma gonderimi

| Alan | Deger |
|---|---|
| Varyant | `NM_022113.6(KIF13A):c.3364A>G (p.Lys1122Glu)` |
| Kayit | VCV002570672 / SCV004009653, olusturma 2023-07-08, gonderim 2023-07-07 |
| Siniflandirma | Uncertain significance, "criteria provided, single submitter" |
| **Fenotip** | **Apraxia; Speech apraxia; Seizure; Moderate intellectual disability; Chiari malformation** |
| Yontem | **research** (klinik test degil), NumberTested 1 |
| Gonderici | **HudsonAlpha Institute for Biotechnology**, OrgID 505530 |
| Iletisim | Susan Hiatt `(adres ClinVar başvuru sahibi sayfasında)` +1 256-327-0932; Greg Cooper (PI) `(adres ClinVar başvuru sahibi sayfasında)` +1 256-327-9490; Don Latner `(adres ClinVar başvuru sahibi sayfasında)`. 601 Genome Way Suite 2005, Huntsville AL 35806, ABD |

**Neden bu kayit degerli:** HudsonAlpha'nin ClinVar'a toplam **1,293** gonderimi var ve
**KIF13A'da tam olarak 1 gonderimi** var. Ambry'nin 217, Rogan Lab'in 164 toplu gonderimiyle
kiyaslayin. Yani bu bir toplu VUS yatirimi degil, **fenotiplenmis bir probandda bilincli bir aday gen
cagrisi**. Ustelik ticari bir tani laboratuvari degil, bir **arastirma** grubu: isbirligi teklifine
cevap verme olasiligi cok daha yuksek ve hasta duzeyinde konusma cercevesi cok daha musait.

**Durust sinir:** 1122. pozisyon **motor alanin disinda** ve bu varyant **bizim tier A listemizde
degil** (bizimkilerin hepsi 1-370 arasinda). Yani bu kayit KIF13A'yi bir norogelisimsel aday gen
olarak destekler, **bizim motor-alan missense hipotezimizi dogrudan desteklemez.**

#### Lead 2: Yalnizca KIF13A'yi silen delesyon, etkilenmis cocuk

| Alan | Deger |
|---|---|
| Varyant | `GRCh37/hg19 6p22.3(chr6:17758616-18011991)x1`, yaklasik **253 kb heterozigot delesyon** |
| Kayit | VCV000442594 / SCV000585270, gonderim 2017-09-26 |
| Siniflandirma | Uncertain significance, **"no assertion criteria provided"** |
| **Fenotip** | **Delayed speech and language development** |
| **Etkilenme durumu** | **YES** |
| Yontem | clinical testing; Microarray + Fluorescence in situ hybridization |
| Gonderici | **ISCA site 1**, OrgID 505237, Bethesda MD |
| Iletisim | Erin Riggs, Coordinator, `(adres ClinVar başvuru sahibi sayfasında)` |
| Baglanti | Kaminsky ve ark. 2011, PMID 21844811; dbVar nstd101 |

**Delesyonun icerigi bagimsiz olarak dogrulandi.** Ensembl GRCh37 REST ile
`6:17758616-18011991` araligindaki genler sorgulandi ve **tek protein kodlayan gen KIF13A**
(17,759,414-17,987,854); yaninda yalniz iki pseudogen var (ENSG00000203518, ENSG00000220920).

**Neden bu kayit degerli:** Bu, **norogelisimsel fenotipi olan etkilenmis bir bireyde, pratikte
KIF13A'dan baska hicbir seyi silmeyen bir delesyon**. Bolum 5b'deki 13 genli bitisik-gen sendromunun
aksine, burada karisik degisken yok. **KIF13A haploinsufficiency'si icin bugun var olan en guclu
insan kaniti budur ve ClinVar'da 2017'den beri, hicbir kriter uygulanmadan "VUS" olarak, kimsenin
dikkatini cekmeden duruyor.**

**Durust sinir:** kalitim bildirilmemis (de novo mu, ebeveynden mi bilinmiyor), siniflandirma kriteri
yok, ve bu bir **LoF** kanitidir. Lam 2021 ve SysNDD'nin "Limited" otozomal dominant kaydiyla
tutarlidir; **motor-alan missense hipotezimizi test etmez.**

#### Bu iki kaydin projeye etkisi

1. **Vaka arama hedefimiz degisti.** En somut iki hedef R0 rotasidir (bolum 1), bizim ClinVar'daki
   tier A varyantlarimiz degil.
2. **Iki kayit da LoF/haploinsufficiency yonune isaret ediyor**, bizim missense hipotezimize degil.
   Bu, gnomAD kisitlama profiliyle de uyumlu (bolum 7a: guclu LoF kisitlamasi, mutevazi missense
   kisitlamasi). Bu gerilimi rapor ederken saklamayalim.
3. **Her iki kayit da VUS ve kriter uygulanmamis.** Yani KIF13A icin ClinVar'da bekleyen, henuz
   kimsenin toplamadigi bir kanit yigini var.

## 5. Yayinlanmis insan vakalari

Tablo: `data/interim/casefinding/KIF13A_literature_human_cases_20260916.tsv`
Arama: Europe PMC REST, 2026-09-16. `KIF13A` toplam 649 kayit, baslik/ozet'te 649'un 86'si.

### 5a. Lam 2021, dogrulandi

**PMID 33526817**, DOI 10.1038/s10038-020-00883-w, *J Hum Genet* 66 (2021), Lam CW, Chan CY,
Wong KC, Chang ST. *"Postzygotic inactivating mutation of KIF13A located at chromosome 6p22.3 in a
patient with a novel mosaic neuroectodermal syndrome."*

- Hasta: 3 yasinda kiz, karakteristik deri hipopigmentasyonu (hipomelanozis of Ito) + gelisim geriligi.
- Varyant: **NM_022113.6: c.2357dupA**, heterozigot de novo frameshift, NMD ongoruluyor.
- Dusuk mutant allel orani -> postzigotik, embriyonik mozaiklik.

**Bizim adaylarimizla ortusmuyor:**
- Mekanizma farkli: LoF/haploinsufficiency, bizimkiler missense.
- Bolge farkli: c.2357 yaklasik 786. kodon, **motor alanin disinda**; bizim 52 adayin tamami
  1-370 arasinda.
- Yani bu vaka KIF13A'nin norogelisimle ilgili oldugunu destekler ama **bizim motor-alan
  missense hipotezimizi ne dogrular ne de yalanlar.**

### 5b. 6p22.3 delesyon sendromu (baglamsal)

**PMID 34303382**, DOI 10.1186/s13039-021-00557-y, *Mol Cytogenet* 2021, PMC8310580.
2.8 yasinda erkek, `arr[GRCh37] 6p22.3(15,794,379_18,277,334)x1 dn`, 2.5 Mb, de novo, **13 gen**
iceriyor ve KIF13A (GRCh37 chr6:17,759,357-17,987,866) bunlardan biri. Fenotip: gelisim geriligi +
hafif dismorfizm. Makale, KIF13A'yi 6p22.3 delesyon fenotipinde "en tutarli sekilde implike edilen"
genler arasinda sayiyor ve KIF13A icin pLI 1.00 veriyor. Ayni makale, 2021'de literaturde **19**,
DECIPHER'de **18** ortusen 6p22.3 delesyon vakasi bulundugunu yaziyor.

**Sinir:** bu bir bitisik-gen sendromu. 13 genin hangisinin ne yaptigini soylemez ve tam gen
delesyonudur; motor alan missense'lerimiz icin dogrudan kanit degildir.

### 5c. Bizim tier A protein degisikliklerimiz: sifir

Europe PMC'de `KIF13A AND (G283V OR H175R OR A254V OR R305L OR L370P OR G283W OR D306H OR R173P OR
E19K OR H222R)` -> **hitCount 0**. Uc harfli `p.Gly283Val`, `p.His175Arg`, `p.Ala254Val`,
`p.Arg305Leu`, `p.Leu370Pro`, `p.Gly283Trp` formlari -> 17 kayit, **hepsi baska genlerde**.

**Bizim adaylarimizdan hicbiri hicbir yayinda gecmiyor.**

### 5d. Negatif kontrol: KIF13A hala haritada degil

**PMID 42651835** (2026-08-18), *Curr Issues Mol Biol*, *"Traffic Jams in the Brain: How Kinesin
Dysfunction Shapes Neurodevelopmental Disorders"*: 2026 tarihli guncel kinesin-NDD derlemesi hastalik
genlerini **KIF1A, KIF5A, KIF7, KIF11, KIF2A, KIF5C, KIF14, KIF15, KIF16B** diye sayiyor.
**KIF13A listede yok.** Ayni sekilde PMID 40964093 (2025) kinesin-norolojik hastalik derlemesinde de
ozet duzeyinde gecmiyor.

### 5e. Preprint taramasi

`KIF13A AND SRC:PPR` -> 22 preprint, hicbiri KIF13A norogelisimsel insan vakasi bildirmiyor
(cogu kanser, hucre biyolojisi, COVID). En yakin: 2026-01-21 tarihli "Systematic assessment of rare
and de novo structural variants in 57 patient-parent trios using optical genome mapping" preprinti;
KIF13A'ya ozgu bir bulgu icerdigi **DOGRULANMADI**.

---

## 6. Matchmaker Exchange: hangi kapi acik, hangisi kapali

Tablo: `data/interim/casefinding/matchmaking_nodes_20260916.tsv`

`matchmakerexchange.org/participants.html` (erisim 2026-09-16) 12 katilimci listeliyor.
Resmi dugum listesi (`i_am_a_clinician_laboratory.html`) sekiz tanedir: DECIPHER, GeneMatcher, IRUD,
MyGene2, PatientMatcher, PhenomeCentral, RD-Connect GPAP, seqr. Genis katilimci sayfasi bunlara
BCH seqr, Broad seqr, ModelMatcher, Monarch Initiative ve PubCaseFinder'i ekliyor.

MME'nin olcegi (PMID 35537081, *Hum Mutat* 2022, DOI 10.1002/humu.24373): **13,520'den fazla essiz
gen (tum protein kodlayan genlerin yaklasik %68'i)**, 98 ulkede 12,000'den fazla katkici,
120,000'den fazla vaka. Yani **KIF13A'yi zaten birinin kaydetmis olma ihtimali ciddidir ve bunu
ogrenmenin tek yolu kaydolmaktir.**

### 6a. GeneMatcher: ucretsiz, elle onay yok, yalin gen sembolu yeterli

> **Onemli teknik not:** genematcher.org, WebFetch'e Cloudflare HTTP 403 donuyor ama **tarayici
> User-Agent'i ile curl'e HTTP 200 donuyor.** Asagidakilerin tamami **sitenin kendi sayfalarindan,
> 2026-09-16'da canli okundu**; arsiv kopyasi kullanilmadi. FAQ'in dinamik atif satiri
> "Accessed September 16th, 2026" yaziyor, `/statistics` sayfasi "September 16, 2026" tarihli.
> **Hicbir hesap acilmadi, hicbir form gonderilmedi.**

**Kimler uye olabilir** (`genematcher.org/help/faq`, SSS Q2, aynen): *"Who can be a part of the
GeneMatcher community? / Researchers / Model Organism Researchers / Clinicians/Counselors/Other
health care providers / Patients and their families"*. **Arastirmacilar listenin basinda; hasta
sarti yok, klinisyen sarti yok.**

**Kayit** (`genematcher.org/account/`, tarayicida goruldu): e-posta, sifre, opsiyonel API anahtari,
aciklama, rol kutulari (Researcher / Model Organism Researcher / Health Care Provider / Patient),
ad, soyad, derece, kurum, adres, sehir, eyalet, posta kodu, ulke, telefon. Sayfanin kendi ifadesi:
*"All fields are required except for the API key, the description, the degree and the second address
line."* ve *"Once the account is created you will be emailed a link which will allow you to activate
the account."* **Elle onay sureci yok, kurumsal e-posta sarti yok.** Kurum alani zorunlu ama
serbest metin.

**Gen-tek-basina gonderim: EVET** (SSS Q9, aynen): *"At least one of the following pieces of
information is required: MIM number ... Gene Symbol ... Genomic Location ... Phenotypic Features"*.
**Yalin bir gen sembolu yeterlidir.** About sayfasi: *"Users have the option, though are not
required, to provide a variant (or variants) ... as well as to submit clinical features of the
patient/family"*.

**Ucret: yok.** About sayfasi: *"GeneMatcher is a freely accessible web site"*.
**Ama bir celiski var, bilerek not ediyorum:** EULA metni (`genematcher.org/help/eula`)
*"Upon User's payment of the required fees, if any"* diyor ve uygulamayi *"intended solely for
personal, non-commercial use in a clinical research setting by an appropriately trained clinical
researcher ... and is not intended for use by an individual patient"* diye tanimliyor. Bu madde
sitenin kendi SSS'siyle celisiyor (SSS hastalari ve model organizma arastirmacilarini kabul
ediyor). Standart JHU lisans metni olmasi kuvvetle muhtemel ve bir uygunluk engeli olarak
uygulanmiyor, **ama hesabi acacak kisinin bunu gormus olmasi gerekir.**

**Eslesme bildirimi** (SSS Q6, aynen): *"Users will receive an email informing them of the presence
or absence of a match in the database. If there is currently no match in the database, the user will
be notified and the submission will be stored in the database unless or until it is deleted by the
user. In the future, if another user creates a submission that matches, both users will be notified
of the match simultaneously."* ve *"The email about the match includes the information that matched
and the email contact of the users involved."* **Gonderim silinene kadar eslesmeye devam eder.**

MME katilimci sayfasi: *"Matches are automatically made when the user saves their submission, and
can initiate a match against MME at any time, selecting which MME servers to match against."*
**Yani tek bir GeneMatcher hesabi ile DECIPHER, PhenomeCentral, MyGene2, seqr, PatientMatcher,
RD-Connect ve IRUD'a sorgu gonderilebiliyor.** (Ayni sayfa: *"GeneMatcher excludes patients from
matching when receiving MME match requests, and does not allow patients to make MME match
requests."*)

**Guncel olcek** (`genematcher.org/statistics`, 16 Eylul 2026): 117 ulkeden **20,310 gonderici**,
**116,446 gonderim**, **17,144 essiz gen**, **12,831 gen en az bir kez eslesmis**. Insan disi
gonderimler de var (fare 283, Drosophila 36, zebra baligi 32, C. elegans 11, maya 1), yani
platformu klinik disi arastirmacilar fiilen kullaniyor.

**Yayinlanmis olcek** (Hamosh ve ark., *Hum Mutat* 2022;43(6):668-673, **PMID 35170833**,
PMC9133194, Europe PMC tam metninden birebir): *"As of October 1, 2021, there are 12,531 submitters
from 94 countries who have submitted 58,134 submissions with 13,498 unique genes in the database.
Among these genes, 8970 (64%) have matched at least once and the total number of matches is 378,806,
growing by about 10,000 per month."* Ayni makale: *"All submissions in GeneMatcher are available for
matching across the MME."*

**Kurucu makale:** Sobreira N, Schiettecatte F, Valle D, Hamosh A, **PMID 26220891**, *Hum Mutat*
2015, DOI 10.1002/humu.22844.

**Ayni grubun varyant duzeyinde araci:** VariantMatcher, *"developed to enable public sharing of
variant-level data and phenotypic information"* (Wohler ve ark., **PMID 34407837**,
*Orphanet J Rare Dis* 2021, PMC8371856). Bolum 3a'daki 12 varyanti buraya koymak gen kaydini
tamamlar.

### 6b. ModelMatcher: ucretsiz, "in silico" acikca kapsamda, iki hesap gerekiyor

`modelmatcher.net/instructions.html` (erisim 2026-09-16), aynen: *"If you are a scientist willing
and able to perform experiments on a gene, you are eligible to register as a 'Scientist User' in
ModelMatcher."* Sozluk tanimi: *"Scientists who register in ModelMatcher are willing to use in vivo,
in vitro, and/or in silico methods to experimentally study the function of genes and variants linked
to human diseases."* **"in silico" acikca kapsam icinde.**

**Hasta hicbir zaman gerekmiyor: vaka degil GEN kaydedilir.** *"You can list as many genes of
interest as you like for each model organism"*, her gen icin Primary/Secondary kademe ve
Public/Restricted gizlilik. Tutorials sayfasi: *"You can add any gene that you are interested in
establishing a collaboration for."*

**Eslestirme iki yonlu ama iki farkli hesapla isliyor. Ikisini de acmak gerekiyor:**
- **Scientist Registry** uzmanlik genlerinizi tutar ve **gelen** talepleri alir: *"After you register
  your gene of interest, you may receive messages from our Match Dashboard or partner registry users
  who want to discuss a collaboration with you."*
- **Match Dashboard** **giden** MME sorgusunu baslatir: *"Your 'Match with Clinicians and Patients'
  tab allows you to identify clinicians and patients registered in Matchmaker Exchange who share
  common genes of interest. To send match requests, click the Exchange info button next to the gene,
  input your reason for contact, select the partner registries you want to match with, and hit the
  send button."*
- Site oneriyor: *"A different email and a password can be used to register in the Scientist Registry
  and Match Dashboard but we recommend to keep this the same to avoid any confusion."*

MME katilimci sayfasi dogruluyor: *"ModelMatcher users (scientists, clinicians, and other users) can
choose to initiate a match using their registered genes by selecting the external databases they
would like to connect with. Users are notified of results as soon as they are received."*

**Ucret:** Match Dashboard icin aynen *"Anyone is allowed to make one of these accounts. It is
completely free and there are no required qualifications, just an email address verification step."*
Scientist Registry icin hicbir sayfada ucret ibaresi yok; **ucretsiz oldugu ucretin yoklugundan
cikarilmistir, acik bir alinti degildir.**

**Erisim genisligi:** Public isaretli bir gen RDMM (Kanada), Solve-RD (Avrupa), AFGN (Avustralya) ve
MME ortaklari GeneMatcher, PhenomeCentral, MyGene2 ile paylasilir; ayrica hesabi olmayan herkesin
public arama sayfasindan bulunabilir. Restricted genler yalniz Match Dashboard uzerinden
kesfedilebilir ve siz cevap verene kadar kimliginiz gizli kalir.

Kaynak makale: Harnish JM ve ark., **PMID 35224820**, *Hum Mutat* 2022, DOI 10.1002/humu.24364.
Ozetin kendi ifadesiyle mevcut platformlarin sorunu: *"these platforms exclude some stakeholders
based on their role or geography"* -- yani tam bizim durumumuz icin yapilmis.

### 6c. RD-Connect GPAP: hastasiz bilim insanini adiyla davet eden tek dugum

`rd-connect.eu/what-we-do/omics/gpap/faq-2/`, SSS 1.1.7. **Sorunun kendisi aynen:** *"I am a basic
scientist and don't have any human data to share. Can I still register?"* Cevap aynen:

> *"Yes! If you are interested in finding human data analogues for your gene or model of interest or
> information about whether a specific variant you are working on in an animal model causes a similar
> phenotypic presentation in a human, please consider registering. You will have to state your
> research interests during the signup process but we do not require you to submit data."*

**Ucretsiz** (SSS 1.1.3: *"Currently, all RD-Connect services are entirely free for all users"*) ama
**denetim agir** (SSS 1.2.1): *"We'll ask you to provide proof of your credentials as legitimate
researchers, a scan of your ID (e.g. passport), and explain your research interests ... You will
also need to read and accept the RD-Connect Code of Conduct by signing an Adherence Agreement ...
Your application will be considered by the Data Access Committee (DAC) and you will receive a
response within 10 working days."* Erisim **12 aylik** dilimlerle veriliyor ve kullanim izleniyor.
Kurumsal baglanti sart (SSS 1.5.1: *"associated with a recognised academic institution"*).

**Sinir:** verisiz kayitli bir kullanici olarak GPAP'in kendi veri setlerinde **arama** hakki
kazanirsiniz; MME'ye gonderim yapmak hasta kaydi gerektirir (SSS 1.4.3: MME katilimi veri
yatirandan tetiklenir). Bu **CIKARIM**, acik bir alinti degil.

Bu, uc acik kapidan **en agir surtunmeli ama en yuksek kaliteli** olani: GPAP gercek genom-fenotip
verisi olan bir platform, sadece bir ilan panosu degil.

### 6d. DECIPHER: hesap olmadan da klinisyene ulasilabiliyor

Onemli ayrim: DECIPHER'e **hasta yatirmak** icin klinisyen olmak sart
(`deciphergenomics.org/join`): proje *"part of an academic centre affiliated to a hospital"* olmali,
*"led by a senior clinician"*, lab ve klinisyenler meslek orgutu uyesi olmali, kurumsal e-posta
zorunlu, ve her bireysel hesabi sorumlu klinisyen onayliyor.

**Ama okumak ve temas kurmak icin hesap gerekmiyor.** Oturum acmadan goruntulenen bir hasta temas
sayfasi aynen soyluyor:

> *"Please click the 'Initiate contact' button below and fill in the form to request contact with the
> clinician for this patient. Your request will be validated by the DECIPHER team and forwarded to
> the clinician on your behalf. For this reason, it is important that you provide as much information
> as possible for reasons for contacting the clinician (collaboration, request for more information,
> etc.). You can also provide other DECIPHER identifiers in your request if you wish to contact more
> than one clinician."*

ve `/join` sayfasinda: *"It is not necessary to join DECIPHER to find and read open-access patient
data."* Birincil literaturde de var (Chatzimichali ve ark., *Hum Mutat* 2015, PMC4832335):
*"External matchmaking: external users can exploit the search and browsing functionality, and query
through publicly available anonymized patient data in order to discover potential matches and
initiate contact with the depositor via DECIPHER."* MME katilimci sayfasi: butun dis temas
talepleri iletilmeden once *"manually vetted for authenticity"*.

**KIF13A'nin DECIPHER'deki durumu (canli, 2026-09-16):** gen sayfasi *"DECIPHER holds no open-access
sequence variants in this gene"* diyor, ama **55 eslesen hasta varyanti** var: **47 kopya sayisi
degisimi**, 8 "other". Acik erisimli CNV kayitlari gezilebiliyor ve her birinde Contact sekmesi var.

**Yani DECIPHER bugun KIF13A icin gercek bir oku-ve-temas-kur kaynagidir**, hesap acmadan. ClinVar'da
bulunmayan bir aracilik mekanizmasi sunuyor. (JSON API ise kapali: butun `/api/` uc noktalari
HTTP 403 `"Client key not provided"`.)

### 6e. Eslesme ne kadar ise yariyor, gercek sayilar

- **Care4Rare Canada, 2 yil** (PMID **34906465**, *Genet Med* 2022, DOI 10.1016/j.gim.2021.08.014):
  **194 yeni aday gen** gonderildi, **1514 eslesme** cikti, **gonderilen genlerin %15'i isbirligine
  donustu**. Eslesmelerin **%82'si GeneMatcher ile** ve fenotip/kalitim bilgisi olmadigi icin
  **her eslesmede ayrica e-posta trafigi gerekti**. Yazarlarin sonucu: etkili ama **emek yogun**.
- **GeneMatcher, 2.5 yil** (PMID **30563986**, *Genet Med* 2019, DOI 10.1038/s41436-018-0383-z),
  Europe PMC ozetinden dogrulandi: *"We used GeneMatcher over 2.5 years to share 71 novel candidate
  genes identified by ES"*, *"Matches occurred in 60/71 candidate genes allowing to confirm the
  implication of 39% of matched genes as causal and to rule out 6% of them."*
  **Yani aday genlerin %85'i en az bir eslesme uretti, eslesenlerin %39'u nedensel dogrulandi.**
  *(Ayni makaleden sik alintilanan "ilk eslesmeye kadar medyan 4 saat", "medyan yanit suresi 31 saat"
  ve "38 gondericiden 9'u yanit vermedi" rakamlari odeme duvarinin arkasinda kaldi: **DOGRULANMADI**.
  Kullanilacaksa bu kayitla kullanilmali.)*
- **GeneMatcher'in kurumsal etkisi** (PMID **35224800**, *Hum Mutat* 2022): GeneDx tek basina tum
  GeneMatcher gonderimlerinin %20'sinden fazlasini yapiyor ve uc yilda 200'den fazla yeni
  gen-hastalik iliskisi yayini uretmis.

**Bizim icin cikarim:** KIF13A'yi kaydetmek buyuk olasilikla eslesme uretir (17,144 essiz genin
12,831'i eslesmis). Ama eslesmenin bilgiye donusmesi aylarca surer ve karsi tarafla insan-insana
yazismayi gerektirir. **Bu is bir insanin ustlenmesi gereken bir istir, bir scriptin degil.**

### 6f. Bize kapali dugumler

| Dugum | Karar | Gerekce (canli dogrulandi 2026-09-16) |
|---|---|---|
| **PhenomeCentral** | **HAYIR** | Hesaplar *"appropriate for clinicians and scientists who can contribute phenotype and genotype information about Patients ... under their care or being studied in their research environment"*. Kurumsal e-posta zorunlu: *"G-mail/Hotmail/Yahoo/other non-academic or non-institutional e-mail addresses etc. will not be accepted"*. Elle gecmis kontrolu: LinkedIn, hastane siteleri, PubMed. Ve acikca: *"Applicants will not be provided accounts for the purpose of browsing PhenomeCentral as the software is not designed for these types of searches"* |
| **MyGene2** | **HAYIR, ama yazili bir arka kapi var** | Katkici olmak icin *"A de-identified phenotype description and at least one candidate gene/variant"* gerekiyor. **Ama SSS aynen sunu diyor: *"Basic-science researchers without cases may contact us for alternatives."*** Ayrica klinisyen ve arastirmaci profilleri **public**, yani KIF13A kayitlari hesapsiz gezilebilir |
| **seqr** | Giris EVET, MME gonderimi HAYIR | Demo projesi icin AnVIL/Terra uzerinden self-servis hesap, ucretsiz (*"There are currently no costs associated with requests to load data from your AnVIL workspace to seqr"*). Ama MME gonderimi yuklenmis, joint-called VCF'li gercek bir proje gerektiriyor |
| **PubCaseFinder** | **HAYIR** | Kullanim sartlari: *"The PubCaseFinder is intended to be used by qualified and licensed physicians"*. CaseSharing kayitsiz calisiyor ama o HPO tabanli bir tani/vaka arama araci, gen park etme yeri degil |
| **Monarch Initiative** | **HAYIR, MME dugumu olu** | MME katilimci sayfasi: *"Given Monarch contains no patient data, no notifications of matches are sent."* GA4GH wiki'sinde listeli `mme.monarchinitiative.org` uc noktasi **NXDOMAIN** (test edildi 2026-09-16). Monarch, matchmakerexchange.org'un guncel resmi dugum listesinde **yok** |
| **PatientMatcher / Scout** | **HAYIR** | Kendi kendine barindirilan acik kaynak MME sunucu yazilimi; canli Stockholm dugumu hastalari yalniz Karolinska'da *"a registered and authorized collaborator"* uzerinden Scout ile aliyor. Halka acik kayit yok |
| **IRUD** | **HAYIR** | *"Only clinicians and researchers who join IRUD project can access and register patients' information"*. Japon AMED ulusal programi; Japonya disindan bir gruba acik kayit yolu bulunamadi. Kesin kurallar **DOGRULANMADI**, kapali varsayin |

### 6g. MME API'sini dogrudan sorgulayabilir miyiz? Hayir, test edildi

**Canli test (2026-09-16):** `X-Auth-Token` **olmadan** ve **bos JSON govdeyle** (`{}`, hicbir hasta
verisi yok) gonderilen POST istekleri:

| Uc nokta | Kimlik dogrulamasiz POST sonucu |
|---|---|
| `phenomecentral.org/rest/remoteMatcher/match` | **401** |
| `seqr.broadinstitute.org/api/matchmaker/v1/match` | **401** |
| `mygene2.gs.washington.edu/MyGene2/mme/match` | **401** |
| `mme.monarchinitiative.org/match` | baglanti kurulamadi, NXDOMAIN |

**Ulasilabilen her dugum kimlik dogrulamasiz sorguyu reddediyor. Halka acik sorgu uc noktasi yok.**

Spesifikasyon (`github.com/ga4gh/mme-apis/blob/master/search-api.md`, surum 1.1,
DOI 10.5281/zenodo.997379): *"This is a protected resource and as such requires the use of an API key
to access."* Hata tablosu: *"401 | Unauthorized | invalid API key"*. Istek
`<base_remote_url>/match` adresine POST, surum `Content-Type: application/vnd.ga4gh.matchmaker.<version>+json`
basliginda. Govde: `patient.id` **zorunlu**, `patient.contact` **zorunlu** (`name` ve `href` zorunlu;
`roles` degerleri `"clinician" | "researcher" | "patient"` -- **spesifikasyonun kendisi bir
`researcher` rolu tanimliyor**), ve **`features` veya `genomicFeatures`'dan en az biri zorunlu**.

Dugum olmak elle onaya bagli (`join-protocol.md`): *"Allowing an external site (A) to submit match
requests to a local site (B) requires manual acceptance and configuration of the A system into the B
system."* Token `X-Auth-Token` basliginda, 255 karakterden kisa, *"A random 40-chars SHA1 key is
recommended"*, PGP/GPG sifreli e-postayla takas ediliyor, *"For data security, HTTPS is mandatory"*.
GeneMatcher da ayni seyi soyluyor (`genematcher.org/help/mmapi`): *"You need to join the Matchmaker
Exchange project to be able to exchange match requests with GeneMatcher."*

**Dugum olmanin sartlari agir** (Matchmaker Exchange Service Requirements, **Mart 2026**,
matchmakerexchange.org'dan PDF olarak indirildi): 11 sart; aralarinda federe sorgu icin vaka verisi
yatirma zorunlulugu, en az iki nokta-nokta API baglantisi, Metrics API uygulamasi, GA4GH Data
Security Infrastructure Policy beyani ve bir Steering Committee temsilcisi.
**Bizim durumumuz icin kritik olan:** 1 numarali sart su istisnayi tasiyor:

> *"Require users to deposit case data to undertake a federated query across the MME service
> providers. (Note: Not required for MME Connected Knowledge Sources and MME Functional Study
> Nodes.)"*

Yani **Mart 2026 sartlari, hasta verisi tutmayan bir dugum sinifini resmen taniyor**
("MME Functional Study Nodes") ve **ModelMatcher tam olarak bu sinifta**. Teorik olarak bir
hesaplamali grup bu sinifta dugum olabilir, ama tek bir gen icin bu emege girmek mantiksiz.

**Sonuc:** Kabul edilmis dugum olmadan MME sorgulanamaz. Ama **ucretsiz bir GeneMatcher hesabi acip
KIF13A'yi gondererek ve gonderim sirasinda sorgulanacak MME veritabanlarini secerek tum federasyona
dolayli olarak ulasilabilir.**

## 7. Buyuk kohortlar: izin istemeden bakabilecegimiz yerler (yapildi)

Tablo: `data/interim/casefinding/KIF13A_cohort_portals_20260916.tsv`

**Alti kaynagin tamami bugun kayit olmadan, ucret odemeden sorgulandi.** Tek maliyet engeli
Genebass'in toplu indirmesi (requester-pays GCP kovasi). Asagidaki rakamlarin tamami canli API
cagrilariyla goruldu.

### 7a. gnomAD v4 kisitlama, KIF13A vs KIF1A

MANE Select ENST00000259711.11 (NM_022113), kisitlama bayragi yok.

| Olcut | KIF13A | KIF1A (kontrol) |
|---|---|---|
| pLI | **0.999750572233937** | 1 |
| LOEUF (`oe_lof_upper`) | **0.4489627711545366** (11. persentil) | 0.2587130941913371 (3. persentil) |
| `oe_lof` | 0.37233590809106626 | 0.20234013919264657 |
| gozlenen / beklenen LoF | 79 / 212.174 | 46 / 227.340 |
| `lof_z` | 7.756 | 10.203 |
| `mis_z` | **2.79488363965389** | 6.867936501880208 |
| `oe_mis` | **0.860569442734971** (GA 0.8288-0.8938) | 0.6770190021350824 |
| gozlenen / beklenen missense | 1901 / 2209.002 | 1683 / 2485.898 |
| `oe_syn` / `syn_z` | 0.9757 / 0.4117 | - |

**Durust okuma:** KIF13A **guclu LoF intoleransi** gosteriyor (pLI 0.9998, LOEUF 0.449, 11.
persentil) ama **missense kisitlamasi mutevazi** (mis_z 2.79, o/e 0.86) ve **her eksende KIF1A'dan
belirgin sekilde daha az kisitli**. Bu, gen duzeyinde bakildiginda missense-baskin bir modeli
kendiliginden desteklemez. Bizim hipotezimizin gercek testi **motor alana yerel kisitlamadir** ve bu
gen duzeyi rakamlari onu olcmez. Bu gerilimi raporda saklamamak gerekir.

> **Surum tuzagi (CLAUDE.md uyarisiyla ayni cinsten):** gnomAD GraphQL'inde
> `gene(reference_genome: GRCh38)` v4 kisitlamasini, `GRCh37` ise **v2.1.1** kisitlamasini donduruyor.
> GRCh37 degerleri: pLI 0.99999, LOEUF **0.232**, gozlenen/beklenen LoF 15/99.7, mis_z 2.047.
> **PMID 34303382'de gecen pLI 1.00 bu v2 degeridir.** API hicbir surum dizesi tasimadigi icin
> "v4.x" kesin, ama tam nokta surumu (4.1 mi) **DOGRULANMADI**. gnomAD tanim geregi fenotipsizdir:
> GraphQL `Gene` tipinde "phenotype" iceren hicbir alan yok.

### 7b. Genebass, 394,841 UKB eksomu, 4,529 fenotip

Acik, giris yok, onay tiklamasi yok. Per-gen degerler `main.genebass.org/api` uzerinden **okundu**.

| Yuk kumesi | Fenotip sayisi | En kucuk p | p<2.5e-7 | p<1e-6 | p<1e-5 |
|---|---|---|---|---|---|
| pLoF | 4,525 | 2.375e-04 | **0** | **0** | **0** |
| missense\|LC | 4,529 | 5.262e-05 | **0** | **0** | **0** |
| synonymous | 4,529 | 8.151e-05 | **0** | **0** | **0** |

**Uc anotasyon kumesinin hicbirinde eksom capinda anlamliliga yaklasan tek bir iliski yok.**

*Toplu indirme uyarisi:* Hail tablolari `gs://ukbb-exome-public/500k/results/` altinda ama kova
**requester-pays**; anonim HTTPS `HTTP 400 UserProjectMissing` donuyor. Kayit gerekmiyor ama kendi GCP
projeniz gerekiyor ve cikis trafigini odersiniz. **Yani "ucretsiz" degil.**

### 7c. AstraZeneca PheWAS Portal

Acik tarayici, giris yok, ucret yok; `/api` altindaki REST arka ucu kimlik dogrulamasiz okunabiliyor.
KIF13A yalnizca **500k WGS (v2)**, **460k WGS** ve **460k WES** veri setlerinde var, **470k (v5)**
eksom setinde yok.

| Veri seti | Donen sonuc | En kucuk p | p<1e-8 | p<2e-9 |
|---|---|---|---|---|
| 500k WGS v2 (484,111 genom, 15,129 fenotip, 11 collapsing modeli) | 21,887 | **2.576e-05** | **0** | **0** |
| 470k v5 (419,391 eksom, EUR) | 2,697 | 3.684e-05 | **0** | **0** |
| KIF1A kontrolu, 500k WGS | 18,776 | 1.603e-06 | 0 | 0 |

En yuksek KIF13A hiti: p=2.576e-05, `flexdmg`, surekli, *"25845 Volume of grey matter in Cuneal
Cortex (right)"*, **n=179**. Ikincisi p=2.654e-05, *"Mean OD in cingulum cingulate gyrus on FA
skeleton (left)"*, **n=44**, bes modelde tekrarlaniyor. **Ikisi de cok kucuk orneklemli beyin
goruntuleme fenotipi: gurultu sayilmali.**

> **Tuzak, kaydedilmeye deger:** birden fazla `collapsing_model_ids` virgulle birlestirilirse API
> sessizce `"Associations were not found"` donduruyor ve bu gercek bir negatif gibi okunuyor.
> **Her model kimligi tek tek sorgulanmali.**
>
> **Anlamlilik esigi:** portalda hicbir yerde gosterilmiyor. Yayinlanmis AZ esigi p < 1e-8
> (Wang ve ark., *Nature* 2021) ama portalin kendisinden **DOGRULANMADI**. Hem 1e-8 hem 2e-9
> uygulandi; **KIF13A ikisini de gecemiyor.** Toplu ozet istatistik indirme baglantisi bulunamadi
> (var olup olmadigi **DOGRULANMADI**).

### 7d. All of Us public Data Browser

Acik, hesap yok; `public.api.researchallofus.org` kimlik dogrulamasi istemiyor.
CDR **v0-6-rc5** (*"2025q4r5_combined"*), `numParticipants` **747,040**.
KIF13A gen olarak taniniyor; **varyant sayisi 127,713** (KIF1A 79,479). Ornek satir:
`6-17759115-C-G` AC=211, AN=1,071,234, AF=1.97e-4, hom=0. AN yaklasik 1,071,240 olmasi bu lokusta
kabaca **535,600 bireyin WGS'i** oldugunu ima ediyor.

**DOGRULANDI: public Data Browser hicbir genotip-fenotip iliskisi vermiyor.** Tum genomik API yuzeyi
varyant arama ve sayimdan olusuyor (`search-variants`, `variant-details`, `participant-counts`,
`genomic-filter-options`, ...); **association, burden veya PheWAS uc noktasi yok.** Sitenin kendi
tanimi: veriler *"aggregate-level"*, *"rounded aggregate data to counts of 20"*, ve
*"Individual-level data are available for analysis in the Researcher Workbench."*

**Registered/Controlled kademe sartlari (kurumsal baglanti, DUA, egitim, ABD ikametgahi) bu oturumda
kontrol edilmedi: DOGRULANMADI.**

### 7e. Open Targets Platform

623 iliskili hastalik var **ama veri tipi dagilimi aldatici**: `animal_model` 541, `literature` 43,
`genetic_association` 41, `rna_expression` 5. En yuksek genel iliski skoru yalnizca **0.286**
(glomerulonefrit).

**KRITIK NEGATIF:** KIF13A'ya bagli 525 credible set'in `studyType` dagilimi
`eqtl` 343, `tuqtl` 130, `sqtl` 49, `sceqtl` 3 ve **`GWAS` = 0**.
**KIF13A'ya atanmis tek bir GWAS credible set'i yok.**

L2G kaniti yalnizca 13 hastalik satirinda var ve **hicbiri norolojik degil**: glomerulonefrit 0.773
(PMID 39024449), meme adenozu 0.726, lagoftalmi 0.709, karaciger bozuklugu 0.705 (PMID 40069456),
multinodular guatr 0.698, serviks karsinomu 0.651, okuler hipotansiyon 0.651, artropati 0.581,
alkol tuketimi 0.570 (PMID 41606330), mevsimsel alerjik rinit 0.570, diyabetik ketoasidoz 0.549,
ozofagus bozuklugu 0.371.

Tek norolojik genetik iliski **popülasyon kaniti degil, ClinVar kaynakli**: **RCV003311600 /
rs2532452593**, skor 0.32, apraksi + konusma apraksisi + nobet + orta zihinsel yetersizlik + Chiari
malformasyonu. **Bu bizim 4f'deki Lead 1'in ta kendisi (VCV002570672).** Open Targets'in tum
norolojik sinyalinin tek kaynagi o bir kayittir.

> **Uyari:** listedeki 181 "noro" hastaliginin ezici cogunlugu **fare goz/retina `animal_model`**
> kayitlaridir (Leber konjenital amaurozu 0.061, retinitis pigmentosa 0.056). **Bunlari insan
> genetik kaniti gibi sunmayin.**

### 7f. GWAS Catalog

Acik, anahtar yok. `/genes/{symbol}` 404 donuyor;
`/singleNucleotidePolymorphisms/search/findByGene?geneName=KIF13A` kullanilmali. HTTP 429 ile hiz
siniri var.

> **KRITIK TUZAK:** `findByGene` 93 essiz SNP donduruyor ama **yalnizca 41'inin KIF13A'ya uzakligi
> 0**. Diger 52'si "en yakin gen" veya "100 kb icinde" eslemeleri ve bircogu chr6:18.10-18.12 Mb'de,
> yani **KIF13A govdesinin (17,759,183-17,987,635) DISINDA ve TPMT bolgesinde**. Manset niteligindeki
> "yaslanma" ve tiyopurin hitleri (p degerleri 4e-54'e kadar) **o yan bolgeye aittir, KIF13A'ya
> degil.** Bunlari KIF13A sonucu diye bildiren herkes yanlis yapiyor.

Gen govdesi icindeki gercek iliskiler (p'ye gore):

| Ozellik | SNP | p | PMID |
|---|---|---|---|
| TPMT protein duzeyi | rs562345710 | 4e-20 | 39789286 |
| SCARA5 protein duzeyi | rs112683346 | 5e-15 | 39789286 |
| TPMT protein duzeyi | rs748954374 | 5e-14 | 39789286 |
| Dolasimdaki SCARA5 | rs570304138 | 1e-13 | 42097137 |
| Dolasimdaki SCARA5 | rs536447969 | 2e-12 | 42097137 |
| Akciger fonksiyonu FEV1/FVC | rs1014146 | 2e-11 | 30595370 |
| Aort stenozu | rs6942390 | 6e-10 (OR 1.05) | 41419686 |
| **Sulkus derinligi** | **rs75428660** | **1e-09** (beta 6.09 z, N=33,748) | **34910505** |
| KOAH | rs62398171 | 2e-09 | 37069358 |
| **Kortikal kalinlik** | **rs75428660** | **5e-09** (beta 5.85 z) | **34910505** |
| **Kortikal kalinlik** | **rs75428660** | **2e-08** (N=35,657) | **34560273** |
| BMI-duzeltilmis bel cevresi | rs75964046 | 3e-08 | 31669095 |
| Boy | rs3734232 (splice_region) | 3e-08 | 36224396 |
| *(esik alti)* **Progresif supranukleer palsi** | **rs4712314** | **2e-07** (beta 0.27; 1,718 vaka / 2,944 kontrol) | **39152475** |

**Tek genom capinda anlamli norolojik hikaye rs75428660'dir:** KIF13A'nin intronunda,
chr6:17,885,814, uzaklik 0 dogrulandi, iki bagimsiz beyin goruntuleme GWAS'inda kortikal kalinlik ve
sulkus derinligi ile 1e-09 ile 2e-08 arasinda.

**Sarti soylemek zorunlu:** GWAS Catalog **gen eslemesi** bildirir, nedensel atama yapmaz; ve
Open Targets KIF13A'ya **hicbir GWAS credible set'i atamiyor** (7e). Dolayisiyla rs75428660 sinyali
**konumsaldir ve L2G destekli degildir.**

### 7g. Bu bolumun ozeti

Alti kaynagin tamami bugun, kayitsiz, kullanilabilir durumda. Durust manset **kisitlama ile iliski
arasindaki bosluktur**: KIF13A guclu LoF kisitlamasi gosteriyor (pLI 0.9998, LOEUF 0.449, 11.
persentil) **ama ne 394k Genebass eksomunda ne 484k AZ genomunda eksom capinda anlamliliga yaklasan
hicbir nadir varyant yuku sinyali yok**, ve tek genom capinda anlamli norolojik iliskisi
rs75428660'daki yaygin varyant beyin morfolojisi.

**Onemli sinir:** bu portallarin hicbiri bizim spesifik missense adaylarimizi test edemez. Yuk
testleri ya pLoF'lari ya tum missense'leri toplar; 52 ultra-nadir motor-alan missense'i UKB'de
yeterli tasiyici bulamaz. **Buradan cikan negatif sonuc hipotezimizi curutmez**, sadece "kaba gen
duzeyinde acik bir sinyal yok" der. Ustelik gen duzeyi missense kisitlamasinin mutevazi olmasi
(7a) tam da bunu bekletir: sinyal varsa **motor alana yerel** olmali.

## 8. Onerilen sira

1. **EN ONEMLI, klinisyen/PI imzasi gerektirir:** 4f'deki iki fenotipli vakanin gondericisine
   yazmak. **HudsonAlpha'yi once deneyin**: arastirma grubu, tek bir KIF13A gonderimi yapmislar,
   fenotip zengin (apraksi, nobet, orta ID, Chiari), ve ClinVar'da uc ayri iletisim kisisi yayinli.
   Ardindan **ISCA site 1 / Erin Riggs**: yalniz KIF13A'yi silen delesyon, etkilenmis cocuk,
   2017'den beri kriter uygulanmamis VUS olarak duruyor. Taslak: bolum 9b.
2. **Bugun, izinsiz, ucretsiz:** bolum 7'deki alti portal sorgulandi ve rakamlari bu belgeye islendi.
   Yeni yapilacak is yok; tekrar bakilacaksa AZ PheWAS'in model-tek-tek tuzagina dikkat.
3. **Insan karari, en dusuk surtunme:** **GeneMatcher** hesabi, "Researcher" rolu, yalin `KIF13A`
   gonderimi (istege bagli olarak bolum 3a'daki 12 varyantla). Gonderim ekraninda sorgulanacak MME
   veritabanlarini secerek DECIPHER, PhenomeCentral, MyGene2, seqr, PatientMatcher, RD-Connect ve
   IRUD'a tek hamlede ulasilir. Taslak: bolum 9a. **EULA celiskisini (6a) once okuyun.**
4. **Insan karari, ayni gun:** **ModelMatcher**'da hem Scientist Registry hem Match Dashboard hesabi;
   KIF13A'yi Primary gen olarak Public isaretleyin. "in silico" acikca uygunluk kapsaminda.
5. **Hesapsiz, hemen:** **DECIPHER**'de KIF13A'yi kapsayan 47 acik erisimli CNV kaydini gezip,
   norogelisimsel fenotipli olanlar icin "Initiate contact" formunu kullanmak. DECIPHER ekibi
   dogrulayip klinisyene iletiyor.
6. **2 ile 4 hafta, denetimli:** **RD-Connect GPAP** basvurusu (kimlik taramasi, Code of Conduct,
   DAC karari). Uc acik kapinin en kaliteli verisi burada.
7. **Ikincil:** **Ambry**'ye G283V ve H175R icin isbirligi teklifi (taslak 9c).
   Toplu VUS yatirimi oldugu icin R0'dan sonra gelir. **MyGene2**'ye kendi yazili arka kapisindan
   yazmak (*"Basic-science researchers without cases may contact us for alternatives"*).
8. **Uzun vadeli ve en cok kontrolumuzde:** in silico hipotez makalesini yayimlamak. KIF13A'yi
   norogelisimsel aday olarak haritaya koydugumuz an, bir sonraki klinisyen varyanti gordugunde bizi
   bulur. **Bugun literaturde bizim rolumuzu oynayacak hicbir sey yok:** 2026 tarihli kinesin-NDD
   derlemesi KIF13A'yi hastalik genleri arasinda saymiyor (5d).

**Yapilmayacaklar:**
- **Rogan Lab'a yazmayin.** Orada hasta yok, TCGA tumor ornekleri var (4b).
- **PhenomeCentral'a basvurmayin.** Kendi politikasi gezinme amacli hesap vermedigini yaziyor (6f).
- **MME dugumu olmaya calismayin.** Tek bir gen icin 11 sartlik bir surec mantiksiz (6g).
- **GWAS Catalog'un TPMT bolgesi hitlerini KIF13A sonucu diye bildirmeyin** (7f).
- **Open Targets'in fare retina `animal_model` satirlarini insan kaniti gibi sunmayin** (7e).

## 9. Taslak metinler

> ### UYARI
> Asagidaki metinlerin hicbiri gonderilmedi. Bunlar **bir insanin (danisman klinisyen veya
> laboratuvar sorumlusu) kendi adiyla, kendi kurum imzasiyla gozden gecirip gonderecegi taslaklardir.**
> Koseli parantezli alanlarin doldurulmasi gerekir. Hicbir metin bir hasta, bir tani veya bir
> patojenisite iddiasi ileri surmuyor; hepsi acikca hesaplamali bir hipotez oldugunu soyluyor.

### 9a. Eslestirme platformu "ilgi beyani" (kisa, GeneMatcher / ModelMatcher gen notu icin)

**Turkce:**

> KIF13A (HGNC:14566, MANE NM_022113.6) genine hesaplamali bir ilgi duyuyoruz. Paralogu KIF1A'nin
> klinik varyantlari uzerinde kalibre ettigimiz bir onceliklendirme ile, KIF13A kinezin motor alaninda
> KIF1A'da patojenik varyantlarin bulundugu rezidulere denk gelen ve gnomAD'da gorulmeyen bir dizi
> missense aday belirledik. Elimizde hasta yok. KIF13A'da nadir missense varyant, ozellikle motor
> alanda (yaklasik aminoasit 1-370) de novo bir varyant goren ve norogelisimsel bir fenotip
> (gelisim geriligi, zihinsel yetersizlik, spastisite, ataksi, noropati) bulunan klinisyenlerle
> iletisim kurmak istiyoruz. Yapisal ve in silico degerlendirmeyi karsiliksiz saglayabiliriz.
> Bugun itibariyla KIF13A'nin kurulmus bir Mendelyan hastalik iliskisi yoktur; bilinen tek
> fenotiplenmis vaka postzigotik inaktive edici bir varyanttir (PMID 33526817).
> Iletisim: [ad, unvan, kurum, e-posta].

**English:**

> We have a computational interest in KIF13A (HGNC:14566, MANE NM_022113.6). Using a prioritisation
> calibrated on clinical variants of its paralogue KIF1A, we have identified a set of missense
> candidates in the KIF13A kinesin motor domain that fall on residues where KIF1A carries pathogenic
> variants and that are absent from gnomAD. We do not hold any patients. We would like to hear from
> clinicians who have seen a rare missense variant in KIF13A, particularly a de novo variant in the
> motor domain (approximately residues 1-370), in an individual with a neurodevelopmental phenotype
> (developmental delay, intellectual disability, spasticity, ataxia, neuropathy). We can offer
> structural and in silico assessment at no cost. As of today KIF13A has no established Mendelian
> disease association; the single phenotyped case in the literature is a postzygotic inactivating
> variant (PMID 33526817).
> Contact: [name, title, institution, email].

### 9b. EN ONEMLI TASLAK: fenotipli vakanin gondericisine (HudsonAlpha)

> **Gonderen bir klinisyen veya laboratuvar sorumlusu olmalidir.** Hasta kimligi istemiyor, isbirligi
> teklif ediyor. Alici: ClinVar gonderici sayfasinda yayinlanan HudsonAlpha iletisim adresleri
> (Susan Hiatt `(adres ClinVar başvuru sahibi sayfasında)`, PI Greg Cooper `(adres ClinVar başvuru sahibi sayfasında)`).

**Turkce taslak:**

> Konu: KIF13A p.Lys1122Glu (VCV002570672) uzerine isbirligi teklifi
>
> Sayin Dr. Cooper ve Dr. Hiatt,
>
> [Kurum]'da KIF13A (HGNC:14566) icin paralog-bilgili bir hesaplamali varyant onceliklendirme
> calismasi yurutuyoruz. Yaklasimimiz, paralog KIF1A'nin ClinVar'daki klinik varyantlari uzerinde
> kalibre edilmis olup KIF13A kinezin motor alanindaki missense varyantlari siralamaktadir.
>
> Sistematik bir taramada, KIF13A'nin 293 tek-gen ClinVar kaydi icinde **adlandirilmis bir norolojik
> fenotip tasiyan tek kaydin sizin gonderiminiz oldugunu** gorduk:
> NM_022113.6:c.3364A>G (p.Lys1122Glu), VCV002570672 / SCV004009653, apraksi, konusma apraksisi,
> nobet, orta duzey zihinsel yetersizlik ve Chiari malformasyonu ile birlikte. Kayitlarin geri
> kalaninin 217'si "not specified", 13'u "not provided" durumundadir.
>
> KIF13A'nin bugun kurulmus bir Mendelyan hastalik iliskisi yoktur: ClinGen, Gene2Phenotype, PanelApp
> ve OMIM fenotip kayitlarinda yer almaz; yalnizca SysNDD'de tek bir yayina dayanan "Limited"
> siniflandirmasi vardir (PMID 33526817, postzigotik inaktive edici varyant, hipomelanozis of Ito ve
> gelisim geriligi). Buna ek olarak ClinVar'da, yalnizca KIF13A'yi kapsayan yaklasik 253 kb'lik bir
> delesyonun (VCV000442594) konusma ve dil gelisiminde gecikme ile birlikte bildirildigini tespit
> ettik. Yani KIF13A icin dagilmis halde duran, henuz kimsenin bir araya getirmedigi bir insan kaniti
> kumesi var.
>
> Su konularda isbirligine acik olup olmadiginizi sormak istiyoruz: bu probandin fenotipinin ve
> varyant kalitim durumunun (de novo olup olmadigi dahil) ortak bir vaka serisi cercevesinde
> degerlendirilmesi; ve sizin gonderdiginiz varyant ile bizim motor-alan adaylarimizin ayni
> mekanizmaya isaret edip etmediginin incelenmesi. Tarafimizdan karsiliksiz olarak yapisal modelleme,
> paralog kalibrasyonu, kisitlama analizi ve in silico degerlendirme sunabiliriz.
>
> Elimizde hasta yoktur; katkimiz tamamen hesaplamalidir. Ilgilenirseniz aday listemizi ve yontem
> ozetimizi memnuniyetle paylasiriz.
>
> Saygilarimla,
> [Ad, unvan, kurum, kurumsal e-posta, telefon]

**English draft:**

> Subject: Collaboration on KIF13A p.Lys1122Glu (VCV002570672)
>
> Dear Dr. Cooper and Dr. Hiatt,
>
> We run a paralogue-informed computational variant prioritisation study on KIF13A (HGNC:14566) at
> [institution]. Our approach is calibrated on ClinVar clinical variants of the paralogue KIF1A and
> ranks missense variants in the KIF13A kinesin motor domain.
>
> In a systematic sweep of all 293 single-gene ClinVar records for KIF13A, we found that **yours is
> the only record carrying a named neurological phenotype**: NM_022113.6:c.3364A>G (p.Lys1122Glu),
> VCV002570672 / SCV004009653, reported with apraxia, speech apraxia, seizure, moderate intellectual
> disability and Chiari malformation. Of the remainder, 217 are "not specified" and 13 "not provided".
>
> KIF13A has no established Mendelian disease association today: it is absent from ClinGen
> gene-disease validity, Gene2Phenotype, PanelApp and OMIM phenotype records, and appears only in
> SysNDD as "Limited" on the strength of a single publication (PMID 33526817, a postzygotic
> inactivating variant with hypomelanosis of Ito and developmental delay). Separately, we identified a
> ClinVar record for an approximately 253 kb deletion containing KIF13A as its only protein-coding
> gene (VCV000442594), reported with delayed speech and language development. There is therefore a
> scattered body of human evidence for KIF13A that nobody has yet assembled.
>
> We would like to ask whether you would be open to collaborating on: assessment of this proband's
> phenotype and the inheritance status of the variant (including whether it is de novo) within a
> co-authored case series; and whether your variant and our motor-domain candidates point to the same
> mechanism. We can contribute structural modelling, paralogue calibration, constraint analysis and in
> silico assessment at no cost.
>
> We hold no patients; our contribution is entirely computational. We would be glad to share our
> candidate list and a methods summary if that is useful.
>
> Kind regards,
> [Name, title, institution, institutional email, phone]

**ISCA site 1 / Erin Riggs icin kisa varyant (odak delesyon):**

> Konu: KIF13A tek-gen delesyonu (VCV000442594) hakkinda soru
>
> Sayin Dr. Riggs,
>
> KIF13A (HGNC:14566) uzerine hesaplamali bir calisma yurutuyoruz. ClinVar'da ISCA site 1 tarafindan
> bildirilen `GRCh37 6p22.3(chr6:17758616-18011991)x1` kaydini (VCV000442594, SCV000585270)
> inceledik. Ensembl GRCh37 ile dogruladigimiz kadariyla bu aralikta **tek protein kodlayan gen
> KIF13A'dir** ve kayit, konusma ve dil gelisiminde gecikme olan etkilenmis bir bireyde bildirilmis,
> kriter uygulanmadan VUS olarak siniflandirilmistir.
>
> KIF13A'nin kurulmus bir hastalik iliskisi yoktur ve bu kayit bildigimiz kadariyla **izole KIF13A
> haploinsufficiency'si icin var olan en temiz insan kanitidir.** Bu delesyonun kalitim durumu
> (de novo olup olmadigi) veya daha ayrintili fenotipi hakkinda paylasilabilir bilgi olup olmadigini,
> ya da bu kaydin ClinGen doz duyarliligi kurasyonuna dahil edilmesinin uygun gorulup gorulmeyecegini
> sormak istiyoruz. Karsiliginda kisitlama analizi, paralog kalibrasyonu ve yapisal degerlendirme
> sunabiliriz.
>
> Saygilarimla,
> [Ad, unvan, kurum, kurumsal e-posta]

### 9c. ClinVar gonderici laboratuvarina isbirligi mektubu (Ambry, G283V ve H175R icin)

> **Gonderen bir klinisyen veya laboratuvar sorumlusu olmalidir.** Hasta duzeyinde bilgi talep
> etmiyor; isbirligine davet ediyor. Alici: ClinVar gonderici sayfasinda yayinlanan Variant HelpDesk
> adresi.

**Turkce taslak:**

> Konu: KIF13A VUS'lari uzerine arastirma isbirligi (VCV004092057, VCV003533722)
>
> Sayin Variant Science ekibi,
>
> [Kurum]'da KIF13A (HGNC:14566) uzerine hesaplamali bir varyant onceliklendirme calismasi
> yurutuyoruz. Yaklasimimiz, paralog KIF1A'nin ClinVar'daki klinik varyantlari uzerinde kalibre
> edilmis olup KIF13A kinezin motor alanindaki missense varyantlari siralamaktadir.
>
> ClinVar'da laboratuvariniz tarafindan bildirilen iki KIF13A varyanti bizim en ust siramizda yer
> aliyor:
> - NM_022113.6:c.848G>T (p.Gly283Val), VCV004092057 / SCV006347025
> - NM_022113.6:c.524A>G (p.His175Arg), VCV003533722 / SCV005608632
>
> Her ikisi de motor alanda ve KIF1A'da patojenik varyantlarin bulundugu rezidulere denk geliyor.
> KIF13A'nin bugun kurulmus bir Mendelyan hastalik iliskisi olmadiginin ve bu varyantlarin VUS olarak
> siniflandirildiginin farkindayiz.
>
> Hasta duzeyinde hicbir bilgi talep etmiyoruz. Sormak istedigimiz sudur: bu varyantlarin bildirildigi
> vakalarda norogelisimsel veya norolojik bir endikasyon bulunup bulunmadigini, uygun bir cerceve
> icinde (ornegin ortak yazarlikli bir vaka serisi, veya laboratuvarinizin uygun gordugu herhangi bir
> veri paylasim mekanizmasi) degerlendirmeye acik olur musunuz? Yapisal modelleme, paralog
> kalibrasyonu ve in silico degerlendirmeyi karsiliksiz sunabiliriz.
>
> Ilgilenirseniz, KIF13A gen listemizi ve yontem ozetimizi memnuniyetle paylasiriz.
>
> Saygilarimla,
> [Ad, unvan, kurum, kurumsal e-posta, telefon]

**English draft:**

> Subject: Research collaboration on KIF13A variants of uncertain significance (VCV004092057,
> VCV003533722)
>
> Dear Variant Science team,
>
> We are running a computational variant prioritisation study on KIF13A (HGNC:14566) at
> [institution]. Our approach is calibrated on ClinVar clinical variants of the paralogue KIF1A and
> ranks missense variants in the KIF13A kinesin motor domain.
>
> Two KIF13A variants submitted by your laboratory rank at the top of our list:
> - NM_022113.6:c.848G>T (p.Gly283Val), VCV004092057 / SCV006347025
> - NM_022113.6:c.524A>G (p.His175Arg), VCV003533722 / SCV005608632
>
> Both lie in the motor domain at residues where KIF1A carries pathogenic variants. We are aware that
> KIF13A has no established Mendelian disease association at present and that both variants are
> classified as VUS.
>
> We are not requesting any patient-level information. Our question is whether you would be open to
> considering, within whatever framework you deem appropriate (for example a co-authored case series,
> or any data sharing mechanism your laboratory supports), whether the cases in which these variants
> were reported carried a neurodevelopmental or neurological indication. We can contribute structural
> modelling, paralogue calibration and in silico assessment at no cost.
>
> We would be glad to share our KIF13A candidate list and a methods summary if that is useful.
>
> Kind regards,
> [Name, title, institution, institutional email, phone]

### 9d. Kisa "gen ilgisi" duyurusu (SysNDD kuratorlerine veya bir konsorsiyuma)

**Turkce:**

> KIF13A SysNDD'de "Limited" kategorisinde (entity 3406, MONDO:0001071, otozomal dominant, 2021),
> tek bir yayina dayali (PMID 33526817, postzigotik inaktive edici varyant). Paralog-bilgili bir
> hesaplamali onceliklendirme ile KIF13A motor alaninda, KIF1A'nin patojenik rezidulerine denk gelen
> ve gnomAD'da bulunmayan missense adaylari belirledik. Bu kategoriyi ilerletebilecek yeni vaka
> bildirimleri veya devam eden kurasyonlardan haberdar olmak isteriz; karsiliginda aday listemizi ve
> yapisal degerlendirmeyi paylasabiliriz. Iletisim: [ad, kurum, e-posta].

**English:**

> KIF13A is curated in SysNDD as "Limited" (entity 3406, MONDO:0001071, autosomal dominant, 2021),
> resting on a single publication (PMID 33526817, postzygotic inactivating variant). Using a
> paralogue-informed computational prioritisation we have identified missense candidates in the
> KIF13A motor domain that fall on residues where KIF1A carries pathogenic variants and that are
> absent from gnomAD. We would be glad to hear of any new case reports or ongoing curation that could
> advance this classification, and can share our candidate list and structural assessment in return.
> Contact: [name, institution, email].

---

## 10. Bu belgenin sinirlari ve dogrulanmayanlar

**Dogrulanamayanlar:**
- **gnomAD kisitlamasinin tam nokta surumu** (4.1 mi, 4.x mi): API hicbir surum dizesi tasimiyor.
- **AZ PheWAS'in portal icinde belirtilen calisma capinda anlamlilik esigi**: portalda hicbir yerde
  gosterilmiyor. p<1e-8 yayinlanmis degerdir (Wang ve ark. *Nature* 2021) ama portaldan okunmadi.
  AZ'nin toplu ozet istatistik indirmesi olup olmadigi da bilinmiyor.
- **All of Us Registered/Controlled kademe sartlari** (kurumsal baglanti, DUA, egitim, ABD ikametgahi):
  basvuru sayfalari acilmadi.
- **KIF13A'nin DECIPHER'deki hasta sayisinin tam dokumu**: gen sayfasi 55 varyant (47 CNV, 8 other)
  gosteriyor ama bunlarin kac ayri hastaya karsilik geldigi ve fenotiplerin dagilimi okunmadi.
- **IRUD'un kesin uygunluk kurallari**: Japonya disindan bir gruba acik kayit yolu bulunamadi.
- **ModelMatcher Scientist Registry'nin ucretsiz oldugu**: ucret ibaresinin yoklugundan CIKARIM,
  acik bir alinti degil.
- **RD-Connect GPAP'ta verisiz kullanicinin MME'ye gonderim yapamayacagi**: CIKARIM, acik alinti degil.
- **OMIM'in tam morbidmap/genemap2 dosyalari** (lisans gerekiyor). `mim2gene.txt`'te KIF13A'nin tip
  alaninin `gene` olmasi guclu bir isaret ama tam kanit degildir.
- **PMID 30563986'nin "medyan 4 saat / 31 saat / 38'de 9 yanitsiz" rakamlari**: odeme duvarinin
  arkasinda. Kullanilacaksa bu kayitla kullanilmali.
- **Ambry'nin hasta duzeyinde bilgi paylasmayacagi beklentisi** makuldur ama yayinlanmis bir Ambry
  politikasindan alinti degildir.

**Yontemsel notlar:**
- Her kurumsal kaynak sorgusunda **KIF1A pozitif kontrol** olarak kullanildi, boylece bos donen
  sonuclarin API hatasi degil gercek yokluk oldugu gosterildi.
- genematcher.org, variantmatcher.org ve phenodb.org WebFetch'e Cloudflare 403 donuyor; **tarayici
  User-Agent'i ile curl HTTP 200 donduruyor.** GeneMatcher bulgularinin tamami canli sitedendir,
  arsivden degil.
- DECIPHER'in `/api/` uc noktalari Client Key olmadan HTTP 403 donuyor; gen ve hasta sayfalari ise
  tarayici ile okunabiliyor.

**Bu oturumda yapilmayanlar:** hicbir hesap acilmadi, hicbir form gonderilmedi, hicbir yere kayit
yapilmadi, hicbir e-posta gonderilmedi, hicbir DECIPHER "Initiate contact" butonuna basilmadi.
MME uc noktalarina yapilan testler **bos JSON govdeyle** yapildi; **hicbir hasta verisi
gonderilmedi.** Mevcut hicbir proje dosyasi (CLAUDE.md, progress.md, scripts/) degistirilmedi.

## 11. Uretilen veri dosyalari

Hepsi `atlas_mcp.core.save_tsv` ile yazildi, yaninda `.prov.json` provenance dosyasi var ve
`data/raw/PROVENANCE.jsonl` defterine islendi.

| Dosya | Icerik |
|---|---|
| `data/interim/casefinding/KIF13A_phenotyped_leads_20260916.tsv` | **En onemli dosya.** Norolojik fenotipi olan iki KIF13A kaydi, gonderici iletisim bilgileriyle |
| `data/interim/casefinding/KIF13A_curated_resources_20260916.tsv` | 9 kurumsal kaynagin canli durumu, KIF1A pozitif kontrolu ile |
| `data/interim/casefinding/KIF13A_clinvar_tierA_submitters_20260916.tsv` | ClinVar'daki 4 tier A adayi, tam gonderici provenansi ve TCGA tuzagi |
| `data/interim/casefinding/KIF13A_clinvar_cnv_overlaps_20260916.tsv` | KIF13A ile ortusen 20 tek-gen-disi ClinVar kaydi + ozet sayimlar |
| `data/interim/casefinding/KIF13A_literature_human_cases_20260916.tsv` | Yayinlanmis insan gozlemleri + tier A negatif sonucu |
| `data/interim/casefinding/matchmaking_nodes_20260916.tsv` | 12 eslestirme dugumunun uygunluk, ucret, denetim ve bildirim kosullari |
| `data/interim/casefinding/KIF13A_cohort_portals_20260916.tsv` | 6 kohort/iliski portali: erisim kosullari ve canli KIF13A gen duzeyi sonuclari |
