# VCV000442594 / SCV000585270 dogrulama raporu

6p22.3 253 kb delesyonunun "KIF13A haploinsufficiency icin en guclu insan kaniti" olarak kullanilabilirligi.

Tarih: 2026-09-16. Tum sayilar bu tarihte canli olarak cekildi. Hicbir PMID, erisim numarasi veya URL canli dogrulanmadan yazilmadi.

---

## Karar (verdict)

**Bu kayit KIF13A haploinsufficiency'sini desteklemiyor; olcek uzerinde en fazla "suggestive" (telkin edici), yani verilen skalanin en alt basamagi.** Kaydin kendisi teknik olarak dogru: ClinVar XML'i iddiada belirtilen her alani teyit ediyor (253.376 bp, GRCh37 chr6:17.758.616-18.011.991, tek gonderim, ISCA site 1 / OrgID 505237, "Uncertain significance, no assertion criteria provided", etkilenmis cocuk, HP:0000750). Ancak kanit degeri uc bagimsiz noktadan cokuyor. Birincisi ve belirleyicisi: gnomAD v4 ekzom CNV veri setinde (464.297 birey) **436229__DEL**, chr6:17.760.033-18.208.304 (448 kb), ClinVar araliginin **%99,3'unu** ve KIF13A gen govdesinin **%99,2'sini** iceriyor; yani hasta olarak siniflanmamis bir bireyde bu vakadan daha buyuk ve daha yikici bir delesyon zaten mevcut. Buna gnomAD v4 genom SV setindeki (63.046 birey) PASS filtreli **DEL_chr6_33713a79** (164 kb, genin %71,3'u, "lof") ve 84 kb'lik **436227__DEL** ekleniyor. Ikincisi: gozlem tamamen kalitim bilgisizdir. ClinVar XML'inde `<Origin>unknown</Origin>` yazar; de novo/maternal/paternal ayrimi, aile verisi, yas, cinsiyet, test edilen birey sayisi, es-varyant bilgisi ve gercek bir gonderen yorumu **yoktur**; gonderilen tek "yorum" tum ISCA gonderimlerinde ayni olan bir sablon metindir. De novo teyidi olmayan tek bir VUS, ACMG/ClinGen CNV cerceveisinde neredeyse hicbir puan tasimaz. Ucuncusu: ClinGen bugune kadar KIF13A'yi doz duyarliligi acisindan **hic kurate etmemistir** (1.540 genlik listede yoktur) ve 6p22.3 icin kurate edilmis bir bolge yoktur; nedensellik icin gosterilen Kaminsky 2011 / nstd101 baglantisi da guvenilir degildir, cunku nstd101 dbVar tarafindan 2021-11-03'te **obsolete** ilan edilmistir ve zaten "Case-Set", yani tasarimi geregi kontrol kolu olmayan bir calismadir. Ustune, delesyon KIF13A'dan baska 37 Ensembl enhancer'i ve 9 CTCF baglanma bolgesini de cikarmaktadir ve kirilma noktasi, KIF13A'dan daha kisitli bir genin (NUP153, pLI ~1,0, LOEUF 0,277) promotorune yalnizca 52 kb uzaktadir. Lehte kalan tek somut sey KIF13A'nin kendi kisitlilik profilidir (pLI 0,9998, LOEUF 0,449) ve bu, delesyonun degil genin ozelligidir.

**Guncelleme (ayni gun, EK BOLUM sonrasi):** DECIPHER ve DDD sorgulandiktan sonra karar seviyesi degismedi (hala "suggestive"), ama gerekce agirlasti. DECIPHER'da KIF13A'nin tek gen oldugu ikinci bir hasta delesyonu bulundu (hasta 451718, 146 kb, gen ici, VUS, "Autistic behavior"), yani bu ClinVar kaydi artik turunun tek ornegi degil; buna karsilik ayni veri setinde genin %73,0'unu silen ucuncu bir KIF13A-tek-gen delesyonu **anneden kalitilmis ve "likely benign"** olarak siniflanmistir (hasta 504413), ve DDD'nin KIF13A'daki iki "Likely LOF" varyantinin (bir stop_gained, bir frameshift) **ikisi de maternal, ikisi de "Uncertain"**dir. KIF13A **DDG2P listesinde yoktur** ve DECIPHER'in %HI tahmini geni **%42,3** persentiline koyar. Dolayisiyla ilk karardaki "kalitim bilinmiyor" maliyeti artik yalnizca bir bosluk degil, aktif olarak ters yonde bir sinyaldir. Ayrintilar icin EK BOLUM'e bakiniz.

---

## 1. Bu delesyon gorunuste saglikli insanlarda goruluyor mu?

**Evet, ve bu tek basina iddiayi ciddi olarak zayiflatiyor.**

Cekilen kaynak: `POST https://gnomad.broadinstitute.org/api` (GraphQL), `dataset: gnomad_sv_r4` ve `dataset: gnomad_cnv_r4`, bolge sorgusu GRCh38 chr6:17.700.000-18.060.000.
Tarayici karsiligi: https://gnomad.broadinstitute.org/region/6-17700000-18060000?dataset=gnomad_sv_r4

Referans koordinatlar:
- KIF13A GRCh38: chr6:17.758.198-17.988.274 (230.076 bp)
- ClinVar delesyonu GRCh38 (dbVar yerlesimi): chr6:17.758.385-18.011.760 (253.375 bp)

### 1a. gnomAD v4 ekzom CNV (gnomad_cnv_r4, SN = 464.297 birey)

Bolgede 11 CNV cagrisi var; KIF13A ile ortusen delesyonlar:

| Cagri | Koordinat (GRCh38) | Uzunluk | SC | SN | SF | Genin % kadari | ClinVar araligi ile ortusme |
|---|---|---|---|---|---|---|---|
| **436229__DEL** | 17.760.033-18.208.304 | **448.271 bp** | **1** | 464.297 | **2,15e-06** | **%99,2** | **251.727 bp = ClinVar delesyonunun %99,3'u** |
| 436227__DEL | 17.765.150-17.849.244 | 84.094 bp | 1 | 464.297 | 2,15e-06 | %36,6 | 84.094 bp |
| 436318__DEL | 17.825.671-17.837.681 | 12.010 bp | 2 | 464.297 | 4,31e-06 | %5,2 | 12.010 bp |
| 436327__DEL | 17.833.864-17.837.681 | 3.817 bp | 1 | 464.295 | 2,15e-06 | %1,7 | 3.817 bp |

Ayni bolgede kazanclar da var ve bazilari sik: 436092__DUP (255.971 bp, SC=67, SF=1,44e-04) ve 436227__DUP (61.970 bp, SC=44, SF=9,48e-05). Yani KIF13A lokusu kopya sayisi acisindan sessiz bir bolge degildir.

**436229__DEL kritik bulgudur.** ClinVar vakasinin delesyonunun %99,3'unu icine alir, 196 kb daha uzaga uzanir ve KIF13A'nin tum kodlayan bolgesini fiilen kaldirir (genin genomik ayak izinin yalnizca 3' ucundaki 1.835 bp disarida kalir). Bu birey siddetli pediatrik hastalik icin secilmemistir.

### 1b. gnomAD v4 genom SV (gnomad_sv_r4, AN = 126.092, yani 63.046 birey)

Bolgeden donen 445 SV'nin chr6 uzerinde KIF13A ile ortusen DEL'leri arasinda anlamli olanlar:

| variant_id | Koordinat | Uzunluk | AC | AN | AF | Hom | Genin %'si | Filtre | Sonuc |
|---|---|---|---|---|---|---|---|---|---|
| **DEL_chr6_33713a79** | 17.820.000-17.984.000 | **164.000 bp** | **1** | 126.090 | **8e-06** | 0 | **%71,3** | **PASS** | **lof** |
| DEL_chr6_5a9a59a9 | 17.771.790-17.774.520 | 2.730 bp | 1 | 126.092 | 8e-06 | 0 | %1,2 | PASS | lof |
| **DEL_chr6_81c144f1** | 17.774.225-17.775.574 | 1.349 bp | **43** | 126.092 | **3,41e-04** | **4** | %0,6 | PASS | lof |
| DEL_chr6_7e8bd330 | 17.785.167-17.785.750 | 583 bp | 175 | 120.420 | 1,45e-03 | 0 | %0,3 | LOWQUAL | lof |
| DEL_chr6_e943fc27 | 17.779.246-17.779.713 | 467 bp | 44 | 106.388 | 4,14e-04 | 0 | %0,2 | LOWQUAL | lof |

Iki nokta:
- DEL_chr6_33713a79 PASS filtrelidir ve genin %71,3'unu siler. KIF13A eksi ipliktedir (strand = -1), yani transkripsiyon baslangici yuksek koordinatta (~17.987.635) yer alir; 17.820.000-17.984.000 araligi promotor-proksimal ekzonlarin neredeyse tamamini kaldirir. Koordinatlarin yuvarlak olmasi bunun derinlik tabanli bir gCNV cagrisi oldugunu dusundurur, bu bir dogruluk cekincesidir, ama cagri filtreden gecmistir.
- DEL_chr6_81c144f1 "lof" olarak etiketlidir ama AF 3,41e-04 ve **4 homozigot** tasiyicisi vardir. DGV bunu bagimsiz olarak dogrular (asagida). Bu, gnomAD SV'nin KIF13A icin "lof" anotasyonunun her zaman gercek fonksiyon kaybi anlamina gelmedigini gosterir ve tabloyu her iki yone de dikkatle okumayi gerektirir.

### 1c. DGV

Cekilen kaynak: `https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dgvGold;chrom=chr6;start=17700000;end=18060000` ve ayni cagrinin `track=dgvMerged` hali.

**DGV Gold Standard:** bolgede 5 kayit, KIF13A ile ortusen 3 tanesi:

| Ad | Koordinat | Uzunluk | Tip | Genin %'si |
|---|---|---|---|---|
| gssvG31523 | 17.679.390-17.949.222 | 269.832 bp | **gain** | %83,0 |
| gssvL104913 | 17.949.314-17.951.700 | 2.386 bp | loss | %1,0 |
| gssvL104909 | 17.774.222-17.775.586 | 1.364 bp | loss | %0,6 |

**DGV Gold Standard'da tum geni kapsayan bir delesyon yoktur.** En buyuk kayip 2.386 bp'dir. Bu, iddia lehine olan ve durustce belirtilmesi gereken tek populasyon bulgusudur. Ancak DGV Gold, gnomAD v4'ten cok daha kucuk ve cok daha dusuk cozunurluklu dizi (array) verisine dayanir, yani orada yokluk gnomAD'daki varligi gecersiz kilmaz.

**DGV Merged:** 57 ortusen kayit. KIF13A'yi ilgilendiren en buyuk kayip nsv830597 (17.953.707-18.123.601, genin %15,0'i, Wong et al 2007, n=95, 1 kayip). DGV Merged'daki 1.352-1.364 bp'lik delesyon 1000 Genomes Faz 3'te 5/2.504 ve Faz 1'de 3/1.151 gozlenmistir (esv3608224, esv2658979), ki bu gnomAD'daki DEL_chr6_81c144f1 ile ayni polimorfizmdir.

### 1d. Toplam

gnomAD v4'te taranan toplam birey sayisi yaklasik **527.343** (63.046 genom + 464.297 ekzom). Bunlarin arasinda KIF13A'nin buyuk bolumunu veya tamamini kaldiran en az **3 bagimsiz delesyon** vardir (436229__DEL, 436227__DEL, DEL_chr6_33713a79). Kaba bir tahminle KIF13A'yi ciddi bicimde silen bir CNV sikligi ~1/175.000 mertebesindedir. Bu, "hicbir saglikli bireyde yok" ifadesini gecersiz kilar.

Cekince, iddianin lehine olacak sekilde belirtilmelidir: gnomAD v4 saf bir "saglikli" kohort degildir; yaygin hastalik kohortlarini icerir ama siddetli pediatrik hastaliklari disarida birakmayi hedefler. Yine de, tam penetransli bir NDD geni icin bu lokuste sifir yikici delesyon beklenir.

---

## 2. dbVar nstd101 ve bu cagri

Cekilen kaynaklar:
- https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd101/
- https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd37/
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=dbvar&term=nsv2774102`
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=dbvar&id=48453301`
- https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz

**Bulgu 1: nstd101 artik gecerli bir kayit degil.** dbVar sayfasinin basligi birebir soyle: "nstd101 (obsolete - ClinGen Kaminsky et al. 2011)". Sayfadaki alanlar: `Study Type: Case-Set`, `Date Obsoleted: 2021-11-03`, `Replacement Study: nstd102`, `3105 Variant Regions, 3833 Variant Calls`. Kuratoryal not: "PLEASE NOTE: Replaced nstd101 variant regions with variants in nstd102". E-utilities `db=dbvar` uzerinde `nstd101` sorgusu artik **0 sonuc** dondurur.

**Bulgu 2: nstd101 tasarimi geregi kontrolsuzdur.** `Study Type: Case-Set`. Kontrol kolu yoktur. Dolayisiyla "bu delesyon nstd101 kontrollerinde gorulmedi" seklinde bir argumanin hicbir zaman kurulamayacagi anlamina gelir. Kayit ile nstd101 arasindaki bag, ClinVar yorumundaki sablon bir cumleden ibarettir ve ozel bir vaka-kontrol analizini temsil etmez.

**Bulgu 3: gonderinin asil calismasi nstd101 degil, nstd37'dir.** ClinVar XML'inde `<SubmissionName>nstd37</SubmissionName>` yazar. nstd37 de "obsolete - ClinGen Laboratory-Submitted", `Study Type: Case-Set`, `Date Obsoleted: 2021-11-17`, `33378 Variant Calls`, yayin olarak Miller et al. 2010. Yani iddiada Kaminsky 2011'e (PMID 21844811) yapilan atif, bu gozlem icin degil, veri setinin genel kunyesi icindir.

**Bulgu 4: bu spesifik cagri.** `nsv2774102` artik `nsv3889946` / study `nstd102` altindadir, GRCh37 17.758.616-18.011.991, GRCh38.p12 17.758.385-18.011.760, `variant_call_count: 1`, listelenen genler KIF13A, RPS12P12, LOC105374953, klinik anlam "Uncertain significance". **Bu bolgede yalnizca tek bir cagri vardir.** Tekrarlayan bir 6p22.3 mikrodelesyonu degildir.

**Bulgu 5: nstd102'de (ClinVar CNV'lerinin dbVar aynasi) KIF13A ile ortusen tum cagrilar.** chr6 uzerinde 3.853 cagriden 21'i KIF13A ile ortusur: 14 kopya kazanci, 1 duplikasyon, 6 kopya kaybi. Kayiplarin dagilimi:

| Cagri | Koordinat (GRCh38) | Uzunluk | Klinik anlam | Kalitim (origin) |
|---|---|---|---|---|
| nssv15124161 | 5.354.169-17.949.848 | 12,60 Mb | Pathogenic | unknown |
| nssv18329282 | 12.005.397-22.849.418 | 10,84 Mb | Pathogenic | unknown |
| nssv20435443 | 12.020.008-22.846.454 | 10,83 Mb | Pathogenic | unknown |
| nssv15151127 | 13.693.620-24.225.287 | 10,53 Mb | Pathogenic | **de novo** |
| nssv15773065 | 13.909.894-21.999.975 | 8,09 Mb | Pathogenic | unknown |
| **nssv15150680** | **17.758.385-18.011.760** | **0,25 Mb** | **Uncertain significance** | **not provided** |

Yani dbVar'da 6p22.3'te benzer boyutta bagimsiz cagri **yoktur**; tek olan bizim kaydimizdir. Patojenik olan her kayip 8 Mb'in uzerindedir. Hicbir kontrol/saglikli ornek dbVar'in bu case-set'lerinde zaten bulunmaz, cunku bu veri setleri kontrol icermez.

---

## 3. Kalitim

Cekilen kaynak: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=442594&from_esearch=true`

XML'de bu gozleme dair **bulunan her alan** (9.278 bayt, tam kayit):

| Alan | Deger |
|---|---|
| VariationID / Accession | 442594 / VCV000442594, Version 2 |
| VariationName | GRCh37/hg19 6p22.3(chr6:17758616-18011991)x1 |
| VariantType | copy number loss, `AbsoluteCopyNumber = 1` |
| variantLength | 253.376 |
| NumberOfSubmissions / NumberOfSubmitters | **1 / 1** |
| DateCreated / DateLastUpdated | 2017-10-26 / 2024-05-08 |
| SCV | SCV000585270, Version 1, SubmissionDate **2017-09-26** |
| Submitter | ISCA site 1, OrgID 505237, laboratory; ikincil: ISCA Consortium, OrgID 500029 |
| Classification | Uncertain significance, `no assertion criteria provided`, **DateLastEvaluated = 2014-05-05** |
| Assertion | variation to disease |
| **`<Origin>`** | **`unknown`** |
| `<AffectedStatus>` | `yes` |
| `<Species>` | human |
| ObservedData | `VariantAlleles = 1` (tek allel; test edilen birey sayisi **verilmemis**) |
| Method 1 | Affymetrix CytoScan HD, SNP array, Purpose = Discovery, MethodType = clinical testing |
| Method 2 | FISH, Purpose = Validation, MethodResult = **Pass** |
| Klinik ozellik | "Delayed speech and language development", HP:0000750, `ClinicalFeaturesAffectedStatus = present` |
| TraitSet (kosul) | "See cases" (PhenotypeInstruction) |
| Atif | PubMed 20466091 (Miller et al. 2010), iki kez tekrarlanmis |
| dbVar capraz referanslari | nssv13656023 (call), nsv2774102 (region) |
| SubmissionName | nstd37 |
| Comment | "Copy number variation identified through the course of routine clinical cytogenomic testing in postnatal populations, with clinical assertions as classified by the original submitter. For data from the original published study, [Kaminsky, et al. 2011], please see [nstd101]." |

**XML'de BULUNMAYAN alanlar, acikca:** de novo / maternal / paternal ayrimi (Origin = `unknown`), ebeveyn testi yapilip yapilmadigi, aile oykusu veya segregasyon, **yas**, **cinsiyet**, etnik koken, test edilen birey sayisi, **es-mevcut (co-occurring) varyant** bilgisi, ek fenotipik ozellik, gelisim basamaklari, zeka olcumu, dismorfizm, ve vakaya ozgu herhangi bir **gonderen yorumu**. `<Comment>` alanindaki metin, tum nstd37 gonderimlerinde birebir ayni olan bir sablondur; bu vaka hakkinda hicbir sey soylemez.

**Bunun maliyeti.** Kalitimin bilinmemesi, kaydin kanit degerinin buyuk kismini yok eder:
- ClinGen CNV yorumlama cerceveisinde de novo teyidi (bolum 4) en agirlikli puan kaynaklarindan biridir; `unknown` origin sifir puan verir.
- Delesyon, etkilenmemis bir ebeveynden kalitilmis olabilir, ki gnomAD'daki 448 kb'lik kontrol delesyonu bunun gercekci bir senaryo oldugunu gosterir.
- FISH dogrulamasi yalnizca delesyonun **gercek** oldugunu gosterir, **nedensel** oldugunu degil. Iddiada "microarray + FISH" bir guc unsuru gibi sunulmus; teknik olarak sadece teknik artefakti disliyor.
- Siniflamanin son degerlendirme tarihi 2014-05-05'tir, yani gnomAD SV/CNV var olmadan once. Kayit 12 yildir yeniden degerlendirilmemistir.

---

## 4. ClinGen doz duyarliligi

Cekilen kaynaklar:
- https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv (dosya basligi: "#ClinGen Gene Curation Results", "#16 Sep,2026", 1.540 satir)
- https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv (519 satir)
- https://search.clinicalgenome.org/kb/gene-dosage/KIF13A (sayfa doner, `search.clinicalgenome.org/api/dosage/KIF13A` HTTP 404)

**KIF13A ClinGen doz duyarliligi listesinde YOKTUR.** `^KIF13A\t` ile tam eslesme sifir sonuc verir; buyuk/kucuk harf duyarsiz arama da sifir verir. Yani ne bir haploinsufficiency skoru, ne bir triplosensitivity skoru, ne bir tarih, ne de bir gerekce metni mevcuttur. Kurate edilecek gen kuyrugunda bile degildir.

Kalibrasyon icin, ayni dosyadan canli cekilen karsilastirmalar:

| Gen | HI skoru | HI aciklamasi | Son degerlendirme |
|---|---|---|---|
| **KIF13A** | **(kayit yok)** | **(kayit yok)** | **(kayit yok)** |
| KIF1A | 0 | "No evidence available" | 2013-08-08 |
| KIF11 | 3 | (yeterli kanit) | 2013-08-01 |
| KIF1C | 30 | (otozomal resesif) | 2016-08-22 |
| KIF7 | 30 | (otozomal resesif) | 2016-08-22 |
| NRXN1 | 3 | "Sufficient evidence for dosage pathogenicity" | 2021-06-08 |
| SHANK3 | 3 | "Sufficient evidence for dosage pathogenicity" | 2018-08-14 |

Dikkat cekici olan sudur: **KIF1A bile** (iyi tanimlanmis, klinik olarak karakterize edilmis paralog) haploinsufficiency skoru **0 = "No evidence available"** almistir, cunku KIF1A hastaligi agirlikli olarak missense/dominant-negatif ve resesif mekanizmalarla ilerler, tam gen kaybiyla degil. Bu, paralog kalibrasyonlu bir proje icin dogrudan uyaricidir: kinesin-3 ailesinde LoF-haploinsufficiency modeli varsayilan mekanizma degildir.

**Bolge kurasyonu:** ClinGen bolge listesinde 6p22.3 icin kurate edilmis hicbir bolge yoktur. `6p22` aramasi yalnizca iki 6p22.1 kaydi dondurur, ikisi de populasyon bolgesidir ve skoru **40 = "Dosage sensitivity unlikely"**: ISCA-46608 (chr6:29.869.423-29.881.496) ve ISCA-46607 (chr6:29.924.783-29.938.057). Bunlar KIF13A ile ilgisizdir ve lokusun ~12 Mb sentromerik tarafindadir.

---

## 5. Baska tek genli KIF13A delesyonlari

Cekilen kaynak: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=KIF13A[gene] AND "copy number loss"[Type of variation]` ve esummary.

**Iddia kismen dogru, sayilar duzeltilmeli.**

ClinVar'da KIF13A'ya anote edilmis **kopya sayisi kaybi** kayitlarinin tamami **6 tanedir** (iddiadaki "20" sayisi, kayiplar + kazanclar toplami olan 19 ile karistirilmis olmali; `KIF13A[gene] AND (copy number loss OR copy number gain)` = **19**, `KIF13A[gene]` toplami = **313**).

| VCV | Boyut | ClinVar'in listeledigi gen sayisi | Siniflama | Inceleme durumu | Fenotip |
|---|---|---|---|---|---|
| **VCV000442594** | **253.376 bp** | **1 (KIF13A)** | **Uncertain significance** | no assertion criteria provided | See cases |
| VCV000688193 | 8.090.080 bp | 23 | Pathogenic | no assertion criteria provided | not provided |
| VCV000442457 | 10.531.664 bp | 29 | Pathogenic | no assertion criteria provided | See cases |
| VCV004755369 | 10.826.447 bp | **382** | Pathogenic | criteria provided, single submitter | See cases |
| VCV001809386 | 10.844.018 bp | 34 | Pathogenic | criteria provided, single submitter | not provided |
| VCV000253430 | 12.595.678 bp | 53 | Pathogenic | criteria provided, single submitter | See cases |

**Duzeltme:** "her patojenik olan 23-237 gen kapsiyor" ifadesi yerine dogrusu **23 ila 382 gen** olmalidir (VCV004755369 GRCh38 uzerinde kaydedilmistir ve ClinVar daha eksiksiz bir gen listesi cikarir; ayni fiziksel delesyonun GRCh37 esdegeri VCV001809386 yalnizca 34 gen listeler). Asil nokta degismiyor: **patojenik kayiplarin en kucugu 8,09 Mb'dir, yani bizim kaydimizdan 32 kat buyuktur** ve hicbiri KIF13A'yi tek aday olarak birakmaz. 6p24-p22 araliginda JARID2, ATXN1, SOX4, TFAP2A gibi cok daha iyi kurulmus NDD adaylari bulunur.

**Cok daha sert bir ek bulgu.** ClinVar'daki 313 KIF13A kaydinin tamaminin siniflama dagilimi:

| Siniflama | Sayi |
|---|---|
| Uncertain significance | **221** |
| (siniflanmamis) | 62 |
| Pathogenic | 15 |
| Likely benign | 8 |
| Benign | 7 |
| **Likely pathogenic** | **0** |

Ve **15 patojenik kaydin 15'i de** cok megabazlik kromozom 6p CNV'sidir (chr6:106.431-18.360.595 x3'ten chr6:156.975-170.919.482 x3'e kadar; listenin tamami dogrulandi). **ClinVar'da bugune kadar patojenik veya muhtemelen patojenik olarak siniflanmis tek bir intragenik KIF13A varyanti yoktur.** Nonsense varyant sorgusu sifir sonuc verir. NCBI Gene 63971 kaydinda KIF13A icin listelenmis hicbir fenotip yoktur (MIM 605433 yalnizca gen kaydidir, fenotip MIM'i degildir).

**Adil karsi-agirlik:** PubMed'de KIF13A + delesyon/haploinsufficiency/norogelisimsel sorgusu yalnizca 7 kayit dondurur. Bunlarin icinde tek gercek insan vaka raporu **PMID 33526817** (Lam CW ve ark., J Hum Genet 2021;66(8):825-829): hipomelanozis of Ito ve gelisim geriligi olan 3 yasinda bir kizda **de novo, postzigotik mozaik** KIF13A cerceve kaymasi (NM_022113.6:c.2357dupA, NMD ile sonuclanan). Bu kayit iddiaya yardim eder ama uc cekince tasir: (1) varyant **mozaiktir**, dolayisiyla constitutional heterozigot haploinsufficiency'nin dogrudan kaniti degildir; (2) fenotip agirlikla **noroektodermaldir** (cilt hipopigmentasyonu, melanozom biyogenezi ile tutarli), saf bir NDD degildir; (3) yazarlarin kendi gerekcesi dairesel bir akil yurutmedir, cunku destek olarak **6p22.3-p23 mikrodelesyon sendromunu** gosterirler, ki bu tam olarak burada test ettigimiz cok genli cikarim problemidir.

---

## 6. Duzenleyici (regulatory) karisitirici

Cekilen kaynaklar:
- `https://grch37.rest.ensembl.org/overlap/region/human/6:17300000-18500000?feature=gene`
- `https://rest.ensembl.org/overlap/region/human/6:17300000-18500000?feature=gene`
- `https://rest.ensembl.org/overlap/region/human/6:17758385-18011760?feature=regulatory`
- `https://rest.ensembl.org/overlap/region/human/6:17706942-17758385?feature=regulatory`
- `POST https://gnomad.broadinstitute.org/api` (gene constraint)

**Komsu protein kodlayan genler ve mesafeler (GRCh37; delesyon 17.758.616-18.011.991):**

| Gen | Koordinat (GRCh37) | Ip | Delesyona konum | Mesafe |
|---|---|---|---|---|
| NUP153 | 17.615.269-17.706.656 | -1 | telomerik (once) | **51.960 bp** |
| KIF13A | 17.759.414-17.987.854 | -1 | **delesyonun icinde** | 0 |
| NHLRC1 | 18.120.718-18.122.851 | -1 | sentromerik (sonra) | **108.727 bp** |
| TPMT | 18.128.542-18.155.305 | -1 | sentromerik | 116.551 bp |
| KDM1B | 18.155.560-18.224.084 | +1 | sentromerik | 143.569 bp |

Ensembl GRCh37 ve GRCh38 bagimsiz olarak dogruluyor: NUP153 ile NHLRC1 arasindaki tum aralikta **tek protein kodlayan gen KIF13A'dir**. Yani projenin daha once yaptigi tespit teyit edilmistir.

**Delesyonun icindeki duzenleyici elemanlar (Ensembl Regulatory Build, GRCh38 17.758.385-18.011.760):** toplam **48** eleman:

| Tip | Sayi |
|---|---|
| enhancer | **37** |
| CTCF_binding_site | **9** |
| open_chromatin_region | 1 | 
| promoter | 1 |

Ornekler: ENSR6_993WX (17.815.231-17.815.660, enhancer), ENSR6_993XH (17.817.285-17.817.616, enhancer), ENSR6_5673TL (17.804.357-17.804.378, CTCF), ENSR6_5678HQ (17.831.388-17.831.409, CTCF), ENSR6_8J4SP (17.768.442-17.768.647, acik kromatin).

NUP153'un 3' ucu ile delesyonun baslangici arasindaki 51,4 kb'lik bosluk ayrica **9 duzenleyici eleman** icerir.

**Bunun anlami.** Delesyon 9 CTCF baglanma bolgesini kaldirir. CTCF bolgeleri TAD sinirlarinin yapi taslaridir, dolayisiyla bu delesyonun yerel kromatin alan yapisini degistirmesi mekanik olarak mumkundur. Kirilma noktasinin NUP153 promotoruna 52 kb uzakta olmasi onemlidir, cunku **NUP153 KIF13A'dan daha kisitlidir**:

| Gen | pLI | LOEUF (oe_lof_upper) | gozlenen/beklenen pLoF |
|---|---|---|---|
| **NUP153** | **1,000** | **0,277** | **27 / 134,6** |
| KIF13A | 0,9998 | 0,449 | 79 / 212,2 |
| KIF1A | 1,000 | 0,259 | 46 / 227,3 |
| KIF13B | 1,4e-20 | 0,663 | 127 / 222,0 |
| NHLRC1 | 0,003 | 1,166 | 10 / 14,6 |
| TPMT | 0,003 | 0,781 | 18 / 34,2 |

Yani "delesyonun icindeki tek protein kodlayan gen KIF13A'dir" cumlesi dogru olsa da, "dolayisiyla fenotip KIF13A'dandir" cikarimi gecerli degildir: aralik, 52 kb otede duran ve nukleoporin isleviyle norogelisimle iliskili, KIF13A'dan daha LoF-kisitli bir genin potansiyel uzak duzenleyici bolgesini de kaldirmaktadir.

Spekulasyon siniri: Bu enhancer'lardan hangisinin hangi gene bagli oldugunu gosteren dokuya ozgu bir 3C/Hi-C veya eQTL kaniti **cekilmedi ve iddia edilmiyor**. Soylenebilecek olan sudur: KIF13A-disi bir mekanizma, mevcut verilerle **dislanmamistir**.

Ayrica: NHLRC1 (EPM2B, Lafora hastaligi) delesyonun **disindadir** (108,7 kb uzakta), dolayisiyla resesif bir alelin aciga cikmasi (unmasking) senaryosu bu kayit icin gecerli degildir. KIF13A'nin kendisi icin boyle bir senaryo (delesyon + trans'ta bir noktasal varyant) **test edilmemistir**, cunku gonderide es-mevcut varyant bilgisi yoktur (bkz. bolum 3).

---

## 7. Fenotip olarak konusma ve dil geriliginin ozgullugu

Cekilen kaynaklar:
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&...` (fenotip/siniflama sayimlari)
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=21844811,20466091` (Kaminsky 2011, Miller 2010)
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=40250875` (BMJ Open 2025 DLD scoping review)

**"Delayed speech and language development" (HP:0000750) son derece ozgul olmayan bir gozlemdir.**

1. **Populasyon taban orani.** Canli cekilen PMID **40250875** (BMJ Open 2025;15(4):e089667, Arap ulkelerinde gelisimsel dil bozuklugu epidemiyolojisi, 17 calisma): bildirilen DLD prevalansi **%3,2** (Misir, 3-6 yas) ile **%25,6** (Suudi Arabistan, 1-5 yas) arasinda degismektedir. Yani bu fenotip, secilmemis cocuk populasyonlarinda bile yaygin bir bulgudur ve tek basina hicbir lokusa isaret etmez.

2. **Yonlendirme populasyonundaki taban orani.** Miller ve ark. 2010 (PMID **20466091**, ClinVar kaydinin kendi atfettigi makale) CMA'yi "aciklanamayan gelisim geriligi/entelektuel yetersizlik, otizm spektrum bozuklugu veya coklu konjenital anomali" olan bireyler icin birinci basamak test olarak tanimlar; 21.698 hastalik literatur derlemesinde tanisal verim **%15-20**'dir. Yani **yonlendirilen hastalarin %80-85'inde CMA nedensel bir bulgu vermez**, ama cogunda gelisimsel gecikme zaten vardir. Konusma/dil gecikmesi, CMA'ya yonlendirilme **sebebidir**, delesyonun **sonucu** degil. Bu, klasik bir ascertainment (secilim) yanliligidir: etkilenmis bir cocukta bulunan herhangi bir nadir CNV, tanimi geregi "etkilenmis bir cocukta bulunmus" olacaktir.

3. **Kaminsky 2011'in kendisi bu lokusu incelemedi.** Canli cekilen ozet (PMID **21844811**): 15.749 ISCA vakasi ve 10.118 yayimlanmis kontrol, ancak analiz "**14 tekrarlayan CNV bolgesine**" odaklanmistir ve 14 delesyon ile 7 duplikasyon anlamli bulunmustur. 6p22.3/KIF13A bu 14 bolgeden **biri degildir**. Dolayisiyla ClinVar kaydinin Kaminsky 2011'e atfi, bu varyant icin hicbir vaka-kontrol istatistigi tasimaz.

4. **ClinVar'in kendi CNV havuzunda sinyalin seyrekligi.** ClinVar'daki tum kopya sayisi kaybi + kazanci kayitlari **49.796**'dir; bunlarin **33.625'i (%67,5) "uncertain significance"** olarak siniflanmistir. Bizim kaydimiz bu 33.625'in icindedir ve icinde ayirt edici hicbir ozelligi yoktur. (Not: ClinVar'in fenotip indeksi taban orani hesabi icin kullanilamaz, cunku CNV kayitlarinin 23.821'i "See cases" altinda dosyalanmistir ve yalnizca 44 CNV kaydi hastalik alaninda HP:0000750 tasir; bu, gercek bir taban orani degil, indeksleme artefaktidir ve boyle raporlanmalidir.)

5. **Kaydin fenotip derinligi sifira yakindir.** Gonderide **tek** bir HPO terimi vardir. Entelektuel yetersizlik, otizm, nobet, mikrosefali, dismorfizm, MRG bulgusu, buyume parametresi yoktur. Hipopigmentasyon (PMID 33526817'deki KIF13A vakasinin ayirt edici bulgusu) da yoktur. Tek terimli bir fenotip, gen-hastalik iliskisi kurmak icin ClinGen Gene-Disease Validity cerceveisinde en dusuk agirlikli kanit tipidir.

---

## Bu kayit ne kadar ileri goturulebilir

### Kanitin gercekten izin verdigi en guclu cumle

> "ClinVar'da (VCV000442594 / SCV000585270) KIF13A'yi tek protein kodlayan gen olarak iceren tek bir 253 kb'lik 6p22.3 delesyonu kayitlidir; delesyon konusma ve dil gecikmesi olan bir cocukta mikrodizi ile bulunup FISH ile dogrulanmistir. Kalitim durumu bildirilmemistir, kayit tek bir gonderene aittir, resmi siniflamasi 'uncertain significance / no assertion criteria provided'dir ve 2014'ten beri yeniden degerlendirilmemistir. gnomAD v4'te KIF13A'nin neredeyse tamamini kaldiran en az uc bagimsiz delesyon bulundugu icin (bunlardan 436229__DEL bu vakanin araliginin %99,3'unu icerir), bu kayit KIF13A haploinsufficiency'sinin kaniti degil, ancak KIF13A'yi 6p22.3 delesyonlarinda **aday** olarak tutan minimal-aralik gozlemidir."

Yani izin verilen ifade bir **aday belirleme** ifadesidir, bir **nedensellik** ifadesi degil. KIF13A'nin aday olarak kalmasini destekleyen bagimsiz unsurlar, bu kayit degil, sunlardir: gnomAD LoF kisitliligi (pLI 0,9998, LOEUF 0,449, 79 gozlenen / 212,2 beklenen pLoF) ve PMID 33526817'deki de novo mozaik cerceve kaymasi vakasi.

### Asiri iddia olacak cumle

> ~~"ClinVar VCV000442594, KIF13A'nin bir kopyasinin kaybinin bir norogelisimsel bozukluga yol actigina dair en guclu insan kanitidir."~~

Bu cumle su altı nedenle savunulamaz:
1. gnomAD v4 ekzomlarinda (464.297 birey) **436229__DEL**, bu araligin %99,3'unu ve KIF13A'nin %99,2'sini iceren 448 kb'lik bir delesyondur ve siddetli pediatrik hastalik icin secilmemis bir bireyde bulunur. Buna gnomAD genomlarindaki PASS filtreli 164 kb'lik DEL_chr6_33713a79 (genin %71,3'u) eklenir.
2. Kalitim bilinmiyor (`Origin: unknown`). De novo teyidi olmadan, etkilenmemis bir ebeveynden kalitilmis olma ihtimali dislanamaz.
3. ClinGen KIF13A'yi hic kurate etmemistir; doz duyarliligi skoru **yoktur** (kurate edilmis 1.540 genin icinde degildir). Paralog KIF1A bile HI skoru 0 = "No evidence available" almistir.
4. Baglandigi calisma **nstd101 obsoletedir** (2021-11-03), tasarimi geregi **Case-Set**'tir (kontrol kolu yoktur) ve atfedilen Kaminsky 2011 analizi bu lokusu degil, yalnizca 14 tekrarlayan CNV bolgesini incelemistir.
5. Delesyon KIF13A'nin yaninda **37 enhancer ve 9 CTCF bolgesi** kaldirir ve kirilma noktasi KIF13A'dan daha kisitli olan NUP153'un (pLI 1,000, LOEUF 0,277) promotoruna 52 kb uzaktadir; KIF13A-disi bir mekanizma dislanmamistir.
6. Fenotip tek bir HPO terimidir ve CMA'ya yonlendirilme sebebiyle ozdestir; ClinVar'daki tum CNV kayitlarinin %67,5'i zaten VUS'tur ve ClinVar'da bugune kadar **sifir** intragenik KIF13A varyanti patojenik veya muhtemelen patojenik olarak siniflanmistir.

### Projenin bunu duzeltmek icin yapabilecegi seyler

- ClinVar'dan bu SCV'nin gonderenine (ISCA site 1 / Geisinger, Erin Riggs) kalitim durumu icin dogrudan sorulabilir; `Origin: unknown` alani laboratuvarda cozulmus olabilir.
- DECIPHER ve DDD, bu raporda **sorgulanmadi**; bagimsiz bir tek genli KIF13A delesyonu orada olabilir ve tablosu degistirebilir.
- Asil yuk, CNV'de degil, missense/LoF ekseni uzerindedir: 221 VUS ve 0 patojenik intragenik varyant, projenin paralog kalibrasyonlu onceliklendirme yonteminin tam da doldurmak icin tasarlandigi bosluktur.

---

## Cekilen URL'lerin tam listesi

| Amac | URL / cagri |
|---|---|
| gnomAD SV v4 ve CNV v4 | `POST https://gnomad.broadinstitute.org/api` (GraphQL; `gnomad_sv_r4`, `gnomad_cnv_r4`, bolge 6:17700000-18060000 GRCh38); tarayici: https://gnomad.broadinstitute.org/region/6-17700000-18060000?dataset=gnomad_sv_r4 |
| gnomAD gen kisitliligi | `POST https://gnomad.broadinstitute.org/api` (KIF13A, NUP153, NHLRC1, TPMT, KIF13B, KIF1A) |
| ClinVar tam XML | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&is_variationid&id=442594&from_esearch=true |
| ClinVar arama/ozet | `esearch.fcgi` / `esummary.fcgi`, `db=clinvar`, `KIF13A[gene]` ve tur/siniflama filtreleri |
| dbVar nstd101 | https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd101/ |
| dbVar nstd37 | https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd37/ |
| dbVar varyant ozeti | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=dbvar&id=48453301 |
| dbVar nstd102 tam tablo | https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz |
| ClinGen gen doz listesi | https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv |
| ClinGen bolge doz listesi | https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv |
| ClinGen KIF13A sayfasi | https://search.clinicalgenome.org/kb/gene-dosage/KIF13A (API yolu `.../api/dosage/KIF13A` HTTP 404) |
| DGV Gold Standard | https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dgvGold;chrom=chr6;start=17700000;end=18060000 |
| DGV Merged | https://api.genome.ucsc.edu/getData/track?genome=hg38;track=dgvMerged;chrom=chr6;start=17700000;end=18060000 |
| Ensembl genler GRCh37 | https://grch37.rest.ensembl.org/overlap/region/human/6:17300000-18500000?feature=gene |
| Ensembl genler GRCh38 | https://rest.ensembl.org/overlap/region/human/6:17300000-18500000?feature=gene |
| Ensembl duzenleyici (delesyon ici) | https://rest.ensembl.org/overlap/region/human/6:17758385-18011760?feature=regulatory |
| Ensembl duzenleyici (NUP153 boslugu) | https://rest.ensembl.org/overlap/region/human/6:17706942-17758385?feature=regulatory |
| NCBI Gene KIF13A | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=63971 |
| PubMed ozetleri | `efetch.fcgi?db=pubmed&id=21844811,20466091,33526817,40250875,21792059&rettype=abstract` |

Canli cekilip okunan PMID'ler: **21844811** (Kaminsky 2011), **20466091** (Miller 2010), **33526817** (Lam 2021, KIF13A mozaik vakasi), **40250875** (BMJ Open 2025, DLD prevalansi), **21792059** (Sahoo 2011). DGV Merged kayitlarinda referans olarak gecen PMID'ler (25217958 Coe 2014, 21841781 Cooper 2011, 21293372 1000 Genomes, 32461652 gnomAD SV, 19592680 Shaikh 2009, 17160897 Wong 2007, 23128226 1000G Faz 1, 25118596, 25503493, 24416366, 24896259, 19812545) UCSC track verisinden gelmektedir; bu makalelerin metinleri okunmamistir, yalnizca DGV kaydindaki atif alani raporlanmistir.

---
---

# EK BOLUM: DECIPHER, DDD, nstd102 ve gnomAD alt kohortlari

Tarih: 2026-09-16, ilk raporla ayni gun, koordinatorun "kendi isaret ettigin boslugu kapat" talimati uzerine eklendi.

Not: proje venv yorumlayicisi (`mcp\atlas_server\.venv\Scripts\python.exe`) bu asamada Windows Smart App Control tarafindan engellendi ("Permission denied" / "Unable to create process"), bu yuzden bu bolumun betikleri sistem Python 3.13.14 ile yalnizca standart kutuphane kullanilarak calistirildi. Veri kaynaklari ve sayilar etkilenmedi.

## Karar degisti mi?

**Hayir, karar seviyesi ayni kaliyor: "suggestive" (telkin edici), skalanin en alt basamagi. Ancak gerekce kumesi degisti ve iddianin aleyhine agirlasti.**

Koordinatorun sordugu spesifik soru ("DECIPHER veya DDD'de bagimsiz bir tek genli KIF13A delesyonu var mi, varsa tablo degisir mi") **iki tarafli** bir cevap aldi:

- **Lehte, yeni:** DECIPHER'da KIF13A'nin tek protein kodlayan gen oldugu **ikinci** bir hasta delesyonu vardir (hasta 451718, 146.274 bp, tamamen gen ici, VUS, fenotip "Autistic behavior"). Yani ClinVar kaydi artik turunun tek ornegi degildir.
- **Aleyhte, daha agir:** Ayni DECIPHER veri setinde KIF13A'nin yine tek gen oldugu **ucuncu** bir delesyon vardir (hasta 504413, 182.672 bp, genin **%73,0'unu** siler) ve bu delesyon **anneden kalitilmistir** ve yatiran merkez tarafindan **"likely benign"** olarak siniflanmistir. Buna ek olarak DDD arastirma verisindeki **iki** KIF13A varyantinin **ikisi de "Likely LOF"** (bir stop_gained, bir frameshift), **ikisi de anneden kalitilmis**, **ikisi de "Uncertain"** olarak birakilmistir. Yani gelisimsel bozukluk kohortlarinda bulunan dort bagimsiz KIF13A fonksiyon kaybi olayindan **ucunun kalitim durumu bilinmektedir ve ucu de maternaldir**; hicbiri patojenik olarak siniflanmamistir.

Ayrica KIF13A **DDG2P listesinde yoktur**, DECIPHER'in kendi gen sayfasinda Gene2Phenotype, ClinGen gene/disease, ClinGen Dosage Sensitivity, GenCC ve OMIM Morbid alanlarinin **hepsi bostur**, ve DECIPHER'in tasidigi haploinsufficiency tahmini (%HI) KIF13A'yi **%42,3** persentiline koyar, yani "haploinsufficiency'ye yatkin" esiginin (%10) cok disinda.

Tek gercek karsi-agirlik: DECIPHER'in tasidigi **pHaplo = 0,99** ve pLI = 0,99. Bu metrik iddiayi destekler ve boyle raporlanmalidir; ancak ayni sayfadaki sHet = 0,049 ve %HI = %42,3 ayni yone isaret etmez. Bu ic tutarsizlik, KIF13A icin haploinsufficiency tahminlerinin **uzlasmis olmadigini** gosterir.

**Bolum "Karar"a eklenen cumle:** yukaridaki maternal kalitim bulgulari nedeniyle, ilk karardaki "kalitim bilinmiyor" maliyeti artik yalnizca bir bosluk degil, aktif olarak ters yonde bir sinyaldir.

---

## EK-1. DECIPHER

Cekilen kaynaklar:
- https://www.deciphergenomics.org/gene/KIF13A/overview/clinical-info (DECIPHER v11.40, tarayicida render edildi)
- https://www.deciphergenomics.org/data/gene/ENSG00000137177/patient-overlap (acik erisim JSON, 88.611 bayt)
- https://www.deciphergenomics.org/data/gene/ENSG00000137177/counts
- https://www.deciphergenomics.org/data/gene/ENSG00000137177/dddrv-overlap
- https://www.deciphergenomics.org/files/downloads/HI_Predictions_Version3.bed.gz
- https://www.deciphergenomics.org/files/downloads/population_cnv_grch38.txt.gz

### EK-1a. Gen duzeyinde doz ve hastalik metrikleri

DECIPHER gen sayfasinin "Gene/disease association" blogu, birebir:

| Alan | Deger |
|---|---|
| Gene2Phenotype | **-** |
| OMIM | 605433 (gen kaydi) |
| **Morbid** | **-** |
| GeneReviews | **-** |
| **ClinGen gene/disease** | **-** |
| **ClinGen Dosage Sensitivity** | **-** |
| **GenCC** | **-** |

Bes ayri gen-hastalik kaynagi da bos. Bu, ilk rapordaki ClinGen FTP bulgusunu bagimsiz olarak dogrular ve genisletir.

"Gene predictive scores" blogu:

| Metrik | Deger | Yon |
|---|---|---|
| pLI | 0,99 | iddia lehine |
| LOEUF | 0,47 | notr, kisitli ama ekstrem degil |
| sHet (heterozigot LoF intoleransi) | **0,049** | zayif, iddia aleyhine |
| **pHaplo** | **0,99** | **iddia lehine** |
| pTriplo | 0,60 | ilgisiz |
| Missense Z | 2,27 | notr |

"Protein predictive scores": pDN (dominant-negatif olasiligi) **0,560**, pGOF **0,598**, pLOF **0,616**. Ucu de birbirine yakin, yani mekanizma tahmini KIF13A icin **ayirt edici degildir**; saf LoF mekanizmasi tercih edilmemektedir.

### EK-1b. DECIPHER %HI (Haploinsufficiency Predictions v3)

`HI_Predictions_Version3.bed.gz` dosyasindan birebir satirlar:

| Gen | Ham skor | **%HI** | BED rengi |
|---|---|---|---|
| NUP153 | 0,346309412 | **%24,4** | 175,80,0 |
| **KIF13A** | 0,154932935 | **%42,3** | 125,130,0 |
| KIF1A | 0,102859335 | %50,62 | 100,155,0 |
| NHLRC1 | 0,071688092 | %57,06 | 100,155,0 |
| KIF13B | 0,033997776 | %67,78 | 75,180,0 |

DECIPHER'in kullandigi esikte **%HI < %10** olan genler "haploinsufficiency gosterme olasiligi yuksek" kabul edilir. KIF13A **%42,3** ile bu esigin cok uzagindadir. Ustelik ilk raporda duzenleyici karisitirici olarak isaret edilen **NUP153, KIF13A'dan daha haploinsufficiency'ye yatkin** siralanmaktadir (%24,4 karsisinda %42,3). Bu, ilk rapordaki 52 kb komsuluk argumanini zayiflatmak yerine guclendirir.

### EK-1c. DECIPHER hasta CNV'leri: genel dagilim

`counts` uc noktasi: `{"patient-overlap":55, "ddd-research-variant-overlap":2}`. Sayfa bunu 0 dizi varyanti + **47 kopya sayisi varyanti** + 8 "diger varyant" olarak ayirir. Sayfa ayrica acikca sunu yazar: "**DECIPHER holds no open-access sequence variants in this gene**" ve "Currently there is not enough open-access patient data to present other quantitative information about this gene."

47 CNV'nin dagilimi:

| Ozellik | Dagilim |
|---|---|
| Sinif | 29 delesyon, 16 duplikasyon, 2 triplikasyon |
| Patojenite | 21 kayit yok, 9 likely_pathogenic, 9 pathogenic, 7 uncertain, **1 likely_benign** |
| Kalitim | 27 unknown, 16 de_novo, 1 de_novo_confirmed, 1 inherited_from_affected_parent, **1 maternal**, 1 imbalance_from_balanced_rearrangement |

Yalnizca 29 delesyon icin: patojenite 15 kayit yok, 6 likely_pathogenic, 5 pathogenic, 2 uncertain, 1 likely_benign; kalitim 16 unknown, 11 de_novo, 1 de_novo_confirmed, **1 maternal**.

### EK-1d. KIF13A'yi tek protein kodlayan gen olarak iceren DECIPHER delesyonlari

Referans: KIF13A GRCh38 17.758.198-17.988.274. Komsu protein kodlayan genler NUP153 (biter 17.706.942) ve NHLRC1 (baslar 18.120.440). Yani **17.706.943-18.120.439** araliginda kalan her delesyon KIF13A-tek-gen delesyonudur.

DECIPHER'in 29 delesyonundan **ikisi** bu kosulu saglar:

| Hasta | Koordinat (GRCh38) | Uzunluk | KIF13A'nin %'si | Patojenite | Katki | **Kalitim** | Cinsiyet | Fenotip |
|---|---|---|---|---|---|---|---|---|
| **451718** | 17.825.889-17.972.163 | 146.274 bp | %63,6 (**tamamen gen ici**) | **uncertain** | - | **unknown** | 46_XY | **Autistic behavior** |
| **504413** | 17.743.555-17.926.227 | 182.672 bp | **%73,0** | **likely_benign** | none | **maternal** | 46_XY | **(fenotip kaydi yok)** |

Bu iki kayit tam olarak birbirini dengeler ve net etki iddianin aleyhinedir, cunku ikincisinde **kalitim gercekten olculmustur** ve sonuc maternaldir. KIF13A eksi ipliktedir ve transkripsiyon baslangici yaklasik 17.987.635'tedir; 504413'un delesyonu genin dusuk koordinatli, yani **3' yarisinin tamamini ve kodlayan bolgenin buyuk kismini** siler, dolayisiyla "kucuk, zararsiz bir ucu kirpma" olarak aciklanamaz.

Ucuncu en kucuk delesyon KIF13A-tek-gen degildir ve iddia lehine gibi gorunse de cok genlidir:

| Hasta | Koordinat | Uzunluk | Patojenite | Katki | Kalitim | Fenotip |
|---|---|---|---|---|---|---|
| 358696 | 17.675.564-18.225.293 | 549.729 bp | likely_pathogenic | **full** | **de_novo** | Global developmental delay, IUGR, Short stature, Strabismus |

Bu aralik KIF13A'nin yani sira **NHLRC1, TPMT, KDM1B**'yi tamamen, **NUP153** ve **DEK**'i kismen kapsar. "contribution: full" etiketi yatiran merkezin yorumudur, tek gen atfi degildir.

Geri kalan 26 delesyonun tamami **1,21 Mb ile 23,07 Mb** arasindadir, yani ilk rapordaki ClinVar tablosuyla ayni yapiyi gosterir: patojenik sinyal yalnizca cok genli, megabaz olcekli delesyonlarda vardir. Ornek olarak, "Delayed speech and language development" fenotipi tasiyan iki DECIPHER hastasinin delesyonlari **12,90 Mb** (hasta 396348) ve **23,07 Mb**'dir (hasta 398397).

Duplikasyon tarafinda not edilmesi gereken bir kayit: hasta 270332, chr6:17.556.955-17.825.945 (268.990 bp) duplikasyonu, kalitim **"inherited_from_affected_parent"**, patojenite kaydi yok.

### EK-1e. DECIPHER populasyon CNV izi

`population_cnv_grch38.txt.gz` (57.639 kayit) icinde chr6:17.700.000-18.100.000 araliginda 14 kayit vardir. KIF13A ile ortusen **en buyuk delesyon yalnizca 11.902 bp'dir** (chr6:17.829.938-17.841.840, 1/5.919, Affy6 calismasi). DECIPHER'in populasyon setinde geni yok eden bir delesyon **yoktur**.

Uc bagimsiz kaynagin ayni polimorfizmde bulustugu nokta: chr6:17.774.223-17.775.574 (1.351 bp) delesyonu DECIPHER populasyon setinde **58/1.892 = %3,07** frekansta ("1G del genotyped") gorulur. Bu, gnomAD'daki DEL_chr6_81c144f1 (AC=43, AF=3,41e-04, 4 homozigot, "lof" etiketli) ve DGV'deki esv3608224 (1000 Genomes Faz 3, 5/2.504) ile ayni varyanttir. Uc kaynak birlikte, gnomAD'in bu varyanta verdigi **"lof" etiketinin yanlis yonlendirici** oldugunu teyit eder.

---

## EK-2. DDD (Deciphering Developmental Disorders)

### EK-2a. DDG2P listesi

Cekilen kaynak: https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2026_08_28/DDG2P_2026-08-28.csv.gz (en guncel yayim, 634.487 bayt, **2.869 kayit**).

**KIF13A DDG2P listesinde YOKTUR.** Dosya genelinde buyuk/kucuk harf duyarsiz `KIF13A` aramasi sifir sonuc verir.

DDG2P'de bulunan kinesin genleri (12 tane): KIF11, KIF14, KIF1A, KIF22, KIF2A, KIF3B, KIF4A, KIF5A, KIF5B, KIF5C, KIF7, KIFBP. KIF13A bunlarin arasinda degildir.

**Paralog kalibrasyonu acisindan kritik olan KIF1A kayitlari** (birebir alanlar):

| G2P id | Hastalik | Allelik gereklilik | Guven | Molekuler mekanizma |
|---|---|---|---|---|
| G2P00959 | KIF1A-related NESCAV syndrome | **monoallelic_autosomal** | definitive | **altered gene product structure** |
| G2P00963 | KIF1A-related neuropathy, hereditary sensory | **biallelic_autosomal** | definitive | **absent gene product / loss of function** |

Yani iyi karakterize edilmis paralogda **tek kopya (monoallelic) hastaligi LoF ile degil, "degismis gen urunu yapisi" ile**, yani missense veya dominant-negatif ile olusur; **gercek fonksiyon kaybi ise ciftallelik (biallelic)** gerektirir. Bu, ilk rapordaki ClinGen KIF1A HI skoru 0 bulgusunun bagimsiz dogrulamasidir ve kinesin-3 ailesinde haploinsufficiency'nin **varsayilan mekanizma olmadigini** gosterir.

### EK-2b. DDD arastirma varyantlari

Cekilen kaynak: https://www.deciphergenomics.org/gene/KIF13A/ddd-research-variant-overlap (veri uc noktasi `https://www.deciphergenomics.org/data/gene/ENSG00000137177/dddrv-overlap`).

DDD'nin KIF13A'daki **tum** arastirma varyantlari, **2 tane**, ikisi de dizi varyanti:

| Konum (GRCh38) | Degisim | Sonuc | **Kalitim** | **Patojenite** | Genotip | Fenotip (DDD ust duzey terimler) |
|---|---|---|---|---|---|---|
| 6:17.764.129 | C>T | **stop_gained (Likely LOF)** | **Maternally inherited** | **Uncertain** | Heterozygous | Abnormality of head or neck; limbs; breast; eye; integument; musculoskeletal system; **nervous system** |
| 6:17.787.791 | C>CT | **frameshift_variant (Likely LOF)** | **Maternally inherited** | **Uncertain** | Heterozygous | Abnormality of head or neck; limbs; eye; musculoskeletal system |

**Bu, ek bolumun en sert bulgusudur.** DDD, siddetli gelisimsel bozuklugu olan cocuklarin trio bazli calismasidir ve tasarimi geregi **de novo** varyantlari yakalamaya ayarlidir. KIF13A'da buldugu iki "Likely LOF" varyantin **ikisi de de novo degil, anneden kalitilmistir** ve DDD kendi degerlendirmesinde ikisini de **patojenik degil, "Uncertain"** olarak birakmistir. Dominant, tam penetransli bir haploinsufficiency modeli altinda bu beklenen sonuc degildir.

Cekince, adil olmak icin: maternal kalitim tek basina zararsizligi kanitlamaz; azalmis penetrans, degisken ifade veya annede hafif etkilenmislik mumkundur. DDD kayitlarinda anne fenotipi verilmemistir. Ancak "anneden kalitilmis artı Uncertain" kombinasyonu, "de novo artı Pathogenic" kombinasyonunun **tersi** yonde kanittir.

---

## EK-3. dbVar nstd102: vaka mi, kontrol mu?

Cekilen kaynaklar:
- https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd102/
- https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz

Calisma basligi: "**nstd102 - Clinical Structural Variants**". Sayfadaki alan: **`Study Type: Case-Set`**.

Yani obsolete nstd101'in ve nstd37'nin yerine gecen nstd102 de **tasarimi geregi kontrol kolu tasimayan bir vaka setidir**. Bu, koordinatorun sordugu "vaka karsisinda kontrol sayilari" sorusunun cevabini kesinlestirir: **boyle bir sayi yoktur ve olamaz**.

KIF13A ile ortusen 21 nstd102 cagrisinin ilgili alanlarinin tamami:

| Alan | Dagilim (21 cagri) |
|---|---|
| `sample_id` | **21/21 bos** |
| `sampleset_id` | 21/21 = "1" (tek, ayrimsiz kume) |
| `experiment_id` | 21/21 = "1" |
| `validation` | **21/21 bos** |
| `zygosity` | **21/21 bos** |
| `clinical_source` | 21/21 = "ClinVar" |
| `origin` | 9 "not provided", 8 "unknown", 2 "de novo", 1 "germline", 1 "maternal" |
| `clinical_significance` | 15 Pathogenic, 6 Uncertain significance |
| `phenotype` | 15 "See cases", 5 "MedGen:C3661900,OMIM:104250.0001,OMIM:109630.0001", 1 "MedGen:CN517202" |

Ornek veya kontrol ayrimi yapabilecek her alan (`sample_id`, `sampleset_id` cesitliligi, `validation`, `zygosity`) **bostur veya sabittir**. nstd102, ClinVar'in klinik gonderimlerinin duz bir aynasidir; icinde saglikli birey yoktur.

---

## EK-4. gnomAD alt kohortlari

Cekilen kaynak: `POST https://gnomad.broadinstitute.org/api` (GraphQL sema introspection ve bolge sorgulari).

### EK-4a. v4'te alt kohort yoktur

GraphQL `StructuralVariantDatasetId` enum'unun **tam icerigi**:

```
gnomad_sv_r2_1
gnomad_sv_r2_1_controls
gnomad_sv_r2_1_non_neuro
gnomad_sv_r4
```

**gnomAD SV v4 icin `controls` veya `non_neuro` alt kumesi yoktur.** Dolayisiyla ilk raporda bulunan **436229__DEL** (448 kb, ekzom CNV) ve **DEL_chr6_33713a79** (164 kb, genom SV) tasiyicilarinin norogelisimsel fenotipleri disarida birakan bir kohortta olup olmadigi **API uzerinden belirlenemez**. Bu, o bulgunun ne kadar ileri goturulebilecegine dair gercek bir sinirdir ve boyle raporlanmalidir: gnomAD v4 siddetli pediatrik hastaliklari disarida birakmayi hedefler, ama tasiyici duzeyinde bunu dogrulayamiyoruz.

### EK-4b. v2.1 alt kohortlari (GRCh37, KIF13A 17.759.414-17.987.854)

| Veri seti | Donen SV | Gorulen max AN | Yaklasik birey | KIF13A ile ortusen DEL | **En buyuk ortusme** |
|---|---|---|---|---|---|
| gnomad_sv_r2_1 | 57 | 21.694 | ~10.847 | 13 | 5.000 bp (**genin %2,2'si**) |
| **gnomad_sv_r2_1_controls** | 34 | 10.384 | ~5.192 | 8 | 2.244 bp (**genin %1,0'i**) |
| **gnomad_sv_r2_1_non_neuro** | 48 | 16.684 | ~8.342 | 9 | 5.000 bp (**genin %2,2'si**) |

Uc veri setinin **hicbirinde** KIF13A'yi yok eden bir delesyon yoktur. En buyugu 5.000 bp'dir.

**Ancak bu null sonuc kanit degildir, guc yetersizligidir.** v4'te olculen orana gore (yaklasik 527.343 bireyde 3 olay, yani yaklasik 1/175.000), beklenen sayilar:

| Veri seti | Birey | Beklenen KIF13A-yikici delesyon sayisi |
|---|---|---|
| gnomad_sv_r2_1 | ~10.847 | **0,06** |
| gnomad_sv_r2_1_controls | ~5.192 | **0,03** |
| gnomad_sv_r2_1_non_neuro | ~8.342 | **0,05** |

Yani sifir gozlemek, varyant tamamen zararsiz olsa da tamamen patojenik olsa da beklenen sonuctur. **v2.1 kontrol ve non-neuro setleri bu soruyu cozemez** ve iddia lehine kullanilamaz.

### EK-4c. DEL_chr6_33713a79 cagri kalitesi ve tasiyici

Tekil varyant sorgusu (`structural_variant(dataset: gnomad_sv_r4, variantId: "DEL_CHR6_33713A79")`):

| Alan | Deger |
|---|---|
| pos - end | 17.820.000 - 17.984.000 (164.000 bp) |
| filters | **[] (PASS)** |
| **algorithms** | **['depth']** |
| **evidence** | **['BAF', 'RD']** |
| qual | 67 |
| consequences | `lof`, gen: **KIF13A** |
| AC / AN / AF | 1 / 126.090 / 8e-06, homozigot 0 |
| Tasiyici populasyonu | **afr**: AC=1 / AN=33.816. Diger tum populasyonlarda AC=0 (nfe 0/59.088, amr 0/12.594, sas 0/4.406, eas 0/4.054, fin 0/6.474, asj 0/3.180, ami 0/842, mid 0/64, rmi 0/1.572) |
| Tasiyici cinsiyeti | **XX** (AC=1/65.800 XX, 0/60.290 XY) |

Ilk raporda yuvarlak koordinatlar nedeniyle one surulen cekince **dogrulandi**: cagri yalnizca derinlik tabanlidir (`algorithms: ['depth']`), yani kirilma noktalari bin sinirlarina yuvarlanmistir. Ancak kanit tipi yalnizca okuma derinligi degildir; **BAF (B-allel frekansi)** destegi de vardir, ki bu gercek bir kopya sayisi degisimi icin bagimsiz bir isarettir. Filtreden gecmistir ve qual 67'dir. Ozetle: kirilma noktasi koordinatlari yaklasiktir, ama delesyonun varligi makul olcude desteklenmektedir.

---

## EK-5. Ek bolumun kanit tablosuna net etkisi

Ilk raporun 5. bolumundeki "ClinVar'da KIF13A-tek-gen delesyonu yalnizca bir tanedir" ifadesi **artik genisletilmelidir**. Tum kaynaklar birlestirildiginde, KIF13A'nin tek protein kodlayan gen oldugu bilinen delesyonlar:

| Kaynak | Kayit | Uzunluk | Genin %'si | Siniflama | **Kalitim** | Fenotip |
|---|---|---|---|---|---|---|
| ClinVar | VCV000442594 | 253.376 bp | %100 | Uncertain significance | **unknown** | Delayed speech and language development |
| DECIPHER | hasta 451718 | 146.274 bp | %63,6 | uncertain | **unknown** | Autistic behavior |
| DECIPHER | hasta 504413 | 182.672 bp | %73,0 | **likely_benign** | **maternal** | (yok) |
| gnomAD v4 ekzom | 436229__DEL | 448.271 bp | %99,2 | (siniflama yok) | (bilinmiyor) | (kohort bireyi) |
| gnomAD v4 genom | DEL_chr6_33713a79 | 164.000 bp | %71,3 | (siniflama yok) | (bilinmiyor) | (kohort bireyi) |

Buna DDD'nin iki maternal "Likely LOF" dizi varyanti eklenir.

Yani KIF13A-tek-gen fonksiyon kaybi olaylari **nadir ama tek ornek degildir**; bilinen kalitim durumlarinin **tamami maternaldir**; ve bugune kadar hicbir kaynak bunlardan herhangi birini **patojenik** olarak siniflamamistir.

### Guncellenmis "izin verilen en guclu cumle"

> "KIF13A'yi tek protein kodlayan gen olarak silen heterozigot delesyonlar hem etkilenmis bireylerde (ClinVar VCV000442594, DECIPHER hasta 451718) hem de siniflandirilmamis veya zararsiz kabul edilen bireylerde (DECIPHER hasta 504413, anneden kalitilmis, likely benign; gnomAD v4'te en az iki tasiyici) gozlenmistir. Bilinen kalitim durumlarinin tamami maternaldir, hicbiri patojenik olarak siniflanmamistir, KIF13A ne ClinGen doz duyarliliginda ne DDG2P'de ne GenCC'de yer alir, ve DECIPHER'in %HI tahmini geni %42,3 persentiline koyar. Bu kayitlar KIF13A'yi 6p22.3 delesyonlarinda bir aday olarak tutar; haploinsufficiency'yi kurmaz."

### Guncellenmis "asiri iddia olacak cumle"

> "DECIPHER ve DDD, KIF13A haploinsufficiency'sini bagimsiz olarak dogrular."

Tam tersi gecerlidir: DDD'nin iki KIF13A LoF varyantinin ikisi de maternal ve Uncertain'dir, DECIPHER'in tek maternal KIF13A-tek-gen delesyonu likely benign'dir, ve KIF13A DDG2P'de yoktur.

---

## EK bolumun cekilen URL'leri

| Amac | URL |
|---|---|
| DECIPHER gen sayfasi (tarayici) | https://www.deciphergenomics.org/gene/KIF13A/overview/clinical-info |
| DECIPHER hasta CNV JSON | https://www.deciphergenomics.org/data/gene/ENSG00000137177/patient-overlap |
| DECIPHER sayimlar | https://www.deciphergenomics.org/data/gene/ENSG00000137177/counts |
| DECIPHER DDD arastirma varyantlari | https://www.deciphergenomics.org/gene/KIF13A/ddd-research-variant-overlap ve `/data/gene/ENSG00000137177/dddrv-overlap` |
| DECIPHER %HI | https://www.deciphergenomics.org/files/downloads/HI_Predictions_Version3.bed.gz |
| DECIPHER populasyon CNV | https://www.deciphergenomics.org/files/downloads/population_cnv_grch38.txt.gz |
| DDG2P | https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2026_08_28/DDG2P_2026-08-28.csv.gz |
| G2P yayim dizini | https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/ |
| dbVar nstd102 calisma sayfasi | https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd102/ |
| dbVar nstd102 cagri tablosu | https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz |
| gnomAD GraphQL (sema, v2.1 alt kumeler, v4 varyant detayi) | `POST https://gnomad.broadinstitute.org/api` |

Bu ek bolumde **hicbir yeni PMID** kullanilmamistir; tum ifadeler yukaridaki veritabani uc noktalarindan canli cekilmistir.
