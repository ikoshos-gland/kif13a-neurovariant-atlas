# NSD1 -> NSD3: yöntemin örneklem dışı doğrulaması

Tarih: 17 Eyl 2026. Betik: `scripts/nsd_pair_calibration.py`. Veri: `data/interim/nsd_pair/`,
aday listesi `data/processed/NSD3_candidates_20260917.tsv`, her indirmenin kaydı
`data/raw/PROVENANCE.jsonl` ve dosya başına `.prov.json` yan dosyası.

Sınanan iddia: çok çiftli tarama (`scripts/paralog_pairs_scan.py`, `docs/olcekleme_fizibilite.md`)
paralog aktarımının genel bir özellik **olmadığını** buldu (3.198 çiftin 575'i dört zayıflık
filtresini geçiyor, kazanç medyanı 0,048, %35,5'inde kazanç sıfır ya da negatif) ama "ucuz ölçütlerle
hangi çiftte işe yarayacağı önceden kestirilebilir" dedi. Yayınlanabilir olan ikinci iddiadır ve
örneklem dışı hiç sınanmamıştı. Bu belge, taramanın 1. sıradaki nörogelişimsel çifti olan NSD1/NSD3
üzerinde tüm boru hattını çalıştırıp ölçülen aktarımı öngörülenle karşılaştırır.

Körlüğün tam kapsamı (abartmamak için): çalışmaya başlarken bilinen tek şey, çiftin taramanın en iyi
nörogelişimsel adayı olduğuydu (`HANDOFF.md` bölüm 4.1 bu beş çifti sayıyor). Bilinmeyenler şunlardı ve
1 ile 8 arasındaki adımlar bunlara bakılmadan yürütüldü: öngörülen kazanç, aktarım sinyali, pencere
arka planı, sıcak nokta penceresi, zayıflık bayrakları, fizibilite puanı, kesin sıra ve taramanın
beklediği negatif sayısı. `docs/olcekleme_fizibilite.md` dosyasındaki NSD1/NSD3 satırı ile
`paralog_scan_ranked_20260916.tsv` ancak 9. adımda açıldı. Hiçbir eşik, hiçbir domain sınırı, hiçbir
set tanımı sonuca göre ayarlanmadı; tüm eşik ve tier kuralları KIF1A/KIF13A vakasından birebir
kopyalandı (`scripts/kif1a_clinvar_curation.py`, `scripts/kif1a_build_calibration_sets.py`,
`scripts/evaluate_calibration_sets.py`, `scripts/build_kif13a_candidates.py`). Paylaşılan modülün
sınırı (NSD1 1539-2210), tek bir etiket yerleştirilmeden önce UniProt ve InterPro anotasyonundan
sabitlendi; sonuca bağımlı olmadığı 7.1'de ayrıca sınandı.

---

## 0. Hüküm

**Öngörü tuttu.** Tarama NSD1/NSD3'ü 3.198 çiftin 1. sırasına koymuş, fizibilite puanı 0,931 vermiş ve
sıfır zayıflık bayrağı bulmuştu; KIF1A/KIF13A ise 152. sırada, puan 0,7305. Ölçüm bu sıralamayı doğru
yönde doğruluyor: NSD1'in paylaşılan katalitik modülü içinde AlphaMissense AUROC **0,973**
(%95 GA 0,942-0,995, %90 özgüllükte duyarlılık 0,95, gen-özgü eşik 0,9732), KIF1A'nın motor
domainindeki 0,925 (0,878-0,964, duyarlılık 0,79, eşik 0,968) değerinin üzerinde; aynı yön, etiketle
eğitilmemiş üç bağımsız modelde de görülüyor (CADD 0,781 -> 0,820, ESM1b 0,851 -> 0,902, PrimateAI
0,801 -> 0,857), yani sonuç AlphaMissense'in popülasyon-frekansı eğitiminden gelen döngüsellikle
açıklanamaz. Taramanın ucuz ölçütleri de birebir yeniden üretildi (aktarım sinyali öngörü 0,914 /
ölçüm 0,915, pencere arka planı 0,664 / 0,664, kazanç 0,250 / 0,251, boşluk 0,010 / 0,009, sıcak nokta
penceresi 1968-2219 / 1968-2219, pencere içinde AC >= 3 olan negatif sayısı öngörü 31 / ölçüm 31).
Buna karşılık öngörünün **sıralamada kullandığı** niceliğin kendisi, yani kazanç (lift), beklenen şeyi
vermedi: kazancı iki katı olan NSD1/NSD3'te paralog-özdeşlik referans çizgisi 0,740 çıktı ve bu
KIF1A/KIF13A'nın 0,744'ünden ayrılmıyor. Yani "yüksek kazanç = daha güçlü paralog kanıtı" okuması
yanlış; işe yarayan şey paketin bütünü (iyi incelenmiş 125 etiket, korunmuş tek bir katalitik modül,
kısıtlı kardeş) oldu. En önemli sınırlama sayıdır: bu, örneklem dışı **tek** noktadır ve taramanın
"başarısız olacak" dediği bir çift hâlâ hiç çalıştırılmadı, dolayısıyla kestirim iddiası şu an yalnız
üst uçta desteklenmiş durumdadır.

---

## 1. Kimlik (canlı doğrulandı, 17 Eyl 2026)

| Gen | Ensembl geni | MANE transkript | RefSeq | Ensembl proteini | UniProt | Uzunluk |
|---|---|---|---|---|---|---|
| NSD1 (öğretmen) | ENSG00000165671, chr5:177131787-177300213, + | **ENST00000439151.7** | NM_022455.5 | ENSP00000395929 | Q96L73 | **2696 aa** |
| NSD3 (öğrenci) | ENSG00000147548, chr8:38269704-38383021, - | **ENST00000317025.13** | NM_023034.2 | ENSP00000313983 | Q9BZ95 | **1437 aa** |

Kaynaklar: gnomAD v4.1 GraphQL gene sorgusunun `mane_select_transcript` bloğu
(`https://gnomad.broadinstitute.org/api`, `data/interim/nsd_pair/NSD1_NSD3_gnomad_gene_20260917.json`);
Ensembl REST `lookup/symbol/homo_sapiens/NSD1` ve `lookup/id/ENST00000439151?expand=1`
(`https://rest.ensembl.org`); protein dizileri `rest.ensembl.org/sequence/id/<ENST>?type=protein`
(`data/raw/sequences/NSD1_ENST00000439151.fasta`, `NSD3_ENST00000317025.fasta`); UniProt kayıtları
`https://rest.uniprot.org/uniprotkb/Q96L73.json` ve `Q9BZ95.json` (`data/raw/uniprot/`).
Rapor: `data/interim/nsd_pair/NSD1_NSD3_identity_20260917.json`.

**İzoform tuzağı: bu çiftte yok.** Her iki gende MANE Select = Ensembl kanonik transkript, ve UniProt
kanonik dizisi ile MANE çevirisi **karakter karakter aynı** (NSD1 2696 = 2696, sha256 öneki
`3050432b8784dfca`; NSD3 1437 = 1437, `1740a00df4b91dc1`). Dolayısıyla `config/targets.yaml` içinde
her iki gen için `isoform_offsets: []` yazılıdır; KIF1A'daki +9 / +101 kayması burada yoktur ve
AlphaMissense'in protein başına tablosu ile VEP eklentisi aynı numaralandırmayı kullanır.

Önceki proje çalışmasında NSD1 için görülen ENST00000439151 **doğrulandı** (sürüm 7). Yine de iki
tuzak kayda geçirildi:

- NSD3'ün eski sembolü WHSC1L1'dir ve Ensembl `lookup/symbol` bu sembole 400 döner; sorgular NSD3 ile
  yapılmalıdır.
- NSD3'ün literatürde çok çalışılan **kısa izoformu** Q9BZ95-3 (NM_060248 / NP_060248.2, 645 aa) SET
  domaini hiç taşımaz. Kısa izoform üzerinden verilmiş bir kalıntı numarası MANE ile karşılaştırılamaz.

Her iki gen `config/targets.yaml` dosyasına doğrulanmış kimlikleriyle eklendi.

Kısıtlılık (aynı gnomAD çekimi): NSD1 LOEUF 0,149, pLI 1,00, missense Z 5,80, o/e missense 0,762;
NSD3 LOEUF 0,210, pLI 1,00, missense Z 5,93, o/e missense 0,667. İki gen de fonksiyon kaybına
dayanıksız; kardeş (NSD3) taramadaki en kısıtlı kardeştir ve bu yönüyle KIF13A'dan (LOEUF 0,449,
missense Z 2,79) belirgin biçimde daha elverişlidir.

## 2. Domain mimarisi ve **domain başına** özdeşlik

UniProt öznitelikleri (Q96L73, Q9BZ95) ve InterPro girdileri
(`https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/<acc>/`) canlı çekildi:
`data/interim/nsd_pair/NSD1_NSD3_domains_20260917.tsv` (183 satır),
`data/raw/interpro/Q96L73_domains_20260917.json`, `Q9BZ95_domains_20260917.json`.
İki protein birebir aynı modül dizilimini taşıyor: PWWP1, beş PHD tipi çinko parmağı, PWWP2, AWS,
SET, Post-SET ve NSD ailesine özgü Cys-His zengin domain.

Özdeşlik **tek bir küresel sayı olarak verilmemelidir**: hizalama 1.319 kalıntı çifti veriyor, küresel
özdeşlik 0,507 (kısa proteine göre 0,466) ve NSD1'in **1.377 kalıntısının NSD3'te karşılığı yok**
(çoğu 409-1538 arası ara bölge ile 2211-2696 arası C ucu).

| Modül | NSD1 | NSD3 | hizalanan | özdeşlik (hizalananlar üzerinden) |
|---|---|---|---|---|
| PWWP1 | 321-408 | 268-349 | 82/88 | 0,463 |
| PHD1 | 1539-1589 | 697-747 | 51/51 | 0,569 |
| PHD2 | 1590-1638 | 748-797 | 49/49 | 0,592 |
| PHD3 | 1639-1706 | 798-851 | 64/68 | 0,484 |
| PHD4 | 1707-1751 | 852-955 | 45/45 | 0,711 |
| PWWP2 | 1752-1850 | 956-1048 | 98/99 | 0,643 |
| AWS | 1890-1941 | 1093-1144 | 52/52 | **0,808** |
| **SET** | **1942-2065** | **1145-1268** | **124/124** | **0,790** |
| Post-SET | 2066-2082 | 1269-1285 | 17/17 | 0,765 |
| PHD5 (atipik) | 2118-2160 | 1321-1363 | 43/43 | 0,674 |
| CysHis | 2161-2210 | 1364-1413 | 49/50 | 0,612 |
| N ucu | 1-320 | 1-267 | 223/320 | 0,327 |
| ara bölge | 409-1538 | 350-696 | 339/1130 | 0,366 |
| C ucu | 2211-2696 | 1414-1437 | 9/486 | (hizalanmıyor) |

Dosya: `data/interim/nsd_pair/NSD1_NSD3_domain_identity_20260917.tsv`.

**Patojenik missense'in kümelendiği domain** (5. adımın kürasyonundan, yani ölçümden): 125 çekirdek
pozitifin **120'si (%96) PHD1'den CysHis'e uzanan bitişik katalitik modülün içinde**, dağılım
SET 36, PHD5 18, AWS 16, CysHis 15, PHD2 13, PHD4 8, PHD3 8, PWWP2 1, PHD1 2, AWS öncesi 2,
Post-SET 1; modül dışında yalnız 5 (PWWP1 2, C ucu 3). Bu nedenle **NSD1 1539-2210** (672 aa) bu
vakada KIF1A'nın motor domaininin karşılığıdır ve bundan sonraki her şey onun içinde olur. NSD3
karşılığı 697-1413'tür ve modül içi özdeşlik 430/666 = **0,646**, yani KIF1A/KIF13A motor domaininin
%62-70 aralığıyla aynı mertebede.

Haritanın dönüm noktası denetimi: UniProt'un NSD1 için işaretlediği 11 S-adenozil-L-metiyonin bağlama
kalıntısının **10'u** NSD3'te özdeş kalıntıya düşüyor (R1952/G1953/W1954 -> R1155/G1156/W1157,
T1994-Y1997 -> T1197-Y1200, N2020/H2021 -> N1223/H1224, N2065 -> N1268); tek uyuşmazlık Post-SET
içindeki K2071 -> H1274, yani korunumlu bir bazik kalıntı değişimi.
Dosya: `data/interim/nsd_pair/NSD1_NSD3_landmark_check_20260917.tsv`.

## 3. Varyant evrenleri

gnomAD v4.1, **transkript sorgusu** (gen sorgusu değil), `https://gnomad.broadinstitute.org/api`:

| | NSD1 (ENST00000439151.7) | NSD3 (ENST00000317025.13) |
|---|---|---|
| transkriptteki varyant | 6.931 | 4.772 |
| missense | 2.789 | 1.336 |
| **güvenilir gözlem (AC >= 1, AC0 filtresi yok)** | **2.721** | **1.250** |
| AC0 filtreli missense | 64 | 85 |
| farklı missense pozisyonu (güvenilir) | 1.679 | 843 |
| sinonim | 1.227 | 597 |

Dosyalar: `data/interim/nsd_pair/NSD1_ENST00000439151_gnomad_r4_all_20260917.tsv`,
`NSD3_ENST00000317025_gnomad_r4_all_20260917.tsv`, özet `NSD1_NSD3_universe_20260917.json`.

ClinVar (NCBI E-utilities, `esearch.fcgi` + `esummary.fcgi`, `<gen>[gene] AND single_gene[prop]`,
sınıflandırma alanı `germline_classification`):

| | NSD1 | NSD3 |
|---|---|---|
| tek gen kaydı | 2.694 | 164 |
| missense kaydı | 1.467 | 136 |
| **P/LP missense** | **158** | **0** |
| B/LB missense | 241 | 7 |
| VUS missense | 874 | 120 |
| MANE olmayan transkript başlıklı kayıt | 0 | 0 |

Dosyalar: `data/interim/nsd_pair/NSD1_clinvar_20260917.tsv`, `NSD3_clinvar_20260917.tsv`.
NSD3 tam olarak gereken öğrenci geni: **hiç P/LP missense'i yok**, yani araçların o gende çalışıp
çalışmadığı doğrudan ölçülemez. (Not: fizibilite belgesi NSD3 için "247 alel ve 1 P/LP alel" diyor;
fark `single_gene[prop]` filtresinden gelir. Bu çekimde 164 tek gen kaydının hiçbiri P/LP değil,
sınıflandırma dağılımı VUS 120, Likely benign 10, Benign 7, sınıflandırması boş 27.)

**Tuzak 14 bağımsız olarak yeniden üretildi.** Aynı gen ve missense filtreleriyle:

- doğru biçim, `("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`: **158**,
  bizim ayrıştırdığımız sayıyla birebir aynı;
- `pathogenic[clinsig]` biçimi: **371**, yani %135 şişmiş. Entrez bu alanı sessizce
  `pathogenic[All Fields]` olarak çeviriyor.

Her iki sayı da `stage universe` tarafından canlı sorgulandı ve
`NSD1_NSD3_universe_20260917.json` içinde `esearch_plp_missense_properties` /
`esearch_plp_missense_wrong_field` alanlarında kayıtlıdır.

## 4. Paralog kalıntı haritası

`atlas_mcp.structure.paralog_map` ile **aynı** hizalama parametreleri (Biopython PairwiseAligner,
global, BLOSUM62, boşluk -10/-0,5, uç boşlukları serbest), böylece iki vaka çalışması
karşılaştırılabilir kalıyor: `scripts/nsd_pair_calibration.py::align_pair`.
Çıktı `data/interim/paralog_maps/NSD1_to_NSD3_pairwise_MANE.tsv` (1.319 satır, hedef pozisyon,
öğretmen pozisyonu ve özdeşlik bayrağı).

- hizalanan çift: 1.319 (NSD3'ün 1.437 kalıntısının %92'si)
- özdeş: 669 (hizalananların 0,507'si)
- paylaşılan modül içinde (NSD1 1539-2210): 666 hizalanan, 430 özdeş, oran **0,646**
- NSD1'de karşılığı olmayan kalıntı: 1.377

## 5. NSD1 üzerinde etiket kürasyonu

Kurallar `docs/KIF1A_ClinVar_kurasyonu.md` dosyasından değiştirilmeden alındı; karar ağacı
`scripts/nsd_pair_calibration.py::assign_tier` içinde ve `scripts/kif1a_clinvar_curation.py`
ile aynıdır. 2.694 kaydın 1.121'i aday, her biri için VCV XML çekildi
(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`, `data/raw/clinvar/vcv_xml/NSD1_vcv_*`),
tek bir Ensembl VEP çalışmasıyla MANE protein pozisyonu ve skorlar eklendi.
Kürasyon tablosu: `data/interim/nsd_pair/NSD1_clinvar_curated_20260917.tsv`,
özet `NSD1_curation_summary_20260917.json`.

| Set | n | kalıntı | modül içi | NSD3'te özdeş kalıntı |
|---|---|---|---|---|
| `missense_positive_core` (80 gold + 45 silver) | **125** | 87 | **120** | 118 |
| `heldout_pp3_dependent_positive` (in siliko kanıta dayananlar) | **26** | 25 | 25 | 24 |
| `missense_negative_module_population` (gnomAD AC >= 3, modül içi) | **72** | 64 | 72 | 35 |
| `missense_negative_module_population_ac2` (duyarlılık) | 47 | 45 | 47 | 15 |
| `missense_negative_clinvar_outside_module` (ClinVar B/LB) | 184 | 172 | **6** | 23 |
| `splice_positive` | 46 | - | - | - |
| `reference_lof` (missense değil) | 452 | 134 | 51 | 43 |
| `sensitivity_conflict_P_vs_VUS` | 23 | 23 | 20 | 16 |

Setler: `data/interim/nsd_pair/NSD1_calibration_sets_20260917.tsv` (975 satır).

Pozitiflerin incelenme durumu kayıtlıdır: 95'i "criteria provided, single submitter" (1 yıldız),
30'u "criteria provided, multiple submitters, no conflicts" (2 yıldız); 50 Pathogenic, 59 Likely
pathogenic, 16 Pathogenic/Likely pathogenic. Yorumlarında 70'i de novo, 17'si fonksiyonel kanıttan
söz ediyor. Gönderenler dağılmış durumda (Labcorp/Invitae 60, University of Chicago 33, GeneDx 25,
Genome-Nilou 13), yani tuzak 26'nın uyardığı "tek gönderenin alışkanlığı" durumu yok. Fenotip alanı
ağırlıkla Sotos sendromu (69 kayıt) ve "not provided" (74).

**Döngüsellik ayrımı.** 158 P/LP missense adayının yorum metinlerinde 84'ü in siliko kanıta, 73'ü
de novo, 18'i fonksiyonel çalışmaya atıf yapıyor. Sınıflandırması PP3'e dayanmış olabilecek
**26 tek gönderenli LP ayrı tutuldu** ve ana metriğe girmedi.

**Negatifler pozitiflerle aynı yapısal bölgeden alındı.** ClinVar'ın B/LB missense setinde 184
varyantın yalnız **6'sı** paylaşılan modülün içinde (108 ara bölge, 41 C ucu, 27 N ucu), yani KIF1A'da
görülen konum karışıklığı burada da geçerli. Bu nedenle ana negatif set, modül içinde gnomAD v4.1'de
AC >= 3 ile görülen ve ClinVar'da P/LP veya çelişkili olarak etiketlenmemiş 72 missense varyanttır;
duyarlılık için AC = 2 olan 47 varyant ayrı tutuldu. Modül içindeki kontrol adayları alel sayısına
göre: AC1 210, AC2 47, AC3-4 43, AC >= 5 29. Her kontrolün referans aminoasidi MANE dizisine karşı
doğrulandı (uyuşmazlıkta betik hata verir).

Etiket sayısı yeterlidir: 125 varyant / 87 kalıntı, bunların 120'si / 82'si modül içinde. Yöntemin
"yaklaşık 30 veya daha fazla küre edilmiş patojenik missense" ölçütü rahatça karşılanıyor, yani
"kalibre edilemeyecek kadar az etiket" bulgusu bu çift için geçerli değil.

## 6. Skorlama

Tek bir Ensembl VEP çalışması (`https://rest.ensembl.org/vep/human/region`, parametreler
`atlas_mcp.variants.VEP_PARAMS`): AlphaMissense, CADD v1.7, dbNSFP ESM1b ve PrimateAI değerlendirme
hedefi; REVEL ve MetaRNN yalnız referans satırı (klinik etiketle eğitilmişler, ne hedef ne konsensüs
üyesi olabilirler). Girdi VCF biçiminde kromozom ve konuma göre sıralandı; dbNSFP'nin transkript
listeleri `Ensembl_transcriptid` üzerinden indekslendi (`atlas_mcp.variants._flatten_vep`).
Pozitifler ve negatifler **aynı çağrıda** skorlandı; skor kapsaması 477/477 (altı modelin hepsinde).
Çıktı: `data/interim/nsd_pair/NSD1_calibration_scores_20260917.tsv`,
ham VEP yanıtı `data/raw/vep/NSD1_calibration_scores_vep_20260917.raw.json.gz`.

## 7. Ölçüm

Ana karşılaştırma: paylaşılan modül içindeki 120 çekirdek pozitif (82 kalıntı) ile modül içindeki
72 gnomAD popülasyon kontrolü (64 kalıntı). Güven aralıkları **kalıntı-gruplu bootstrap** (2.000
tekrar, aynı kalıntının farklı değişimleri birlikte yeniden örneklenir).

| Model | AUROC | %95 GA | %90 özgüllükte duyarlılık | Gen-özgü eşik |
|---|---|---|---|---|
| REVEL (etiketle eğitilmiş, referans) | 0,975 | 0,954-0,991 | 0,92 | 0,755 |
| **AlphaMissense** | **0,973** | **0,942-0,995** | **0,95** | **0,9732** |
| MetaRNN (etiketle eğitilmiş, referans) | 0,972 | 0,945-0,991 | 0,89 | 0,811 |
| ESM1b | 0,902 | 0,845-0,954 | 0,66 | -12,64 |
| PrimateAI | 0,857 | 0,786-0,916 | 0,49 | 0,868 |
| CADD v1.7 | 0,820 | 0,738-0,889 | 0,40 | 32,0 |
| **Referans: kalıntı NSD3'te özdeş mi** | **0,740** | | | |
| **Referans: paylaşılan modülün içinde mi** (modül içi ölçümde tanım gereği) | **0,500** | | | |

AC >= 2 negatifleriyle (n = 119) tablo pratikte aynı: AlphaMissense 0,975 (0,955-0,992), ESM1b 0,915,
CADD 0,843, PrimateAI 0,840, özdeşlik referansı 0,773.

Alt kümeler:

| Karşılaştırma | n+ | AlphaMissense | CADD | ESM1b |
|---|---|---|---|---|
| Çekirdek pozitifler (ana metrik) | 120 | 0,973 | 0,820 | 0,902 |
| PP3'e dayalı, ayrı tutulan pozitifler | 25 | 0,968 | 0,799 | 0,917 |
| Fonksiyonel kanıtlı pozitifler | 17 | 0,969 | 0,877 | 0,945 |

PP3'e dayanan alt kümede AlphaMissense'in **daha yüksek çıkmaması** (0,968 < 0,973) bu çiftte
döngüsellik izinin KIF1A'daki kadar belirgin olmadığını gösteriyor (KIF1A'da 0,961 ve 0,925).

### 7.1 Bölge sınırı duyarlılığı

Paylaşılan modül sınırı sonuca göre seçilmediğini göstermek için aynı ölçüm iki daha dar bölgede
yinelendi (`stage measure` içindeki `SENSITIVITY` karşılaştırmaları):

| Bölge | n+ | n- | AlphaMissense | CADD | ESM1b | paralog özdeşliği |
|---|---|---|---|---|---|---|
| paylaşılan modül 1539-2210 (kullanılan) | 120 | 72 | **0,973** (0,942-0,995) | 0,820 | 0,902 | 0,740 |
| taramanın dar penceresi 1968-2219 | 64 | 25 | 0,949 (0,873-0,997) | 0,800 | 0,850 | 0,744 |
| yalnız SET + Post-SET 1942-2082 | 37 | 14 | 0,986 (0,952-1,000) | 0,854 | 0,863 | 0,714 |

Üç sınırda da AlphaMissense noktası KIF1A'nın 0,925'inin üzerinde ve paralog-özdeşlik referansı
0,71-0,74 aralığında sabit kalıyor. Dar pencerede güven aralığı beklendiği gibi genişliyor
(25 negatif), yani modül sınırını kullanmanın asıl kazancı istatistiksel güçtür, sonucun yönü değil.

### 7.2 Zorunlu referans çizgileri

1. **"Paylaşılan modülün içinde mi", tüm protein üzerinde.** ClinVar B/LB negatifleriyle yapılan
   karıştırılmış karşılaştırmada (125 pozitif, 184 negatif) bu tek kural **AUROC 0,964** alıyor;
   aynı karşılaştırmada AlphaMissense 0,992, PrimateAI 0,993, CADD 0,971, ESM1b 0,987. Yani konum
   karışıklığı KIF1A'da olduğu gibi burada da neredeyse her modelle yarışıyor (KIF1A'da 0,995) ve
   ana metriğin modül içi olması zorunludur.
2. **"Kalıntı NSD3'te özdeş mi".** Modül içi AUROC **0,740** (AC >= 2 setinde 0,773). Pozitiflerin
   %96,7'si, popülasyon kontrollerinin %48,6'sı özdeş kalıntıda.
3. **Bölgede gözlenen/beklenen sinonim.** Samocha trinükleotid hızları ve gnomAD v4.1 ekzomlarıyla,
   beklenti iki genin havuzlanmış sinonim bölgeleri üzerinde kalibre edilerek
   (`scripts/regional_constraint.py` yardımcıları, hız tablosu
   `data/raw/mutation_rates/samocha2014_trinucleotide_rates_20260916.txt`):

| | aa | o/e missense | %90 aralık | o/e sinonim (kalibrasyon kontrolü) |
|---|---|---|---|---|
| NSD1 tüm protein | 2696 | 0,879 | 0,850-0,908 | 1,026 |
| **NSD1 paylaşılan modül** | 672 | **0,433** | 0,395-0,475 | **1,052** |
| NSD1 modül dışı | 2024 | 1,032 | 0,996-1,069 | 1,017 |
| NSD1 SET | 124 | 0,322 | 0,252-0,416 | 0,969 |
| NSD3 tüm protein | 1437 | 0,767 | 0,731-0,805 | 0,950 |
| **NSD3 paylaşılan modül** | 717 | **0,760** | 0,710-0,814 | **0,921** |
| NSD3 modül dışı | 720 | 0,774 | 0,723-0,829 | 0,980 |
| NSD3 SET | 124 | 0,641 | 0,539-0,767 | 1,011 |

Tam tablo (16 bölge x 2 gen, sinonim-normalize edilmiş missense o/e dahil):
`data/interim/nsd_pair/NSD1_NSD3_regional_constraint_20260917.tsv`.

Sinonim o/e değerleri ana bölgelerde 0,92-1,05 aralığında, yani beklenti modeli bu bölgelerde
kalibre. (Uyarı: NSD1 PHD1 1,51 ve PHD5 1,43 gibi küçük ve CpG zengini bölgelerde sinonim o/e
şişiyor; o domainler için sinonim-normalize edilmiş sütun okunmalıdır.)

Bu tablonun yöntem açısından iki sonucu var. Birincisi olumlu: NSD1'de patojenik missense'in
kümelendiği modül, seçilimin de tam olarak yoğunlaştığı yerdir (0,433'e karşı modül dışında 1,032),
yani "mekanizma modülde" varsayımı öğretmen gende doğrulanıyor. İkincisi uyarıcı: **NSD3'te modülün
özel bir kısıtlılığı yok** (0,760'a karşı modül dışı 0,774); NSD3 gen boyunca eşit ölçüde mütevazı
kısıtlı. Karşılaştırma ölçeği KIF çiftiyle şaşırtıcı biçimde paralel: öğretmen modülü / öğrenci
modülü oranı NSD'de 0,433 / 0,760 = 1,75 kat, KIF'te 0,450 / 0,814 = 1,81 kat sığ.

## 8. NSD3'e aktarım ve katman tablosu

Aktarılan eşik: AlphaMissense 0,9732 (NSD1'de %90 özgüllük). Karşılaştırma için AlphaMissense'in
genel "likely pathogenic" sınırı 0,564. Katman kuralları `scripts/build_kif13a_candidates.py`
dosyasından değiştirilmeden alındı (A: eşiği geçiyor + kalıntı özdeş + NSD1'de o kalıntıda patojenik
missense var; B: eşiği geçiyor + özdeş + paylaşılan modülde; C: eşiği geçiyor, başka yerde;
D: eşiğin altında ama paralog desteği var).

| | tüm missense | güvenilir gözlemli |
|---|---|---|
| NSD3 gnomAD missense | 1.336 | 1.250 |
| **genel eşik 0,564 ile çağrılan** | **399** | 359 |
| **kalibre eşik 0,9732 ile çağrılan** | **101** | 89 |
| Katman A | 23 | 19 |
| Katman B | 45 | 43 |
| Katman C | 33 | 27 |
| Katman D | 30 | 22 |
| Katman yok | 1.205 | 1.139 |

Aday listesi: `data/processed/NSD3_candidates_20260917.tsv` (26 sütun, her gnomAD missense varyantı
için skorlar, paralog kanıtı ve katman). Özet: `data/interim/nsd_pair/NSD3_transfer_summary_20260917.json`.
Kalibre eşik, genel eşiğin çağırdığı varyantların dörtte birini bırakıyor; KIF13A'da aynı oran
295/982 idi.

Katman A'nın dağılımı SET 10, PHD2 5, AWS 4, PHD5 3, PHD3 1. En yüksek skorlu örnekler ağırlıkla
çinko koordine eden sistein kalıntıları ve SAM cebi kalıntıları: C794R (AM 0,9996; NSD1'de C1635G
patojenik), C1108R ve C1108Y (NSD1 C1905G), C1341R (NSD1 C2138S/C2138R), D1211V ve D1211E
(NSD1 D2008H), C1324S (NSD1 C2121G/C2121Y), G1176C (NSD1'de G1973R/G1973D/G1973V), W1157S
(NSD1 W1954C, yani SAM bağlama kalıntısı).

## 9. Öngörü ile karşılaştırma (tarama dosyası burada açıldı)

`docs/olcekleme_fizibilite.md` bölüm 4 ve 7.1 ile `data/interim/paralog_pairs/paralog_scan_ranked_20260916.tsv`
satır 1 (NSD1/NSD3) ve satır 152 (KIF1A/KIF13A).

### 9.1 Ucuz ölçütler: birebir yeniden üretildi

| Ölçüt | Taramanın öngörüsü | Bu çalışmanın ölçümü |
|---|---|---|
| küresel özdeşlik (hizalananlar) | 0,507 | 0,507 |
| aktarım sinyali (benzersiz P/LP kalıntısı, ham ClinVar) | 0,914 | **0,915** |
| pencere arka planı | 0,664 | **0,664** |
| **kazanç (lift)** | **0,250** | **0,251** |
| P/LP kalıntılarında boşluk oranı | 0,010 | 0,009 |
| sıcak nokta penceresi (yarısını içeren en dar pencere) | 1968-2219 | **1968-2219** |
| benzersiz P/LP kalıntısı | 105 | 106 |
| pencere içinde gnomAD AC >= 3 (negatif seti kurulabilir mi) | 31 | **31** |
| pencere içinde gnomAD AC >= 2 | 43 | **43** |
| sıcak nokta domaini | SET + Post-SET | SET en büyük tek domain (36/125), ama kümelenme PHD1-CysHis modülünün tamamına yayılıyor |

Dosya: `data/interim/nsd_pair/NSD1_NSD3_measured_lift_20260917.json`. Bu satırlar aynı kaynaklardan
(ClinVar, gnomAD, Ensembl) türetildiği için burada sınanan şey taramanın **aritmetiğinin doğruluğudur**,
kestirim gücü değil. Aritmetik doğru çıktı: tek fark benzersiz kalıntı sayısındadır (106'ya karşı 105) ve
taramanın 158 P/LP missense kaydının 157'sini MANE'ye yerleştirmesinden geliyor
(`n_plp_positions_on_mane` = 157); buradaki ayrıştırıcı 158'inin hepsini yerleştirdi. P/LP missense
kayıt sayısı iki çalışmada da 158.

Küre edilmiş sete geçildiğinde kazanç yükseliyor: modül içi 82 pozitif kalıntının 0,963'ü NSD3'te
özdeş, modül arka planı 0,646, yani **kazanç 0,317**; popülasyon kontrollerine karşı ise
0,963 - 0,486 = **0,477**.

### 9.2 Pahalı uç nokta: hiç ölçülmemişti, öngörülen yönde çıktı

| | NSD1/NSD3 | KIF1A/KIF13A |
|---|---|---|
| taramadaki sıra | **1 / 3.198** | 152 / 3.198 |
| fizibilite puanı | 0,931 | 0,7305 |
| zayıflık bayrağı | 0 | 0 |
| öngörülen kazanç | 0,250 | 0,128 |
| öngörülen aktarım sinyali | 0,914 | 0,826 |
| kardeş LOEUF / missense Z | 0,21 / 5,93 | 0,449 / 2,79 |
| **paylaşılan domain** | PHD1-CysHis modülü, 672 aa | motor domain + boyun bağlayıcı, 365 aa |
| **paylaşılan domain özdeşliği** | 0,646 (SET 0,790, AWS 0,808) | 0,62-0,70 |
| çekirdek pozitif (domain içi) | 120 (82 kalıntı) | 110 (70 kalıntı) |
| popülasyon negatifi AC >= 3 | 72 | 61 |
| **AlphaMissense AUROC** | **0,973 (0,942-0,995)** | 0,925 (0,878-0,964) |
| AlphaMissense %90 özg. duyarlılık | **0,95** | 0,79 |
| kalibre eşik | 0,9732 | 0,968 |
| CADD AUROC (en az döngüsel referans) | 0,820 (0,738-0,889) | 0,781 (0,696-0,858) |
| ESM1b AUROC | 0,902 | 0,851 |
| PrimateAI AUROC | 0,857 | 0,801 |
| **paralog-özdeşlik referansı** | **0,740** | **0,744** |
| "domain içinde mi" referansı (tüm protein) | 0,964 | 0,995 |
| PP3 ayrılmış set (AlphaMissense) | 0,968 (çekirdeğin altında) | 0,961 (çekirdeğin üstünde) |
| öğretmen domain o/e missense | 0,433 | 0,450 |
| **öğrenci domain o/e missense** | **0,760** (modül dışı 0,774) | 0,814 (kuyruk 0,927) |
| kardeşte ClinVar P/LP missense | 0 | 0 |
| aday listesi: gnomAD missense | 1.336 | 2.671 |
| genel eşikle çağrılan | 399 | 982 |
| kalibre eşikle çağrılan | 101 | 295 |
| katman A / B / C / D | 23 / 45 / 33 / 30 | 52 / 73 / 170 / 33 |

### 9.3 Öngörünün hangi kısmı tuttu, hangi kısmı tutmadı

**Tuttu:**

1. Sıralama yönü. 1. sıradaki çift, 152. sıradaki çiftten daha yüksek ölçülmüş ayrım gücü verdi ve bu
   dört modelde de aynı yönde (AlphaMissense +0,048, CADD +0,039, ESM1b +0,051, PrimateAI +0,056).
   Tek bir modelin kendine özgü davranışı değil.
2. "Kalibre edilebilir etiket var" öngörüsü. Tarama 158 P/LP missense ve pencere içinde 31 popülasyon
   kontrolü bekliyordu; kürasyon 125 kullanılabilir pozitif ve modül içinde 72 kontrol verdi.
3. "Sıcak nokta tek bir katalitik domain" öngörüsü. Pozitiflerin %96'sı bitişik tek bir modülde.
4. "Kardeş kısıtlı" öngörüsü, gen düzeyinde: NSD3 LOEUF 0,21 ve missense Z 5,93, taramadaki en
   kısıtlı kardeş; bu ölçüde KIF13A'dan gerçekten daha iyi.
5. Yöntemin kendi uyarıları da taşındı: konum karışıklığı yine neredeyse her modeli yakalıyor (0,964),
   yani ClinVar B/LB ile tüm protein ölçümü yapan bir yaklaşım burada da sahte temiz sonuç verirdi.

**Tutmadı ya da yanıltıcı çıktı:**

1. **Kazanç, paralog kanıtının gücünü kestirmiyor.** Taramanın sıralama ağırlığının %25'i kazançtan
   geliyor ve NSD1/NSD3'ün kazancı KIF1A/KIF13A'nın iki katı (0,250'ye karşı 0,128). Buna karşılık
   ölçülen paralog-özdeşlik AUROC'u iki çiftte ayırt edilemez (0,740 ve 0,744). Kazanç, taramada
   "paralog kanıtı bir şey katar mı" sorusunun göstergesi olarak kullanılıyordu; bu sınavda o rolü
   yerine getirmedi. Yüksek kazanç, modelin başarısını yükseltti (arka planın düşük olması pozitifleri
   daha ayırt edilebilir yapıyor olabilir) ama paralog özelliğinin kendi katkısını yükseltmedi.
2. **Kardeşin kısıtlılığı gen düzeyinde ölçülüyor, oysa aktarım bölge düzeyinde yapılıyor.** NSD3 gen
   düzeyinde taramanın en iyi kardeşi, ama bölgesel ölçüm paylaşılan modülün proteinin geri kalanından
   daha kısıtlı **olmadığını** gösteriyor (0,760'a karşı 0,774). Taramanın hiçbir ölçütü bunu
   göremiyordu; KIF13A'da aynı kör noktanın bulunması bir oturum sürmüştü.
3. "Sıcak nokta domaini SET + Post-SET" tahmini fazla dar. Pozitiflerin yalnız %29,6'sı SET +
   Post-SET aralığında; kalanı PHD parmakları, AWS ve CysHis'te. Pencere ölçütü domain yerine
   kullanıldığında negatif seti gereksiz daralıyor (31'e karşı 72).

## 10. Kestirim iddiasını ne çürütür, bu deney ne kadar yaklaştı

Çürütücü gözlemler, önceden yazılmış hâli:

1. **Doğrudan çürütme:** taramanın üst sırasındaki ve sıfır zayıflık bayraklı bir çiftte, paylaşılan
   domain içi AUROC'un etiketle eğitilmemiş modellerde 0,85'in altına düşmesi, ya da güven
   aralığının paralog-özdeşlik referans çizgisiyle (yaklaşık 0,74) örtüşmesi. Bu olmadı: AlphaMissense
   0,973 (alt sınır 0,942), CADD 0,820 (alt sınır 0,738).
2. **Kalibrasyonun anlamsızlaşması:** gen-özgü eşiğin genel eşikten ayrılmaması, yani kalibrasyonun
   hiçbir şey eklememesi. Bu da olmadı: 0,9732'ye karşı 0,564, çağrılan varyant sayısı 101'e karşı 399.
3. **Etiket yetersizliği:** üst sıradaki bir çiftte kürasyonun paylaşılan domain içinde 30'dan az
   kullanılabilir pozitif bırakması. Bu da olmadı (120).
4. **Sıralamanın ters dönmesi:** 152. sıradaki çiftin 1. sıradakinden daha yüksek ölçülmüş ayrım
   vermesi. Olmadı, ama fark küçüktür (0,925'e karşı 0,973) ve KIF1A'nın güven aralığı
   (0,878-0,964) ile NSD1'in aralığı (0,942-0,995) **örtüşüyor**. İki nokta arasındaki sıralama bu
   nedenle istatistiksel olarak kesin değil.

Deney iki noktada çürütmeye yaklaştı. Birincisi 9.3'teki birinci maddedir: kazanç ölçütü, kendisine
atfedilen işi yapmadı, yani "hangi ucuz ölçüt kestiriyor" sorusunun cevabı taramanın varsaydığından
farklıdır. İkincisi bölgesel kısıtlılıktır: NSD3'ün paylaşılan modülü özel olarak kısıtlı olmadığı
için, aktarılan eşiğin gerçek bir hastalık mekanizmasına oturduğu **gösterilmedi**; gösterilen şey
eşiğin öğretmende ölçülebilir ve öğrenciye uygulanabilir olduğudur. Bir çift, eşik iyi çalıştığı hâlde
kardeşinde hastalık olmadığı için yayın değeri taşımayabilir, ve tarama bu ayrımı yapmıyor.

Gerçek çürütme testi hâlâ yapılmadı: tarama iki çifti de üst %5'te sıralamıştı. Kestirim iddiası
"düşük sıralı ya da bayraklı bir çift gerçekten başarısız olur" kısmını içeriyor ve **o yön hiç
sınanmadı**. Bir yanlış-pozitif kontrolü olmadan gözlemler, "30'dan fazla küre edilmiş etiketi ve
korunmuş bir domaini olan her çift çalışır" şeklindeki daha zayıf ve daha az ilginç hipotezle de
uyumludur.

## 11. Yöntem makalesi için ne anlama geliyor

1. **Merkezi iddia ayakta, ama tek örneklem dışı noktayla.** "Paralog aktarımı genel bir özellik
   değil, ancak hangi çiftte işe yarayacağı ucuz ölçütlerle önceden kestirilebilir" cümlesinin ikinci
   yarısı, ilk kez örneklem dışı sınandı ve doğru yönde sonuç verdi. Makale bu deneyi n = 1 doğrulama
   olarak sunmak zorundadır; iki nokta bir eğri değildir ve iki AUROC güven aralığı örtüşüyor.
2. **Yayınlanabilir asıl ürün, tek bir sayı değil, kalibrasyon protokolüdür.** İki bağımsız çiftte
   aynı üç desen çıktı: (i) klinik benign setiyle tüm protein ölçümü konumu kazandırıyor (0,995 ve
   0,964), (ii) aynı bölgeden çekilmiş popülasyon negatifleriyle gerçek ayrım 0,92-0,97'ye iniyor,
   (iii) gen-özgü eşik genel eşikten çok uzağa düşüyor (0,968 ve 0,9732'ye karşı 0,564) ve aday
   sayısını üçte bire indiriyor. Bu üç desenin yeniden üretilmesi, tek bir genin aday listesinden
   daha sağlam bir katkıdır.
3. **Sıralama ölçütü revize edilmeli.** Kazanç (lift) ağırlığı, ölçülen paralog katkısını kestirmediği
   için savunulamaz durumda. Makale ya ağırlığı düşürmeli ya da kazancı "model performansı
   göstergesi" olarak yeniden tanımlamalıdır. Ayrıca kardeşin **bölgesel** kısıtlılığı taramaya
   eklenmelidir: iki vakada da kardeşin paylaşılan domaini, gen düzeyindeki kısıtlılığın vaat ettiği
   kadar kısıtlı değildi ve bu, aktarımın biyolojik zeminiyle ilgili tek ciddi zayıflıktır.
4. **Sonraki deney, başarısızlık öngörülen bir çift olmalıdır.** Örneğin 3-4 zayıflık bayraklı,
   kardeşi kısıtsız ve kazancı sıfır ya da negatif bir çift aynı boru hattından geçirilmelidir.
   Kestirim iddiası ancak orada düşük AUROC ölçülürse bir iddia olur. `docs/olcekleme_fizibilite.md`
   içindeki 1.737 `no_conservation_lift` bayraklı çift bu kontrolü seçmek için yeterli havuzdur.
5. **Maliyet öngörüsü tutuyor.** Fizibilite belgesi NSD1 için "1 ile 2 gün elle inceleme" demişti;
   bu çalışma, KIF1A vakasından kopyalanan kurallar ve hazır betikler sayesinde tek bir oturumda
   tamamlandı. Yöntemin ölçeklenebilir olduğu iddiası, en az bu ölçüde, ölçülmüş durumdadır.

## 12. Yeniden üretim

```powershell
python -u scripts\nsd_pair_calibration.py --stage identity
python -u scripts\nsd_pair_calibration.py --stage domains,map
python -u scripts\nsd_pair_calibration.py --stage universe
python -u scripts\nsd_pair_calibration.py --stage curate
python -u scripts\nsd_pair_calibration.py --stage score
python -u scripts\nsd_pair_calibration.py --stage measure
python -u scripts\nsd_pair_calibration.py --stage constraint
python -u scripts\nsd_pair_calibration.py --stage lift
python -u scripts\nsd_pair_calibration.py --stage transfer
```

Sistem `python` 3.13.14 kullanıldı; proje venv'i (`mcp\atlas_server\.venv`) bu makinede hâlâ Windows
tarafından bloklanıyor. Her aşama `--out` ile başka bir dizine yazdırılabilir, ilerleme ve kaynak
bilgisi stderr'e gider, her dosya `atlas_mcp.core` yardımcılarıyla yazıldığı için `.prov.json` yan
dosyası ve `data/raw/PROVENANCE.jsonl` satırı üretir.

Ensembl REST bu oturumda aralıklı olarak 500 ve 503 döndürdü; betik Ensembl çağrılarında 9 denemeye
kadar yeniden dener (`ens_get`). VEP toplu çağrıları 200'lük gruplar hâlinde gider ve sıralanmamış
girdi reddedildiği için girdi önceden sıralanır.
