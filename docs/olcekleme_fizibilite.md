# Ölçekleme fizibilitesi: KIF1A -> KIF13A yöntemi kaç paralog çiftinde işler?

Tarih: 16 Eyl 2026. Betik: `scripts/paralog_pairs_scan.py`.
Tablolar: `data/interim/paralog_pairs/paralog_scan_{rich_genes,pairs_all,pairs_qualified,metrics,ranked}_20260916.tsv`
(her birinin yanında `.prov.json`, ayrıca `data/raw/PROVENANCE.jsonl`).

Sorulan soru: atlasın yöntemi (klinik olarak iyi tanımlanmış bir gende varyant etki modellerini kalibre et,
sonra eşikleri ve paralog korunum kanıtını az çalışılmış kardeşine aktar) tek bir vaka incelemesi mi, yoksa
çok sayıda paralog çiftinde işleyen bir yöntem mi? Tam boru hattı hiçbir çiftte çalıştırılmadı; yalnızca ucuz
fizibilite ölçütleri hesaplandı.

## 1. Özet

- **3.510 paralog çifti** biçimsel ölçütü karşılıyor (zengin gende >= 30 ClinVar P/LP missense, kardeşte <= 2).
  Bu sayı **abartılıdır**: tek bir zengin gen çoğu zaman aynı ailedeki onlarca boş kardeşle eşleşiyor
  (USP9X 63, KRT1/KRT14/DES/GFAP/LMNA/KRT5 48'er, KIF1A 26). Anlamlı birim çift değil, **337 zengin gen**.
- Ucuz metrikler 3.198 çiftte hesaplandı. **Tipik çift zayıf**: patojenik kalıntıların kardeşte özdeş kalma oranı
  medyan 0,412; yerel arka plana göre kazanç (lift) medyan yalnızca **0,048**; çiftlerin %35,5'inde kazanç
  sıfır veya negatif; kardeş genlerin LOEUF medyanı 0,912, yani yarısı kısıtlanmış bile değil.
- Dört zayıflık ölçütünün hiçbirine takılmayan çift sayısı **575** (3.198'in %18'i), bunlar **183 farklı zengin gen**
  ve 384 farklı kardeşe dağılıyor.
- **KIF1A/KIF13A 3.198 çift içinde 152., zengin gen başına en iyi kardeş listesinde 337 içinde 71.** Yani yöntemin
  en iyi örneği değil, ama üst dörtte birinde ve dört zayıflık bayrağının hiçbirini almıyor.
- **Karar: çok çiftli çalışma yapılmaya değer, ama "paralog aktarımı genel olarak işler" iddiası kurulamaz.**
  Kurulabilecek iddia şudur: çiftlerin yaklaşık altıda birinde aktarım için gereken koşullar sağlanır ve bu
  koşullar tam boru hattı çalıştırılmadan, ucuz ölçütlerle önceden ayırt edilebilir. Negatif sonuç kısmı gerçektir
  ve saklanmamalıdır: çiftlerin çoğunda paralog kanıtı bir şey katmaz.

## 2. Tarama nasıl yapıldı

### 2.1 Kaynaklar ve tarihler

| Kaynak | Sürüm / tarih | Kullanım |
|---|---|---|
| ClinVar `gene_specific_summary.txt` | dosyanın kendi başlığı: 13 Eyl 2026; indirme 16 Eyl 2026 | gen başına P/LP alel sayısı, ön eleme |
| HGNC complete set | indirme 16 Eyl 2026 | simge <-> ENSG <-> Entrez eşlemesi, onaylı protein kodlayan filtresi |
| NCBI E-utilities (ClinVar esearch/esummary) | sorgu 16 Eyl 2026 | kesin P/LP missense sayısı ve protein konumları |
| Ensembl REST `/homology`, `/sequence/id` | sorgu 16 Eyl 2026 (GRCh38) | insan paralogları, MANE protein dizileri |
| gnomAD GraphQL, `gnomad_r4` (v4.1) | sorgu 16 Eyl 2026 | MANE Select transkripti, LOEUF / pLI / missense Z |
| UniProt REST (reviewed, 9606) | sorgu 16 Eyl 2026 | sıcak noktanın hangi domaine düştüğü |

92.938 ClinVar gen satırı, 19.254 onaylı protein kodlayan HGNC geni tarandı.

### 2.2 Aşamalar

1. **`sources`**: iki toplu dosya indirildi ve önbelleğe alındı.
2. **`rich`**: ClinVar özet tablosunda >= 30 P/LP aleli olan, HGNC'de onaylı protein kodlayan ve ENSG'si olan
   **2.012 gen** ön elemeden geçti. Bu ön eleme kayıpsızdır: P/LP missense sayısı toplam P/LP alel sayısını aşamaz.
   Bu 2.012 gen için kesin esearch sayımı yapıldı, **445 gen** >= 30 P/LP missense taşıyor.
3. **`pairs`**: 445 zengin gen için Ensembl Compara paralogları çekildi (23 gen için Ensembl tekrarlanan
   denemelere rağmen yanıt vermedi, bunlar dışarıda kaldı). Kardeş adayları önce özet tabloyla elendi
   (<= 8 P/LP alel), sonra 1.847 farklı kardeş için kesin esearch sayımı yapıldı. **3.678 aday çiftten 3.510'u**
   "kardeşte <= 2 P/LP missense" ölçütünü geçti.
4. **`metrics`**: 3.510 çiftin tamamı puanlandı (312'si puanlanamadı: 291'inde gnomAD'de MANE Select transkripti
   yok, 19'unda hiçbir P/LP kaydı MANE RefSeq'e oturmuyor, 2'sinde Ensembl proteini dönmedi). Kalan **3.198 çift**
   için metrikler hesaplandı.
5. **`rank`**: bileşik puan ve zayıflık bayrakları.

Hız notu: çift başına iş tek bir yerel hizalamadan ibaret olduğu için (1.000 aa için 0,02 s) bütün çiftleri
puanlamak, bir kısmını seçmekten daha ucuza geldi. Pahalı olan kısım ağ trafiğidir; gnomAD GraphQL takma adlarıyla
25 gen tek istekte, Ensembl `POST /sequence/id` ile 50 protein tek istekte çekildi.

### 2.3 Tuzak: `pathogenic[clinsig]` diye bir ClinVar alanı yok

Planlanan sorgu biçimi `GENE[gene] AND single_gene[prop] AND pathogenic[clinsig]` sessizce yanlış çalışıyor.
Entrez bunu `pathogenic[All Fields]` diye çeviriyor ve "not pathogenic", "Conflicting" gibi metinleri de yakalıyor.
KIF1A missense üzerinde ölçüldü (16 Eyl 2026):

| Sorgu | Sonuç |
|---|---|
| `pathogenic[clinsig] AND "missense variant"[molecular consequence]` | 456 kayıt (yanlış) |
| `("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties]) AND "missense variant"[molecular consequence]` | **135 kayıt** |

135, `docs/KIF1A_ClinVar_kurasyonu.md` içindeki kürasyon sayısıyla birebir aynı. Tarama boyunca yalnızca
`[Properties]` biçimi kullanıldı.

### 2.4 Doğrulama: boru hattı bilinen sayıları üretiyor mu?

Aynı kod yolu KIF1A/KIF13A üzerinde çalıştırıldığında projenin mevcut sayıları yeniden üretildi:

- hizalanan kalıntı çifti **1.405**, mevcut `data/interim/paralog_maps/KIF1A_to_KIF13A_pairwise_MANE.tsv`
  dosyasındaki kayıt sayısıyla aynı;
- motor domain (1-365) özdeşliği **0,622**, `docs/KIF1A_ClinVar_kurasyonu.md` içindeki %62 ile aynı;
- ClinVar P/LP missense **135**, hepsi MANE NM_001244008 üzerinde.

## 3. Ölçülen metrikler ve ne anlama geldikleri

| Sütun | Tanım |
|---|---|
| `identity_over_aligned` | Global hizalamada (BLOSUM62, açılış -10 / uzatma -0,5, uç boşlukları serbest) hizalanmış kalıntı çiftleri üzerinden özdeşlik. Projenin `paralog_map` aracıyla aynı parametreler. |
| `transfer_signal_residues` | Zengin genin benzersiz patojenik kalıntılarının kardeşte **özdeş** olma oranı. Aktarım sinyali budur. |
| `background_identity_in_window` | Aynı bölgedeki **bütün** hizalı kalıntıların özdeşliği. Dürüst karşılaştırma tabanı. |
| `transfer_lift` | Aktarım sinyali eksi yerel arka plan. Asıl bilgi taşıyan sayı: patojenik kalıntılar komşularından daha mı korunmuş? |
| `gap_fraction_at_plp_residues` | Patojenik kalıntıların kardeşte hizalı karşılığı olmayan oranı. Yüksekse bölge kardeşte yok demektir. |
| `window50_start/end`, `window50_frac_of_protein` | Patojenik kalıntıların yarısını içeren en dar pencere ve bunun protein uzunluğuna oranı. Kümelenme ölçütü. |
| `best_decile_occupancy` | Protein uzunluğunun %10'u genişliğindeki en yoğun pencerenin içindeki patojenik kalıntı oranı. |
| `sibling_loeuf`, `sibling_pli`, `sibling_mis_z` | gnomAD v4.1 kısıt ölçütleri. Kardeş seçilim altında değilse aktarımın dayanacağı zemin yok. |

Bileşik puan (ağırlıklar keyfidir ama açıktır, bütün bileşenler tabloda olduğu için okuyucu yeniden sıralayabilir):

```
feasibility = 0,30 * transfer + 0,25 * min(lift / 0,30, 1) + 0,15 * (1 - gap)
            + 0,20 * sibling_constraint + 0,10 * min(benzersiz_kalıntı / 100, 1)
```

## 4. Sıralama: zengin gen başına en iyi kardeş, ilk 12

`n+` = zengin gendeki ClinVar P/LP missense kaydı, `kal.` = benzersiz patojenik kalıntı,
`özd.` = global özdeşlik, `akt.` = aktarım sinyali, `arka` = pencere arka planı, `boşl.` = boşluk oranı.

| # | Çift | n+ | kal. | uzunluk | özd. | akt. | arka | lift | boşl. | kardeş LOEUF | kardeş mis_z | pencere | sıcak nokta domaini |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NSD1/NSD3 | 158 | 105 | 2696/1437 | 0,507 | 0,914 | 0,664 | 0,250 | 0,010 | 0,21 | 5,93 | 1968-2219 | SET 1942-2059, Post-SET 2066-2082 |
| 2 | COMP/THBS1 | 152 | 77 | 757/1170 | 0,556 | 0,909 | 0,655 | 0,254 | 0,013 | 0,31 | 3,77 | 397-515 | UniProt'ta "Disordered 298-503" |
| 3 | PRKCG/PRKCE | 48 | 32 | 697/737 | 0,532 | 0,938 | 0,667 | 0,271 | 0,000 | 0,28 | 5,49 | 101-139 | (C1B bölgesi, UniProt domain kaydı yok) |
| 4 | MTM1/MTMR4 | 65 | 52 | 603/1209 | 0,404 | 0,750 | 0,450 | 0,300 | 0,000 | 0,32 | 5,18 | 179-380 | Myotubularin fosfataz 163-538 |
| 5 | PCDH19/PCDH10 | 125 | 78 | 1148/1040 | 0,469 | 0,885 | 0,573 | 0,312 | 0,013 | 0,31 | 1,49 | 99-344 | Cadherin 1-3 (22-346) |
| 6 | DYSF/FER1L5 | 133 | 96 | 2119/2057 | 0,457 | 0,823 | 0,503 | 0,320 | 0,010 | 0,75 | 4,21 | 1010-1974 | C2 4-7 (1136-1943) |
| 7 | FLT4/EPHB2 | 34 | 27 | 1363/986 | 0,320 | 0,815 | 0,494 | 0,321 | 0,000 | 0,28 | 4,24 | 1035-1111 | Protein kinaz 845-1173 |
| 8 | RYR1/RYR3 | 279 | 199 | 5038/4870 | 0,685 | 0,920 | 0,678 | 0,242 | 0,005 | 0,40 | 2,59 | 3074-5006 | EF-hand 4074-4102 ve çevresi |
| 9 | KRT14/KRT18 | 47 | 29 | 472/430 | 0,502 | 0,862 | 0,432 | 0,430 | 0,000 | 0,36 | 3,04 | 268-415 | IF rod 115-426 |
| 10 | CSF1R/EPHB2 | 43 | 31 | 972/986 | 0,304 | 0,742 | 0,431 | 0,311 | 0,000 | 0,28 | 4,24 | 776-833 | Protein kinaz 582-910 |
| 11 | LMX1B/LHX5 | 43 | 27 | 402/402 | 0,399 | 0,852 | 0,531 | 0,321 | 0,000 | 0,34 | 2,01 | 223-271 | Disordered 176-229 |
| 12 | KCNJ11/KCNJ3 | 46 | 33 | 390/501 | 0,449 | 0,788 | 0,511 | 0,277 | 0,000 | 0,23 | 6,08 | 164-301 | (UniProt domain kaydı yok) |

Tam liste: `paralog_scan_ranked_20260916.tsv` (3.198 satır, 28 sütun).

## 5. Zayıf vakalar ve neden başarısız oldukları

Bayraklar 3.198 çift üzerinde:

| Bayrak | Kaç çift | Ölçüt | Anlamı |
|---|---|---|---|
| `no_conservation_lift` | 1.737 (%54) | lift <= 0,05 | Patojenik kalıntılar komşularından daha korunmuş değil. Paralog kanıtı bir şey katmaz. |
| `sibling_unconstrained` | 1.722 (%54) | kısıt puanı < 0,35 | Kardeş gen popülasyonda seçilim altında değil. Aktarılan eşiğin dayanacağı zemin yok. |
| `region_absent_in_sibling` | 927 (%29) | boşluk > 0,25 | Patojenik kalıntıların dörtte birinden fazlası kardeşte hizalanamıyor. |
| `low_identity` | 837 (%26) | pencere arka planı < 0,30 | Sıcak nokta bölgesi zaten tanınmayacak kadar uzaklaşmış. |
| `thin_labels` | 157 (%5) | benzersiz kalıntı < 25 | Kalıntı gruplu çapraz doğrulama için etiket yetersiz. |

Adıyla sanıyla başarısız örnekler:

- **LDLR/LRP12, LDLR/LRP10, LDLR/LRP3**: LDLR 746 P/LP missense ile listenin en zengin genlerinden biri, ama
  patojenik kalıntıların %28 ila %48'i kardeşlerde hizalanmıyor ve lift negatif (-0,040, -0,093, -0,093).
  LDLR'nin LDL-reseptör tekrarları büyük ailelerde paylaşılıyor, ama patojenik konumlar tekrar sayısı ve
  aralayıcılar farklı olduğu için karşılık bulmuyor. **Zengin olmak yeterli değil.**
- **COL4A5/COL26A1, COL3A1/COL26A1, COL3A1/EMID1**: kolajenlerde global özdeşlik yanıltıcı biçimde yüksek
  görünüyor (0,43 ila 0,50), çünkü Gly-X-Y tekrarı her yerde eşleşiyor; ama patojenik kalıntıların %76 ila %83'ü
  kardeşte hizalanamıyor ve pencere arka planı hesaplanamıyor. **Tekrar bölgeleri özdeşlik ölçütünü şişiriyor.**
- **PAX2/PRRX2, PAX2/DRGX, PAX2/SEBOX**: PAX2'nin 41 P/LP missense'i paired box içinde; bu kardeşlerde yalnızca
  homeodomain var, paired box yok. Boşluk oranı 0,92 ila 1,00, aktarım sinyali **0,000**. Sıralamanın en dibi.
  **Patojenik varyantlar kardeşte bulunmayan bir bölgede.**
- **COL6A3/VIT, COL6A3/VWA2**: boşluk oranı 1,00, yani hiçbir patojenik kalıntının karşılığı yok; üstelik
  kardeş LOEUF 1,50 ve 1,26, kısıtlanmamış. Dört bayrağın dördü birden.

Ayrıca dikkat: **UniProt'un "Disordered" etiketi bir uyarı işaretidir.** COMP/THBS1 (2. sıra) ve COL1A2/COL16A1
gibi çiftlerde sıcak nokta düzensiz bölgeye düşüyor; orada kalıntı düzeyinde aktarım varsayımı yapısal olarak
zayıftır, puan yüksek olsa bile.

## 6. KIF1A/KIF13A bu sıralamanın neresinde

| Ölçüt | Değer |
|---|---|
| Genel sıra | **152 / 3.198** |
| Zengin gen başına en iyi kardeş listesinde | **71 / 337** |
| Bileşik puan | 0,7305 (1. sıra 0,931) |
| ClinVar P/LP missense | 135, benzersiz kalıntı 86 |
| Protein uzunlukları | 1791 / 1805 |
| Global özdeşlik | 0,425 |
| Aktarım sinyali (kalıntı) | **0,826** (medyan 0,412) |
| Aktarım sinyali (varyant) | 0,867 |
| Pencere arka planı (210-344) | 0,698 |
| Lift | 0,128 (medyan 0,048) |
| Boşluk oranı | 0,000 |
| Sıcak nokta penceresi | 210-344, proteinin %7,5'i, en yoğun ondalıkta kalıntıların %60,5'i |
| Sıcak nokta domaini | UniProt "Kinesin motor 5-354" |
| Zayıflık bayrağı | yok |
| KIF13A kısıtı | LOEUF 0,449, pLI 0,9998, missense Z 2,79 |

Okuma: KIF1A/KIF13A **aktarım sinyali bakımından üst dilimde** (0,826, medyanın iki katı), **kapsama bakımından
kusursuz** (hiç boşluk yok), **kümelenme bakımından örnek** (patojenik kalıntıların yarısı proteinin %7,5'ine
sığıyor ve bu tam olarak motor domaindir). Puanı 1. sıraya taşımayan tek şey **lift'in 0,128'de kalması**: KIF13A'nın
motor domaini zaten %62 ila %70 özdeş olduğu için patojenik kalıntıların ekstra korunumu görece küçük görünüyor.
Bu bir kusur değil, ölçütün bir özelliğidir: yüksek arka plan hem aktarımı kolaylaştırır hem de lift'i bastırır.
`docs/kalibrasyon_olcumleri.md` bunu bağımsız olarak doğruluyor: "KIF13A'da aynı kalıntı" özelliği tek başına
motor içi AUROC 0,744 veriyor, yani lift küçük görünse de kullanılabilir bilgi taşıyor.

KIF1A'nın 26 boş kardeşi içinde KIF13A **birinci sırada** (152.), onu CENPE (410.), KIFC3 (619.) ve KIFC2 (863.)
izliyor. Yani projenin kardeş seçimi, bu bağımsız taramaya göre doğru seçimdir.

## 7. İki aday çift ve somut plan

Tam boru hattı hiçbir çiftte çalıştırılmadı. Aşağıdaki iki çift, hem ölçütleri hem de projenin nörolojik
odağını karşıladığı için seçildi.

### 7.1 NSD1 -> NSD3 (sıralamada 1.)

- **Klinik zemin**: NSD1, Sotos sendromu (aşırı büyüme ve zihinsel yetersizlik). 158 P/LP missense, 105 benzersiz
  kalıntı, ClinVar'da toplam 2.842 alel.
- **Kardeş**: NSD3 (WHSC1L1), ClinVar'da 247 alel ve **1 P/LP alel**, 0 P/LP missense. gnomAD LOEUF **0,21**,
  pLI 1,0, missense Z **5,93**: taramadaki en kısıtlı kardeş.
- **Aktarım**: patojenik kalıntıların %91,4'ü NSD3'te özdeş, pencere arka planı 0,664, lift 0,250, boşluk 0,010.
- **Kümelenme**: pencere 1968-2219, proteinin %9,3'ü, UniProt SET (1942-2059) ve Post-SET (2066-2082) domainleri.
  KIF1A'nın motor domaini ile aynı biçimde tek bir katalitik domain.
- **Negatif seti kurulabilir mi?** gnomAD v4.1 (16 Eyl 2026, ENST00000439151): NSD1 üzerinde 2.789 missense,
  sıcak nokta penceresi içinde AC >= 3 olan **31**, AC >= 2 olan 43 varyant. KIF1A'nın aynı dar penceresinde
  AC >= 3 sayısı 20'dir (proje motor domainin tamamını kullanarak 61'e çıkmıştı), yani NSD1 en az KIF1A kadar
  elverişli. Domain sınırı pencere yerine SET+Post-SET olarak alınırsa sayı daha da artar.
- **Aday listesi büyüklüğü**: NSD3 üzerinde gnomAD'de 1.336 missense; doygunluk matrisi 1437 x 19 = 27.303
  (KIF13A'da 34.295).

### 7.2 PCDH19 -> PCDH10 (sıralamada 5.)

- **Klinik zemin**: PCDH19, X'e bağlı gelişimsel epileptik ensefalopati (DEE9). 125 P/LP missense, 78 benzersiz
  kalıntı, ClinVar'da 1.772 alel.
- **Kardeş**: PCDH10, ClinVar'da 165 alel ve **hiç P/LP aleli yok**. LOEUF 0,311, pLI 1,0.
- **Aktarım**: patojenik kalıntıların %88,5'i özdeş, arka plan 0,573, lift **0,312** (ilk beşin en yükseği),
  boşluk 0,013.
- **Kümelenme**: pencere 99-344, proteinin %21,4'ü, UniProt Cadherin 1, 2 ve 3 tekrarları (22-346). KIF1A'ya göre
  daha geniş ama yine de tanımlı ardışık domainler.
- **Negatif seti**: gnomAD v4.1 (ENST00000373034): PCDH19 üzerinde 1.105 missense, pencere içinde AC >= 3 olan
  **41**, AC >= 2 olan 73. KIF1A'dan daha rahat.
- **Aday listesi**: PCDH10 üzerinde 1.524 gnomAD missense; doygunluk matrisi 1040 x 19 = 19.760, taramadaki en
  ucuz hesaplardan biri.
- **Açık risk**: PCDH10'un missense Z değeri yalnızca **1,49**. LoF'a karşı kısıtlı ama missense'e karşı kısıtı
  zayıf. "Kardeş missense seçilimi altında" varsayımı burada NSD3'teki kadar sağlam değil ve raporlanmalı.

### 7.3 Seçilmeyenler ve nedeni

- **RYR1/RYR3** (8.): ölçütleri en güçlü çiftlerden biri (özdeşlik 0,685, aktarım 0,920, 279 P/LP missense) ama
  proteinler 5038 ve 4870 aa. Doygunluk taraması 4870 x 19 = 92.530 varyant, AlphaFold3 üzerinde tetramerik
  kanal pratikte çalıştırılamaz ve sıcak nokta penceresi proteinin %38'ine yayılıyor (RYR1'in bilinen üç ayrı
  sıcak bölgesi, tek domain değil). Hesap maliyeti kazancı haklı çıkarmıyor.
- **COMP/THBS1** (2.): sıcak nokta UniProt'ta düzensiz bölge olarak işaretli ve COMP'un mekanizması baskın-negatif;
  "patojenik = işlev kaybı" varsayımı KIF1A'dakinden daha da kırılgan.
- **PRKCG/PRKCE** (3.): en ucuz hesap (697/737 aa) ve en yüksek aktarım sinyali (0,938), SCA14 ile nörolojik.
  Yedek aday olarak tutulmalı; tek eksiği 32 benzersiz kalıntı, yani kalıntı gruplu çapraz doğrulamada kat başına
  çok az örnek.

### 7.4 Çift başına gereken iş ve süre

Aşağıdaki tahminler KIF1A üzerinde ölçülen gerçek sürelere dayanıyor
(`docs/KIF1A_ClinVar_kurasyonu.md`, `docs/kalibrasyon_olcumleri.md`).

| Adım | KIF1A'daki gerçek maliyet | NSD1 ve PCDH19 için tahmin |
|---|---|---|
| ClinVar kürasyonu (VCV XML, tier atama, döngüsellik ayıklaması) | 3.372 tek gen kaydı, betik çalışması saatler; asıl iş elle inceleme | NSD1 2.842 alel, PCDH19 1.772 alel: her biri için 1 ila 2 gün elle inceleme |
| gnomAD domain içi negatif seti | tek transkript sorgusu, dakikalar | aynı, dakikalar |
| Tek VEP çalışması (AlphaMissense, CADD, REVEL, MetaRNN, ESM1b, PrimateAI, SpliceAI) | 171 varyant tek koşu | 190 ila 210 varyant, tek koşu |
| Kalıntı gruplu bootstrap metrikleri | dakikalar | dakikalar |
| Paralog haritası ve korunum referansı | dakikalar | dakikalar |
| Evo 2 40B (isteğe bağlı) | 171 varyant = 342 NIM çağrısı, 3,3 s/varyant | benzer, kredi harcar |
| Kardeş için aday listesi ve şekiller | mevcut betikler | mevcut betikler |

Gerçekçi takvim: **çift başına yaklaşık 1 hafta**, bunun büyük kısmı ClinVar kürasyonunun elle kısmı.
Hesaplama kısmı (VEP, hizalama, metrikler) çift başına birkaç saat. İki çift için toplam 2 ila 3 hafta, eğer
kürasyon adımları KIF1A betiklerinden genelleştirilirse (gen adı, MANE numaralandırması ve domain sınırları
parametre haline getirilerek) 1,5 ila 2 haftaya iner.

### 7.5 Çok çiftli çalışmanın asıl tasarım sorunu

Üç gen (KIF1A, NSD1, PCDH19) bir yöntem makalesi için hâlâ küçük bir n'dir. Ölçeklenebilir tasarım şudur:
ucuz metrikleri **tahmin edici** olarak kullan, tam boru hattını **doğrulama** olarak kullan. Yani:

1. 3.198 çift için ucuz metrikler zaten var (bu tarama).
2. Aktarım kazancının yüksek, orta ve düşük olduğu dilimlerden **rastgele 12 ila 15 çift** seç (yalnız en iyileri
   değil; yalnız en iyileri seçmek sonucu tanım gereği olumlu çıkarır).
3. Her birinde tam boru hattını çalıştır ve tek bir çıktıyı ölç: paralog korunum özelliğinin model üstüne
   **kattığı** AUROC farkı (KIF1A'da referans 0,744, AlphaMissense 0,925).
4. Ucuz `transfer_lift` ile bu gerçek kazanç arasındaki korelasyonu raporla. Yöntem iddiası budur:
   "boru hattı çalıştırılmadan hangi paralog çiftinde aktarımın işleyeceği önceden kestirilebilir."

Bu tasarımda negatif sonuç da yayımlanabilir bir sonuçtur ve tarama zaten negatif tarafın büyük olduğunu
gösteriyor.

## 8. Bu taramanın ölçmediği şeyler

1. **Lift bir AUROC değildir.** Patojenik kalıntıların daha korunmuş olması, korunumun modelden bağımsız katkı
   yapacağını garanti etmez. KIF1A'da bu bağ ölçüldü (0,744), diğer çiftlerde ölçülmedi.
2. **Tespit yanlılığı.** Patojenik varyantlar iyi çalışılmış, korunmuş katalitik çekirdeklerde birikir. Lift'in bir
   kısmı gerçek biyoloji değil, ClinVar'ın nereye baktığıdır. Doğru kontrol yerel pencere arka planıdır ve bu
   kullanıldı, ama ascertainment tamamen giderilemez.
3. **Domain karışıklığı her çiftte yeniden çıkar.** KIF1A'nın ana dersi, ClinVar B/LB negatiflerinin hiçbirinin
   motor domainde olmamasıydı ("motor domainde mi" tek başına AUROC 0,995). Bu tarama her çift için
   `window50` ve `best_decile_occupancy` veriyor ama negatif setin domain dağılımını **ölçmüyor**. Her yeni
   çiftte gnomAD domain içi negatifleri kurulmadan hiçbir AUROC raporlanmamalı.
4. **Mekanizma tek yönlü varsayılıyor.** KIF1A'da baskın-negatif ve aşırı aktive edici varyantlar var; "patojenik"
   tek bir yön değil. Bu tarama mekanizmayı hiç görmüyor.
5. **Hizalama tek yöntemle yapıldı.** Global ikili hizalama, projenin `paralog_map` aracıyla aynı parametrelerle.
   Kesin kalıntı eşlemesi için WP1'in çoklu paralog MSA'sı ve yapısal hizalama gerekir; özellikle tekrar içeren
   ailelerde (kolajenler, keratinler, LDLR tekrarları) ikili hizalama yanılabilir.
6. **23 zengin gen taramaya hiç giremedi**, Ensembl homoloji ucu tekrarlanan denemelerde yanıt vermedi
   (SPAST dahil). 312 çift de MANE veya konum eşleşmesi eksikliğinden puanlanamadı.
7. **Kardeşin "klinik olarak boş" olması gen gerçekten masum olduğu için değil, kimse bakmadığı için olabilir.**
   Tarama bu ikisini ayırt edemez; zaten projenin varsayımı da ikincisidir.

## 9. Yeniden üretim

```powershell
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\paralog_pairs_scan.py --stage sources
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\paralog_pairs_scan.py --stage rich
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\paralog_pairs_scan.py --stage pairs
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\paralog_pairs_scan.py --stage metrics
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\paralog_pairs_scan.py --stage rank
```

Her aşama `data/raw/paralog_scan/` altına önbellek yazar, ikinci çalıştırmada ağa çıkmaz. `--refresh` önbelleği
yok sayar. Ölçülen süreler (16 Eyl 2026, ev bağlantısı): `rich` 303 s, `pairs` 12.892 s (Ensembl homoloji ucu
sık 503 döndürdü ve iki geçiş gerekti), `metrics` 1.925 s, `rank` 1 s'den az.

## 10. Kaynaklar

- ClinVar `gene_specific_summary.txt`, sürüm 13 Eyl 2026, indirme 16 Eyl 2026.
- NCBI E-utilities esearch/esummary, ClinVar veritabanı, sorgu 16 Eyl 2026.
- HGNC complete set, indirme 16 Eyl 2026.
- Ensembl REST (GRCh38), homology ve sequence uçları, sorgu 16 Eyl 2026.
- gnomAD v4.1 (`gnomad_r4`) GraphQL, sorgu 16 Eyl 2026.
- UniProt REST (reviewed, organism 9606), sorgu 16 Eyl 2026.
- Proje içi referanslar: `docs/KIF1A_ClinVar_kurasyonu.md`, `docs/kalibrasyon_olcumleri.md`,
  `data/interim/paralog_maps/KIF1A_to_KIF13A_pairwise_MANE.tsv`.
