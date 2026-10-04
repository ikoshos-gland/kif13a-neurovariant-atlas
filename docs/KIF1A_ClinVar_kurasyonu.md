# KIF1A ClinVar kürasyonu ve kalibrasyon setleri

Tarih: 15 Eyl 2026. Kaynak: ClinVar (NCBI E-utilities, 3.372 tek-gen kaydı), VCV XML gönderim detayları,
Ensembl VEP (MANE NM_001244008.2), gnomAD v4.1, Europe PMC tam metinleri.
Çıktı: `data/interim/clinvar_curation/KIF1A_calibration_sets_20260915.tsv`.

## 1. Özet

ClinVar'daki KIF1A etiketleri, plandaki haliyle (P/LP ve B/LB karşılaştırılarak) kullanılırsa **yanıltıcı**
bir kalibrasyon verir. Pozitiflerin tamamına yakını motor domaininde, negatiflerin **hiçbiri** motor
domaininde değil. Sadece "motor domaininde mi?" diye soran bir kural AUROC 0,995 alıyor. Bu yüzden:

- pozitifler ClinVar'dan (sıkı filtreyle),
- motor domain negatifleri gnomAD'den (popülasyonda tekrar görülen varyantlar),
- ClinVar benign seti yalnızca motor dışı ikincil değerlendirme için

kullanılmalı ve ana metrik **motor domain içi AUROC** olmalı.

## 2. Profil

| Sonuç | P/LP | B/LB | VUS | Çelişkili |
|---|---|---|---|---|
| missense | 135 | 90 | 980 | 213 |
| kanonik splice | 31 | 2 | 5 | 0 |
| nonsense + frameshift | 33 | 0 | 7 | 0 |
| intronik | 0 | 793 | 106 | 48 |
| sinonim | 0 | 646 | 25 | 78 |

İnceleme düzeyi: 2.409 tek gönderen (1 yıldız), 565 çok gönderen (2 yıldız), 346 çelişkili, 23 kriter yok.
P/LP kayıtlarının hastalık alanı çoğunlukla toplu üçlü (SPG30 + ID AD9 + HSN2C): **mekanizma bu alandan
çıkarılamaz**.

## 3. Bulgular

1. **Domain karışıklığı.** Çekirdek 111 pozitifin 110'u motor domaininde (1-365). ClinVar'ın 88 B/LB missense
   varyantının 0'ı motor domaininde (62 sap, 14 PH, 5 CC bölgesi, 2 C-uç, 5 eşlenemedi). Tüm protein AUROC'u: AlphaMissense 0,959,
   CADD 0,887, REVEL 0,993, "motor domaininde mi" 0,995.
2. **Motor içi gerçek ayrım.** Pozitifler (motor, n=109) ile gnomAD motor missense kontrolleri karşılaştırıldı:

   | Negatif seti | n | AlphaMissense | CADD | REVEL |
   |---|---|---|---|---|
   | gnomAD AC ≥ 1 | 230 | 0,896 | (karışık kaynak) | 0,910 |
   | gnomAD AC ≥ 2 | 93 | 0,914 | (karışık kaynak) | 0,921 |
   | gnomAD AC ≥ 3 | 61 | 0,925 | (karışık kaynak) | 0,920 |
   | gnomAD AC ≥ 5 | 42 | 0,945 | (karışık kaynak) | 0,938 |

   > Bu tablodaki CADD değerleri geçersiz: pozitiflerin skoru VEP'ten, negatiflerin skoru gnomAD
   > anotasyonundan geliyordu (iki farklı sürüm). Tek kaynaktan ölçülen nihai değerler ve güven
   > aralıkları: `docs/kalibrasyon_olcumleri.md`. Motor içi CADD 0,781, AlphaMissense 0,925.

3. **Döngüsellik (circularity).** Missense P/LP yorumlarının 75/135'i in silico kanıttan (PP3, REVEL, CADD...)
   söz ediyor. Tek gönderenli, in silico kanıta dayanan ve de novo/fonksiyonel kanıtı olmayan 18 LP ayrı tutuldu.
   REVEL ve MetaRNN klinik etiketlerle eğitildi; AlphaMissense popülasyon frekansıyla eğitildi (gnomAD tabanlı
   negatiflerde avantajlı). CADD etiketsiz eğitildiği için en az döngüsel referans.
4. **ClinVar B/LB gerekçesi belgelenmemiş.** 90 negatifin 85'inin yorumunda ne sıklık ne in silico gerekçesi
   geçiyor; 52'si gnomAD'de AF < 1e-5
   (yani sıklık nedeniyle değil). Motor dışında olmaları büyük olasılıkla belirleyici.
5. **Paralog korunumu ayrıştırıcı.** Çekirdek pozitiflerin 98/111'i KIF13A'da özdeş kalıntıya denk geliyor
   (70 kalıntının 60'ı); motor gnomAD kontrollerinde 24/61. Yalnız bu özellik motor içi AUROC 0,744 verir.
6. **Zor negatif: p.Ile119Thr.** gnomAD'de 1.165 allel (AF 0,0022), homozigot yok, ClinVar'da yok, AlphaMissense
   0,992. Motor domainde AlphaMissense'in aşırı çağırdığı yerlerin göstergesi.
7. **Mekanizma heterojen.** Boyle 2021 hareketlilik kaybı sınıfları, Anazawa 2022 baskın-negatif etki, Chiba 2019
   SPG varyantlarında aşırı aktivasyon tanımlıyor. V8M (aşırı aktive edici) ClinVar'da yok. KIF13A'ya aktarımda
   "patojenik" tek yönlü bir etki varsayılmamalı.
8. **Literatür kurtarması.** A255V yalnız OMIM kaydı (0 yıldız) olduğu için dışlanıyordu; iki fonksiyonel makalede
   çalışıldığı için çekirdeğe alındı. Fonksiyonel literatürde geçen çekirdek pozitif: 14.
9. MaveDB'de KIF1A derin mutasyon taraması yok.

## 4. Kullanılacak setler

| Set | n | Etiket | Kullanım |
|---|---|---|---|
| `missense_positive_core` | 111 | 1 | Ana pozitifler; 70 kalıntı, CV grubu kalıntı |
| `missense_negative_motor_population` | 61 | 0 | Ana negatif (motor içi), gnomAD AC ≥ 3 |
| `missense_negative_motor_population_ac2` | 32 | 0 | Duyarlılık |
| `missense_negative_clinvar_outside_motor` | 88 | 0 | Yalnız motor dışı, domain-tabakalı raporla |
| `splice_positive` | 26 | 1 | Evo 2 / SpliceAI kalibrasyonu (negatif seti ayrıca kurulacak) |
| `heldout_pp3_dependent_positive` | 18 | 1 | Döngüsellik testi; eğitime ve ana metriğe girmez |
| `sensitivity_conflict_P_vs_VUS` | 28 | yok | Duyarlılık |
| `reference_lof` | 33 | yok | HSN2C / haployetmezlik referansı, missense değil |

Dışlananlar: 183 çelişkili B-VUS, P ve B'yi kapsayan 2 çelişki, 4 kriter belirtilmemiş kayıt (A255V kurtarıldı),
8 SNV olmayan kayıt, 7 protein değişimi olmayan kayıt, 2 BP4'e dayanan benign. Referans olarak ayrı duranlar:
5 kodlamayan B/LB, 5 çerçeve içi P/LP.

## 5. Değerlendirme protokolü (öneri)

0. Ölçüm sonuçları: `docs/kalibrasyon_olcumleri.md` (15 Eyl 2026).
1. **Ana metrik:** motor domain içi AUROC, `missense_positive_core` (motor) ve `missense_negative_motor_population`.
2. Her modelle birlikte iki referans raporlanır: "motor domaininde mi" (tüm protein) ve "KIF13A'da özdeş kalıntı" (motor içi).
3. **Çapraz doğrulama kalıntı gruplu** (`cv_group`): aynı kalıntının farklı değişimleri aynı katta.
4. Konsensüse REVEL, MetaRNN, ClinPred gibi klinik etiketle eğitilmiş skorlar **girmez**. Girerse yalnızca ilk
   ClinVar tarihi eğitim verisinden sonraki varyantlarda değerlendirilir.
5. Duyarlılık: fonksiyonel literatür alt kümesi (14), AC ≥ 2 negatifleri, PP3'e bağlı ayrılmış set.
6. Plandaki başarı ölçütü güncellenmeli: "KIF1A AUROC ≥ 0,85" motor içi ve popülasyon kontrolleriyle ölçülmeli.

## 6. Açık işler

- Chiba 2019 tam metni Europe PMC'de yok; aşırı aktive edici varyant listesi elle eklenmeli.
- Fonksiyonel sınıflar (MT bağlanma, hız, işlemsellik, aşırı aktivasyon) Boyle 2021 ekinden tabloya işlenmeli.
- Splice negatifleri: gnomAD'de yaygın, AF eşleşmiş intronik varyantlar (Evo 2 panelindeki yaklaşım).
- Regex tabanlı literatür taraması elle kontrol edilmeli (şekil etiketleri yanlış eşleşebilir).

## 7. Yeniden üretim

```powershell
mcp\atlas_server\.venv\Scripts\python.exe scripts\pull_wp1_core.py --only clinvar,gnomad,alphamissense
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\kif1a_clinvar_curation.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\kif1a_curation_extras.py
mcp\atlas_server\.venv\Scripts\python.exe scripts\kif1a_build_calibration_sets.py
```

## 8. Kaynaklar (Europe PMC'de doğrulandı)

- Boyle L ve ark. HGG Adv 2021. PMID 33880452
- Chiba K ve ark. PNAS 2019. PMID 31455732
- Anazawa Y ve ark. PNAS 2022. PMID 35917346
- Shatarupa A ve ark. Nat Commun 2026. PMID 41912506
- Kaur S ve ark. Hum Mutat 2020. PMID 32652677
- Nicita F ve ark. J Med Genet 2021. PMID 32737135
